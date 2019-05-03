##### This code, is all from gramtools discover ######
import bisect

## Class that allows mapping between a vcf record in one REF coordinate to a vcf record in another REF coordinate system.
#  Call the two references coordinates ref 2 and ref 1.
#  A `_Region` either marks a region in ref 2 within a variant region in ref 1, or a region in ref 2 within an invariant region in ref 1.
class _Region:
    def __init__(self, base_POS, inf_POS, length, vcf_record_REF=None, vcf_record_ALT=None):
        ## Start coordinate in base reference.
        self.base_POS = base_POS
        ## Start coordinate in inferred reference.
        self.inf_POS = inf_POS
        ## Holds the length of the region
        self.length = length

        ## If in a variant site, the only parts of the vcf_record that we need to track are the REF and ALT sequences.
        self.vcf_record_REF = vcf_record_REF

        self.vcf_record_ALT = vcf_record_ALT

    @property
    def is_site(self):
        return self.vcf_record_REF is not None

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return str(self.__dict__)

    def __lt__(self, other):
        if isinstance(other, _Region):
            return self.inf_POS < other.inf_POS
        else:
            return self.inf_POS < other

    def __len__(self):
        return self.end - self.start + 1

## Marks the personalised reference with `Regions` flagging a relationship with the base reference.
#  There are two types of `Regions`: those inside variant sites in the prg, and those outside.
# Important requirements for vcf:
# * the records are sorted by POS for a given ref ID, and all records referring to the same ref ID are contiguous.
# * the chrom_sizes are ordered in same order as the ref IDs in the records.
# The first point is explicitly checked using assertions.
def _flag_personalised_reference_regions(base_records, chrom_sizes):
    all_personalised_ref_regions = {}  # Maps each CHROM to a list of _Regions
    chrom_positions = {}  # Maps each CHROM to a base_ref and an inferred_ref position
    num_chroms = -1
    prev_chrom_key = None
    prev_record = None

    for vcf_record in base_records:
        chrom_key = vcf_record.CHROM
        if chrom_key not in all_personalised_ref_regions:
            all_personalised_ref_regions[chrom_key] = []
            chrom_positions[chrom_key] = [1, 1]

            # We enter a new ref ID; let's make sure we capture the end of the previous ref ID
            if num_chroms >= 0:
                base_ref_pos, inf_ref_pos = chrom_positions[prev_chrom_key][0], chrom_positions[prev_chrom_key][1]
                if base_ref_pos <= chrom_size:
                    regions_list = all_personalised_ref_regions[prev_chrom_key]
                    non_var_region = _Region(base_POS=base_ref_pos, inf_POS=inf_ref_pos,
                                             length=chrom_size - base_ref_pos + 1)
                    regions_list.append(non_var_region)

            num_chroms += 1
            chrom_size = chrom_sizes[num_chroms]

        else:
            # Enforce ref ID contiguity
            assert chrom_key == prev_chrom_key, "Ref IDs not contiguous: {} and {} interspersed"\
                                                          .format(chrom_key, prev_chrom_key)
            # Enforce position sortedness
            assert vcf_record.POS > prev_record.POS, "Records not in increasing POS order: {} and {}"\
                                                            .format(prev_record, vcf_record)
        regions_list = all_personalised_ref_regions[chrom_key]
        base_ref_pos, inf_ref_pos = chrom_positions[chrom_key][0], chrom_positions[chrom_key][1]

        # Build non-variant region
        if vcf_record.POS > base_ref_pos:
            non_var_region = _Region(base_POS=base_ref_pos, inf_POS=inf_ref_pos,
                                     length=vcf_record.POS - base_ref_pos)
            regions_list.append(non_var_region)
            base_ref_pos += non_var_region.length
            inf_ref_pos += non_var_region.length

        # Build variant region
        # Note the following 'guarantee' from `infer`: the first element of GT is the single most likely haploid genotype call.
        # We need to take that, as that is what gets used to produce the inferred reference fasta- which is the reference used for variant calling.
        picked_alleles = vcf_record.samples[0].gt_alleles

        if set(picked_alleles) == {None}:
            picked_allele = 0  # GT of './.'
        else:
            picked_allele = int(picked_alleles[0])  # Get the first one. Can be 0 (REF) !

        # Only make a var region if the inference procedure did not pick the REF.
        if picked_allele != 0:
            var_region = _Region(base_POS=base_ref_pos, inf_POS=inf_ref_pos,
                                 length=len(vcf_record.ALT[picked_allele - 1]), vcf_record_REF=vcf_record.REF,
                                 vcf_record_ALT=str(vcf_record.ALT[picked_allele - 1]))

            regions_list.append(var_region)
            base_ref_pos += len(vcf_record.REF)
            inf_ref_pos += var_region.length

        chrom_positions[chrom_key] = [base_ref_pos, inf_ref_pos]
        prev_chrom_key = chrom_key
        prev_record = vcf_record

    # End game: deal with the last non-var region if there is one.
    if len(all_personalised_ref_regions) == 0:
        raise ValueError("No records in provided vcf.")

    if base_ref_pos <= chrom_size:
        non_var_region = _Region(base_POS=base_ref_pos, inf_POS=inf_ref_pos,
                                 length=chrom_size - base_ref_pos + 1)
        regions_list.append(non_var_region)

    return all_personalised_ref_regions

## Return the marked `_Region` containing the position referred to by `vcf_record`.
def _find_start_region_index(vcf_record, marked_regions):
    record_start_POS = vcf_record.POS

    # Binary search on a sorted array. The '<' comparator is overloaded in `_Region` definition to allow the search.
    # The produced index is:
    # * If a region with start POS == record_start_POS is found, returns that index (specifically, the first occurrence)
    # * If such a region is not found, returns the index of the region with start POS **strictly** greater than record_start_POS
    # IN THAT CASE, we need to return the region index one smaller; so that the region's interval includes record_start_POS.
    region_index = bisect.bisect_left(marked_regions, record_start_POS)

    #  Case: `record_start_index` is larger than any _Region start.
    if region_index > len(marked_regions) - 1:
        region_index = len(marked_regions) - 1

    selected_region = marked_regions[region_index]

    # Case: record_start_POS was not found exactly.
    if selected_region.inf_POS > record_start_POS:
        region_index -= 1

    return region_index


def rebase_mutant_record(old_vcf_record, regions_list, bases):
    """
    This is a simpler version of the rebasing in `discover`
    Note the following strong assumption: the length of the mutated portion of the inferred REF is 1.
    To extend this, use the full rebasing routine in `discover`- where we deal with potentially overlapping multiple regions.
    """
    index_of_region_hit = _find_start_region_index(old_vcf_record, regions_list)

    region_hit = regions_list[index_of_region_hit]


    if region_hit.is_site:
        rebased_POS = region_hit.base_POS
        rebased_REF = region_hit.vcf_record_REF
        rebased_ALT = ""
        # We also straight away pre-pend any preceding variation relative to the base REF
        if old_vcf_record.POS > region_hit.inf_POS:
            record_inset = old_vcf_record.POS - region_hit.inf_POS
            rebased_ALT += region_hit.vcf_record_ALT[:record_inset] + rebased_ALT

        rebased_ALT += bases

        consumed_bases = old_vcf_record.POS - region_hit.inf_POS + len(bases)
        if consumed_bases < region_hit.length:
            rebased_ALT+= region_hit.vcf_record_ALT[consumed_bases:]

    else:
        rebased_POS = region_hit.base_POS + (old_vcf_record.POS - region_hit.inf_POS)
        rebased_REF = old_vcf_record.REF
        rebased_ALT = bases
        
    rebased_record = _modify_vcf_record(old_vcf_record, POS = rebased_POS, REF = rebased_REF, ALT = [rebased_ALT])

    return rebased_record


## Take all vcf_record class attributes, modify those of interest, and make a new vcf record.
# For valid `new_attributes`, see https://pyvcf.readthedocs.io/en/latest/API.html#vcf-model-record
def _modify_vcf_record(vcf_record, **new_attributes):
    for attribute, value in new_attributes.items():
        setattr(vcf_record, attribute, value)

    return vcf_record
