## @file Make new fasta from fasta + vcf.

import vcf
import random
import sys
import os
import collections

DNA = {'A', 'C', 'G', 'T', 'N'}
REF_READ_CHUNK_SIZE = 1000000

## An interface to the module for making new fasta (recombinant) from fasta + vcf.
def make_recombinant(reference_fpath, vcf_fpath, output_fpath, seed = None, offset = 1):

    CallData=collections.namedtuple("CallData",["GT"])
    picked_alleles = []

    random.seed(seed)

    fhandle = open(vcf_fpath)
    vcf_reader = vcf.Reader(fhandle)
    description = "Recombinant built from ref and vcf"

    new_records = []
    alt_records = []

    for vcf_record in vcf_reader:

        allele_index = random.randrange(len(vcf_record.ALT) + 1) # Add 1 to choose between REF and all ALTs
        picked_alleles.append(allele_index)

        new_records.append(vcf_record)

        # Only write the record if we picked it
        if allele_index > 0:
            alt_records.append(vcf_record)
            vcf_record.samples = vcf_record.samples[:1]
            print(vcf_record.samples)
            vcf_record.samples[0] = vcf.model._Call(site = vcf_record, sample = "Picked", data=CallData(GT=str(allele_index)+"/"+str(allele_index))) # Pick the allele
            vcf_record.FORMAT = None; vcf_record.INFO = None
            print(vcf_record.samples)

    if sum(picked_alleles) == 0:  # We did not pick any ALTs; let's pick one
        switch_pos = random.randrange(len(picked_alleles))
        new_records[switch_pos].samples[0].gt_alleles = [1,1]
        alt_records.append(vcf_record)

    # We only need to pass records with ALT choices; records with a REF choice (GT of 0) work, but if absent the procedure will pick REF anyway.
    make_new_ref_using_vcf(reference_fpath, iter(alt_records), output_fpath)

    fhandle.close()

    return alt_records


def make_new_ref_using_vcf(original_fasta, vcf_reader, out_fasta):
    """
    Assumption: the vcf records are sorted in same order as the original_fasta
    :param original_fasta:
    :param vcf_reader:
    :param out_fasta:
    :return:
    """
    in_parse = _refParser(original_fasta)
    out_parse = _FastaWriter(out_fasta)

    vcf_record = next(vcf_reader, None)
    if vcf_record is None:
        raise LookupError("vcf file {} does not contain any records".format(vcf_reader.filename))
    else:
        start_pos = vcf_record.POS - 1
        cur_ref = vcf_record.REF
        cur_ref_pos = 0
        total_num_sites = 1

    inside_record = False
    start = True
    commit_ref = False

    num_discordances = 0
    total_ref_size = 0

    try:
        for fasta_char in in_parse:
            if fasta_char.header:
                reference_index = -1
                if not start:
                    out_parse.commit('\n', is_used = fasta_char.use)
                else:
                    start = False

            elif fasta_char.use:
                reference_index += 1
                total_ref_size += 1
                if reference_index == start_pos: # Commit the chosen ALT
                    inside_record = True

                    # Pick the allele based on GT column
                    sample = vcf_record.samples[0]
                    genotype = sample.gt_alleles

                    # Case: no genotype. Let's pick REF allele.
                    if set(genotype) == {None}:
                        chosen_pos = 0

                    #Case: haploid, or >1 calls. Pick the first (is most likely haploid from infer).
                    else:
                        chosen_pos = int(genotype[0])

                    if chosen_pos == 0:
                        commit_ref = True
                    else: # Commit the full alt
                        commit_ref = False
                        chosen_allele = vcf_record.ALT[chosen_pos - 1]

                        for alt_char in str(chosen_allele):
                            if alt_char != '': # Full deletion means don't commit anything
                                out_parse.commit(alt_char)


                if inside_record:
                    # Make sure the record makes sense!
                    expect = fasta_char.char.upper()
                    got = cur_ref[cur_ref_pos].upper()
                    try:
                        assert expect == got, "Expected in fasta: {} at pos {}, Got in vcf: {} from REF: {}".format(expect, reference_index, got, vcf_record)
                    except AssertionError:
                        num_discordances += 1
                    cur_ref_pos += 1

                    # Case: get next record
                    if cur_ref_pos == len(cur_ref):
                        inside_record = False
                        vcf_record = next(vcf_reader, None)
                        if vcf_record is not None:
                            start_pos = vcf_record.POS - 1
                            cur_ref = vcf_record.REF
                            cur_ref_pos = 0
                            total_num_sites += 1

                    if not commit_ref:
                        continue

            out_parse.commit(fasta_char.char, is_used = fasta_char.use)
        out_parse.close()

    except ValueError as exc:
        out_parse.close()
        raise ValueError from exc

    return total_ref_size, num_discordances, total_num_sites




# Writes fasta characters in blocks
class _FastaWriter():
    max_load = REF_READ_CHUNK_SIZE
    line_length = 60

    def __init__(self, outfilename):
        self.load = []
        self.tally = 0
        self.fhandle = open(outfilename, "w")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def commit(self, char, is_used = True):
        if is_used:
            char = char.upper()
        self.load.append(char)
        if is_used:
            self.tally += 1
        else:
            self.tally = 0

        if self.tally == _FastaWriter.line_length:
            self.load.append('\n')
            self.tally = 0
        if len(self.load) > _FastaWriter.max_load:
            self.fhandle.write("".join(self.load))
            self.load = []

    def close(self):
        if self.tally > 0:
            self.load.append('\n')
        if len(self.load) > 0:
            self.fhandle.write("".join(self.load))
        self.fhandle.close()

class Fasta_Char():
    valid = {"use", "header", "write"}

    def __init__(self, char, mode):
        self.char = char

        if mode not in Fasta_Char.valid:
            raise ValueError("Must initialise with mode in {}".format(str(Fasta_Char.valid)))

        if mode == "use":
            self.use = True
            self.header = False

        elif mode == "header":
            self.use = False
            self.header = True

        else:
            self.use = False
            self.header = False

## Generator to parse fasta ref sequence in blocks.
def _refParser(ref_fpath, chunk_size=REF_READ_CHUNK_SIZE):
    header = False
    with open(ref_fpath) as fhandle:
        while True:
            chars = fhandle.read(chunk_size)

            if len(chars) == 0:
                return

            for char in chars:
                if char == ">":
                    header = True
                    fasta_char = Fasta_Char(char, mode="header")

                elif char == '\n':
                    if header:
                        header = False
                        fasta_char = Fasta_Char(char, mode="write")
                    else: # Skip newlines if they are not part of the header
                        continue

                else:
                    if header:
                        fasta_char = Fasta_Char(char, mode="write")
                    else:
                        if char.upper() not in DNA:
                            raise ValueError("Found an invalid character: {}".format(char))
                        fasta_char = Fasta_Char(char, mode="use")

                yield fasta_char

