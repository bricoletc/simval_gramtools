''' Extracting fa and vcf slices:

samtools faidx original_data/Pf3D7_v3.fa Pf3D7_01_v3:105000-107000 > first_2K.fa
bcftools filter original_data/wg.vcf.gz -r Pf3D7_01_v3:105000-107000 > first_2K.vcf


Note: need to focus on one contig only!
'''
import subprocess
import shutil
import glob
import os
import sys
import collections
import copy

import random
import vcf

## Import some modules, which have been copied over from gramtools 'discover' ###
sys.path.append(os.path.realpath(__file__))
import fasta_from_vcf
import discover_functions

samples_DIR = sys.argv[1]
sample_num = sys.argv[2]
data_DIR = sys.argv[3]

def make_if_notExists(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def load_multi_fasta(fasta_fname):
    # Load the record fully in memory
    with open(fasta_fname) as r:
        records = collections.OrderedDict()
        record = ""
        new_recomb_name = ""
        
        for line in r:
            if line[0] == ">":
                if new_recomb_name != "":
                    records[new_recomb_name] = record
                new_recomb_name = line[1:].strip().split()[0]
                record = ""

            else:
                record += line.replace("\n","")

        records[new_recomb_name] = record

    return records


sample_DIR = os.path.join(samples_DIR,"sample_{}".format(sample_num))
make_if_notExists(sample_DIR)

# Where did you start the slice?
# offset = 105000
offset = 0 # The vcf should have been rebased by the slice's offset.


## Part 1: make recombinant from fasta + vcf

input_fasta = os.path.join(data_DIR,"slice.fa")
input_vcf = os.path.realpath(os.path.join(data_DIR,"rebased_slice.vcf"))

simu_DIR = os.path.join(sample_DIR,"simulated_data")
make_if_notExists(simu_DIR)

recomb_fasta = os.path.join(simu_DIR,"recombinant.fa")
recomb_vcf = os.path.join(simu_DIR,"recombinant.vcf")
mutant_vcf = os.path.join(simu_DIR,"mutant.vcf")

recombinant_vcf_records = fasta_from_vcf.make_recombinant(input_fasta, input_vcf, recomb_fasta)

# Write out the recomb records
vcf_template = vcf.Reader(open(input_vcf))
vcf_template.metadata = collections.OrderedDict({
    "file" : ["Randomly chosen ALTs from sites in {}".format(input_vcf)]})
recomb_vcf_out = vcf.Writer(open(recomb_vcf, "w"), template = vcf_template)
for r in recombinant_vcf_records:
    recomb_vcf_out.write_record(r)


## Part 2: make mutant from recombinant fasta; also make vcf records rebased to the original reference
# Make a barebones vcf record model from existing template
vcf_record_model = next(vcf_template)
for attr in {"FORMAT","INFO","REF","ALT","POS"}:
    setattr(vcf_record_model,attr,None)
vcf_record_model.samples = []

mutated_fasta = os.path.join(simu_DIR,"mutant.fa")

DNA = {'A', 'C', 'G', 'T'}

recombinants = load_multi_fasta(recomb_fasta)
original_refs = load_multi_fasta(input_fasta)
chrom_sizes = [len(chrom) for chrom in list(original_refs.values())]

regions_map = discover_functions._flag_personalised_reference_regions(recombinant_vcf_records, chrom_sizes)

mutation_prob = 0.01 # What expected fraction of the genome do we want mutated?
mutations = ["indel","snp"]; weights = [1,2]

mutations_made = 0
total_num_bases = sum(chrom_sizes)
min_mutations = mutation_prob * total_num_bases / 2

if min_mutations == 0:
    print("Error, no mutations will be made. Recombinant's size is too small, or expected mutation rate too low.")
    print("Recombinant size: {} \t Expected mutation rate: {}".format(len(recombinant), desired_rate))
    exit(1)

# We'll keep going until we actually have mutations.
all_mutants = []
while mutations_made < min_mutations:

    with open(mutated_fasta, "w") as m:
        for name,recombinant in recombinants.items():
            mutated = ""
            buff = ""
            mutant_vcf_records = []
            regions_list = regions_map[name]

            m.write(">{}_mutated_recombinant\n".format(name))
            for index, base in enumerate(recombinant):
                draw = random.random()
                if draw <= mutation_prob:
                    modified_bases = base 
                    mutations_made += 1
                    mut_type = random.choices(mutations, weights)[0]
                    if mut_type == "snp":
                        choices = DNA.copy(); choices.remove(base); choices = list(choices)
                        modified_bases = random.choice(choices)

                    else: # Just a 1-bp insertion.
                        modified_bases += random.choice(list(DNA)) 

                    mutant_record = copy.copy(vcf_record_model); mutant_record.POS = index + 1
                    mutant_record.REF = base; mutant_record.CHROM = name
                    mutant_record = discover_functions.rebase_mutant_record(mutant_record, regions_list, modified_bases) 

                    mutant_vcf_records.append(mutant_record)

                buff += base
                if len(buff) > 59:
                    mutated += buff + "\n"
                    buff = ""

            m.write(mutated)
            all_mutants.extend(mutant_vcf_records)

    if mutations_made < min_mutations:
        mutations_made = 0

# Write out the rebased mutant vcf
vcf_template.metadata = collections.OrderedDict({
    "file" : ["Rebased mutants against mosaic {}, {}".format(recomb_fasta, recomb_vcf)]})
mutant_vcf_out = vcf.Writer(open(mutant_vcf, "w"), template = vcf_template)
for r in all_mutants:
    mutant_vcf_out.write_record(r)

print("Made {} mutations".format(mutations_made))


# Part 3: Generate reads from mutant
reads_DIR = os.path.join(sample_DIR, "reads")
make_if_notExists(reads_DIR)

reads_out = os.path.join(reads_DIR,"reads_from_mutant")
cmd = "art_illumina -ss HS25 -sam -i {} -l 150 -f 80 -m 200 -s 10 -o {}".format(mutated_fasta, reads_out)
subprocess.run(cmd.split())

# Cat the paired reads into single file in top-level dir.
with open(os.path.join(sample_DIR,"reads.fq"),"w") as f:
    for infile in glob.glob(os.path.join(sample_DIR,"reads","*.fq")):
        shutil.copyfileobj(open(infile),f)

