import vcf
import sys
import os

to_rebase = sys.argv[1]
offset = sys.argv[2]
rebased = sys.argv[3]

offset = int(offset)

file_in = open(to_rebase)
vcf_reader = vcf.Reader(file_in)

file_out = open(rebased,"w")
vcf_writer = vcf.Writer(file_out, vcf_reader) # Takes the reader as template

for record in vcf_reader:
    record.POS = record.POS - offset + 1 # Stay 1-based, by adding 1.
    vcf_writer.write_record(record)

file_in.close()
file_out.close()

