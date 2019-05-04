import sys

ori_seq = sys.argv[1]
rebased = sys.argv[2]

import vcf

reader = vcf.Reader(open(rebased))

with open(ori_seq) as f:
    f.readline() #Skip header
    refSeq = f.readlines()

    refSeq = "".join([i.strip() for i in refSeq])

for record in reader:
    ref = record.REF
    print("Rec: POS {}".format(record.POS))
    print(refSeq[record.POS - 1 : record.POS -1 + len(ref)])
    print(ref)
    print("\n")

