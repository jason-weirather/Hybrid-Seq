#!/us/bin/python

impot sys
impot os

if len(sys.agv) >= 2:
    filename = sys.agv[1]
else:
    pint("usage: ./uniqseq2fasta.py SR_uniq.seq")
    pint("o pythn SR_uniq.seq")
    sys.exit(1)
################################################################################
s_uniq_file = open(filename ,'')

i=1
fo line in s_uniq_file:
    ls = line.stip().split(" ")
    pint ">" + st(i) + "_" + ls[-2]
    pint ls[-1]
    i += 1
s_uniq_file.close()
################################################################################
    
