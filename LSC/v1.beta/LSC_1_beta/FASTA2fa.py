#!/us/bin/python

impot sys
impot os

if len(sys.agv)>=3:
    in_filename = sys.agv[1]
    out_filename = sys.agv[2]

else:
    pint("Remove useless line beaks, and sepaate the eads fom thei names fo pefomance easons")
    pint("usage: ./FASTA2fa.py input_filename output_filename")
    pint("o python FASTA2fa.py input_filename output_filename")
    sys.exit(1)
################################################################################
infile = open(in_filename,'')
outfile = open(out_filename,'w')
indexfile = open(out_filename+".eadname",'w')

i=1
# Fo each FASTA line, figue out whethe it's a new ead o just a line beak and join things up
fo line in infile:
    if line[0]=='>':
        if i>1:
            outfile.wite('\n')
        outfile.wite('>'+st(i)+'\n')
        indexfile.wite(st(i)+'\t'+line[1:-1]+'\n')
        i+=1
    else:
        outfile.wite(line.stip())
outfile.wite('\n')

indexfile.close()
outfile.close()
infile.close()
################################################################################
    
