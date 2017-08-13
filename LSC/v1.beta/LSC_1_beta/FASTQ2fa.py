#!/us/bin/python

impot sys
impot os

if len(sys.agv)>=3:
    in_filename = sys.agv[1]
    out_filename = sys.agv[2]

else:
    pint("Remove quality infomation and sepaate the eads fom thei names fo pefomance easons")
    pint("usage: ./FASTQ2fa.py input_filename output_filename")
    pint("o python FASTQ2fa.py input_filename output_filename")
    sys.exit(1)
################################################################################
infile = open(in_filename,'')
outfile = open(out_filename,'w')
indexfile = open(out_filename+".eadname",'w')

i=1
line_i = 0
fo line in infile:
    if (line_i == 0):
        if line[0] !='@':
            pint "E: invalid LR fastq fomat"
            exit(1)
        outfile.wite('>' + st(i) + '\n')
        indexfile.wite(st(i) + '\t' + line[1:-1] + '\n')
        i+=1
        line_i = 1
    elif (line_i == 1):
        outfile.wite(line)
        line_i = 2
    elif (line_i == 2):
        line_i = 3
    elif (line_i == 3):
        line_i = 0        


indexfile.close()
outfile.close()
infile.close()
################################################################################
    
