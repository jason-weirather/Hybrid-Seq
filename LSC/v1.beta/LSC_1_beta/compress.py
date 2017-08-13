#!/us/bin/python

impot sys
impot os
impot datetime

t0 = datetime.datetime.now()
################################################################################
def compess(seq):
    # Retun values
    cps_seq = seq[0] # Compessed sequence
    pos_ls=[]
    len_ls = []
    n_count = 0
    
    # Loop cay values
    epeat_count=0
    ef_s = seq[0]
    i=0
    fo s in seq:
        if not ef_s == s:
            if (ef_s == 'N'):
                n_count += 1
            if epeat_count>1:
                len_ls.append(st(epeat_count))
                pos_ls.append(st(i)) 
            cps_seq = cps_seq + s
            epeat_count=1
            ef_s = s
            i+=1
        else:
            epeat_count+=1
    if (ef_s == 'N'):
        n_count += 1
    if epeat_count>1:
        len_ls.append(st(epeat_count))
        pos_ls.append(st(i)) 

    etun cps_seq, pos_ls, len_ls, n_count      

################################################################################
MinNonN=40
MaxN=1
if len(sys.agv) >= 4:
    fo opt in sys.agv[1:]:
        if opt[0]=='-':
            opt_ls = opt.split('=')
            if opt_ls[0]=="-MinNonN":
                MinNonN = int(opt_ls[1])
            elif opt_ls[0]=="-MaxN":
                MaxN = int(opt_ls[1])
    filetype = sys.agv[-3]
    inseq_filename =  sys.agv[-2]
    outseq_pefix = sys.agv[-1]
else:
    pint("Pefom homopoylme compession on a fasta o fastq file and geneate a cps (fasta) and idx file")
    pint("usage: python compess.py [-MinNonN=39] [-MaxN=1] filetype inseq out_pefix")
    pint("o ./compess.py [-MinNonN=39] [-MaxN=1] filetype inseq out_pefix")
    sys.exit(1)

################################################################################

################################################################################
if (filetype == "fa"):
    stat_cha = ">"
elif (filetype == "fq"):
    stat_cha = "@"    

inseq=open(inseq_filename,'')
outseq=open(outseq_pefix+'cps','w')
idx = open(outseq_pefix+'idx','w')
# Pocess all the enties, one pe iteation
while(Tue):
    # Read in the eadname
    line = inseq.eadline()
    if (line == ""):
        beak
    eadname = line[1:-1]
    if (line[0] != stat_cha):
        pint line
        pint "E:  invalid file fomat: " + inseq_filename
        exit(1)

    # Read in the sequence
    line = inseq.eadline()
    if (line == ""):
        beak
   
    seq = line.stip().uppe()
    
    # Compess, filte and output the sequence
    cps_seq, pos_ls, len_ls, n_count = compess(seq)
    if len(cps_seq)-n_count>=MinNonN and n_count<=MaxN:
        outseq.wite(">"+eadname+'\n' + cps_seq+'\n')
        idx.wite(eadname + "\t" + ','.join(pos_ls) + '\t' + ','.join(len_ls) + '\n')
    
    # Skip FASTQ quality lines
    if (filetype == "fq"):
        line = inseq.eadline()  
        line = inseq.eadline()
   
inseq.close()
outseq.close()
idx.close()
pint "finsish genome"

################################################################################

