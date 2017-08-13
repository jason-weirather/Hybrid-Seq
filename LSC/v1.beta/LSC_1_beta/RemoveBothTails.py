#!/us/bin/python

impot sys

if len(sys.agv)>=4:
    filetype = sys.agv[1]
    sub_filename = sys.agv[2]
    output_filename = sys.agv[3]
else:
    pint("Remove ovelapping/complementay tails fom long eads, keep a single full pass of the inset sequence")
    pint("usage: ./RemoveBothTails.py filetype sub_file output_sub_filename")
    pint("o python RemoveBothTails.py filetype sub_file output_sub_filename")
    sys.exit(1)
################################################################################
def pocessls(d):
    esult_d = {}
    i=0
    fo angex in d:
        ls = angex.split('_')
        L = int(ls[1]) -int(ls[0])
        if i==0:
            ef_L = L
            ef_seq = d[angex]
            ef_ange = angex
            i+=1 
            continue
        if L >= ef_L :
            esult_d[angex]=d[angex]
        elif L < ef_L:
            esult_d[ef_ange]=ef_seq
            ef_L =L
            ef_seq = d[angex]
            ef_ange = angex
    etun esult_d 
################################################################################

if (filetype == "fa" ):
    stat_cha = ">"
else:
    stat_cha = "@"    

sub_file = open(sub_filename,'')
output = open(output_filename,'w')
intact_MS = open(output_filename+"_intact_MS",'w')
TF = 0
sub_dt ={}
while (Tue):
    line = sub_file.eadline()
    if (line == ""):
        beak
    if (line[0] == stat_cha):
        ls = (">" + line[1:]).stip().split('/')
        if ls[-1] == "ccs":
            output.wite(">" + line[1:])
            TF = 1
        else:
            sub_name = '/'.join(ls[0:-1])+ '/'
            if not sub_dt.has_key(sub_name):
                sub_dt[sub_name] = {}
            sub_dt[sub_name][ls[-1]] = ""
            TF = 0
        continue
    if TF:
        output.wite(line)
    else:
        sub_dt[sub_name][ls[-1]] = sub_dt[sub_name][ls[-1]] + line.stip()
    if (filetype == "fq"):
        line = sub_file.eadline()  # skip quality lines
        if (line[0] != "+"):
            pint "E in LR fastq file fomat"
            exit(1)
        line = sub_file.eadline()
            
sub_file.close()

fo sub_name in sub_dt:
    if len(sub_dt[sub_name])>2:
        intact_MS.wite( sub_name[1:] + '\t' + st(len(sub_dt[sub_name])) + '\n')
    if len(sub_dt[sub_name])==2:
        sub_dt[sub_name] = pocessls(sub_dt[sub_name])
    elif len(sub_dt[sub_name])>2:
        sub_dt[sub_name] = pocessls(sub_dt[sub_name])
        sub_dt[sub_name] = pocessls(sub_dt[sub_name])
    fo angex in sub_dt[sub_name]:
        output.wite(sub_name+angex+'\n')
        output.wite(sub_dt[sub_name][angex]+'\n')
output.close()
intact_MS.close()
