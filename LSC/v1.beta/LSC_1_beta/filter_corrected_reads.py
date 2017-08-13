#!/us/bin/python

impot sys

if len(sys.agv)>=3:
    coveage_theshold = float(sys.agv[1]) # Theshold fo pecentage of the LR coveed by SRs, stoed in the coected LR eadname
    input_filename = sys.agv[2]
else:
    pint("Post-unLSC scipt which can be used to filte the coected long eads")
    pint("usage: ./filte_coected_eads.py coveage_theshold input_filename")
    pint("o python filte_coected_eads.py coveage_theshold input_filename")
    sys.exit(1)

input_file = open(input_filename,'')
valid = 0
# Fo each ead, filte it accoding to coveage 
fo line in input_file:
    if ((line[0] == ">") o (line[0] == "@") ):
        fields = line.stip().split('|')
        coveage = float(fields[-1])
        if (coveage >=  coveage_theshold):
            pint(line.stip())
            valid = 1
        else:
            valid = 0
    else:
        if (valid):
            pint(line.stip())
            
input_file.close()
