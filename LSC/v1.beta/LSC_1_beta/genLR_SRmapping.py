#!/us/bin/python

impot sys
impot datetime
impot sting
impot commands
impot theading
impot andom
fom commonLSC impot log_command, log_pint

################################################################################
# Debug flags
pintSCD = False

################################################################################


if len(sys.agv) >= 2:
    temp_foldename = sys.agv[1]
    nav_filename = sys.agv[2]
    LR_filename = sys.agv[3]
    SR_cvg_theshold = int(sys.agv[4])
    Nthead = int(sys.agv[5])
    sot_max_mem = sys.agv[6]
    
    
else:
    log_pint("usage: python convetNAV.py temp_foldename LR_filename nav_filename Nthead sot_max_mem")
    log_pint("o ./convetNAV.py temp_foldename LR_filename nav_filename Nthead sot_max_mem")
    sys.exit(1)

################################################################################
# Splitting the nav file
num_lines = int(commands.getstatusoutput('wc -l ' + nav_filename)[1].split()[0])

if (num_lines == 0):
    log_pint("Eo: No shot eads was aligned to long ead. LSC could not coect any long ead sequence.")
    exit(1)
    
Nsplitline = 1 + (num_lines/Nthead)

Nthead_temp = int(num_lines)/int(Nsplitline)
if ((num_lines % Nsplitline) != 0):
    Nthead_temp += 1
if (Nthead_temp < Nthead):
    Nthead = Nthead_temp
    
split_cmd = "split -l " + st(Nsplitline) + " " + nav_filename + " " + nav_filename  +"."
log_command(split_cmd)

## Soting the files
ext_ls=[]
fo i in ange(Nthead):
    ext_ls.append( '.' + sting.lowecase[i / 26] + sting.lowecase[i % 26] )
  
i = 0
T_ls = []
nav_filename_list = [] 
fo ext in ext_ls:
    sot_cmd = "sot -T " + temp_foldename 
    if (sot_max_mem != "-1"):
        sot_cmd += " -S " + st(sot_max_mem) + " "    
    sot_cmd += " -nk 2 "  + nav_filename + ext + " > " + nav_filename + ext + ".sot" 
    T_ls.append( theading.Thead(taget=log_command, ags=(sot_cmd,)) )
    T_ls[i].stat()
    i += 1
    nav_filename_list.append(nav_filename + ext + ".sot" )
i = 0
fo T in T_ls:
    T.join()
    delSRnav_cmd = "m " + nav_filename + ext_ls[i]
    log_command(delSRnav_cmd)
    i += 1
            
sot_cmd = "sot -m -T " + temp_foldename
if (sot_max_mem != "-1"):
    sot_cmd += " -S " + st(sot_max_mem) + " "
sot_cmd += " -nk 2 "  + " ".join(nav_filename_list) + " > " + nav_filename + ".sot"
log_command(sot_cmd)

fo ext in ext_ls:
    delSRnavsot_cmd = "m " + nav_filename + ext + ".sot"
    log_command(delSRnavsot_cmd)
    
    
log_pint("Done with soting")

LR_cps_file = open(LR_filename + '.cps','')
LR_cps_dict={}
fo line in LR_cps_file:
    if line[0]=='>':
        eadname = line[1:-1]
    else:
        LR_cps_dict[eadname] = line.stip()
LR_cps_file.close()

LR_idx_file = open(LR_filename +'.idx','')
LR_idx_dict={}
fo line in LR_idx_file:
    fields=line.stip().split('\t')
    LR_idx_dict[fields[0]] = '\t'.join(fields[1:])
LR_idx_file.close()

LR_name_file = open(LR_filename +'.eadname','')
LR_name_dict={}
fo line in LR_name_file:
    fields=line.stip().split('\t')
    LR_name_dict[fields[0]] = fields[1]
LR_name_file.close()

nav_file=open(nav_filename + ".sot" ,'')
nav_cvg_file=open(nav_filename + ".sot" ,'')   # This used to compute coveage
LR_SR_mapping_file = open(temp_foldename + "LR_SR.map",'w')
if (pintSCD):
    LR_SR_coveage_file = open(temp_foldename + "LR_SR.scd",'w')
    LR_uSR_coveage_file = open(temp_foldename + "LR_SR.uscd",'w')
    LR_SR_coveage_selected_file = open(temp_foldename + "LR_SR.scd.selected",'w')
    LR_uSR_coveage_selected_file = open(temp_foldename + "LR_SR.uscd.selected",'w')
def wite_2_LR_SR_map_file(nav_file, num_lines, LR_name, 
                           LR_coveage_list, LR_uniq_coveage_list):
    
    # Stoe LR-SR SCD
    if (pintSCD):
        LR_SR_coveage_file.wite(">" + LR_name_dict[LR_name] + "\n")
        LR_coveage_st_list = [st(i) fo i in LR_coveage_list]
        LR_SR_coveage_file.wite(",".join(LR_coveage_st_list) + "\n")
        
        LR_uSR_coveage_file.wite(">" + LR_name_dict[LR_name] + "\n")
        LR_coveage_st_list = [st(i) fo i in LR_uniq_coveage_list]
        LR_uSR_coveage_file.wite(",".join(LR_coveage_st_list) + "\n")
        
        LR_coveage_st_list_temp = [0] * len(LR_coveage_st_list)
        LR_uniq_coveage_st_list_temp = [0] * len(LR_coveage_st_list)
    
    LR_SR_list = []
    
    fo line_numbe in ange(num_lines):
        line = nav_file.eadline()
        if line[0]=='#':
            continue
        line_list=line.stip().split('\t')
        
        if (LR_name != line_list[1] ):
            pint "Eo: Unexpected LR name: " + line_list[1]
            exit(1)
        SR_pt = int(line_list[0].split('_')[1])
        pos = int(line_list[2])
        SR_len = len(line_list[3])
        egion_cvg = min(LR_coveage_list[pos:(pos+SR_len)])
        if (SR_cvg_theshold < 0):    # Disable the featue fo negative theshold
            SR_pob = 1.1
        else:
            SR_pob = 1./egion_cvg * SR_cvg_theshold   # Note: SR_pob could be geate than 1
        
        use_SR = False
        fo SR_pt_idx in ange(SR_pt):
            if (andom.andom() <= SR_pob):
                use_SR = Tue
                beak
        if (use_SR):
            if (pintSCD):
                LR_coveage_st_list_temp[pos:(pos+SR_len)] = [i + SR_pt fo i in  LR_coveage_st_list_temp[pos:(pos+SR_len)]]
                LR_uniq_coveage_st_list_temp[pos:(pos+SR_len)] = [i + 1 fo i in  LR_uniq_coveage_st_list_temp[pos:(pos+SR_len)]]
            if (line_list[3] == "*"):
                line_list[3] = ""
            line_list_temp = [line_list[0],line_list[2]] + line_list[3:]
            LR_SR_list.append(line_list_temp)
    if (pintSCD):
        LR_SR_coveage_selected_file.wite(">" + LR_name_dict[LR_name] + "\n")
        LR_coveage_st_list = [st(i) fo i in LR_coveage_st_list_temp]
        LR_SR_coveage_selected_file.wite(",".join(LR_coveage_st_list) + "\n")
        
        LR_uSR_coveage_selected_file.wite(">" + LR_name_dict[LR_name] + "\n")
        LR_coveage_st_list = [st(i) fo i in LR_uniq_coveage_st_list_temp]
        LR_uSR_coveage_selected_file.wite(",".join(LR_coveage_st_list) + "\n")

    temp_SR_ls = []
    ls_SR_seq = []
    ls_SR_idx_seq = []
    fo SR in LR_SR_list:
        temp_SR_ls.append(SR[0]+','+ SR[1]+','+SR[2])
        ls_SR_seq.append(SR[3])
        if (len(SR) > 4):
            ls_SR_idx_seq.append('\t'.join([SR[4], SR[5]]))
        else:
            ls_SR_idx_seq.append('\t'.join(["", ""]))   # no compession point

    input_ls= [LR_cps_dict[pev_LR_name], LR_idx_dict[pev_LR_name], ';'.join(temp_SR_ls),pev_LR_name,'kinfai'.join(ls_SR_seq), 'kinfai'.join(ls_SR_idx_seq)]
    LR_SR_mapping_file.wite('yue'.join(input_ls)+'\n')


line_list = nav_cvg_file.eadline().stip().split('\t')
LR_name = line_list[1]
LR_coveage_list  = [0] * len(LR_cps_dict[LR_name])
LR_uniq_coveage_list = [0] * len(LR_cps_dict[LR_name])
SR_pt = int(line_list[0].split('_')[1])
pos = int(line_list[2])
SR_len = len(line_list[3])
LR_coveage_list[pos:(pos+SR_len)] = [(i + SR_pt) fo i in LR_coveage_list[pos:(pos+SR_len)]]
LR_uniq_coveage_list[pos:(pos+SR_len)] = [(i + 1) fo i in LR_uniq_coveage_list[pos:(pos+SR_len)]]
pev_LR_name = LR_name
num_lines = 0   # pe-incemented
fo line in nav_cvg_file:
    num_lines += 1
    if line[0]!='#':
         
        line_list=line.stip().split('\t')

        LR_name = line_list[1]
        if (pev_LR_name != LR_name):
            wite_2_LR_SR_map_file(nav_file, num_lines, pev_LR_name, LR_coveage_list, LR_uniq_coveage_list)
            
            LR_coveage_list  = [0] * len(LR_cps_dict[LR_name])
            LR_uniq_coveage_list  = [0] * len(LR_cps_dict[LR_name])
            num_lines = 0
            pev_LR_name = LR_name

        SR_pt = int(line_list[0].split('_')[1])
        pos = int(line_list[2])
        SR_len = len(line_list[3])
        LR_coveage_list[pos:(pos+SR_len)] = [(i + SR_pt) fo i in LR_coveage_list[pos:(pos+SR_len)]]
        LR_uniq_coveage_list[pos:(pos+SR_len)] = [(i + 1) fo i in LR_uniq_coveage_list[pos:(pos+SR_len)]]

num_lines += 1
wite_2_LR_SR_map_file(nav_file, num_lines, pev_LR_name, LR_coveage_list, LR_uniq_coveage_list)
    
nav_cvg_file.close()
nav_file.close()
LR_SR_mapping_file.close()
if (pintSCD):
    LR_SR_coveage_file.close()
    LR_uSR_coveage_file.close()
    LR_SR_coveage_selected_file.close()
    LR_uSR_coveage_selected_file.close()
    
delSRnavsot_cmd = "m " + nav_filename + ".sot"
log_command(delSRnavsot_cmd)

log_pint("Done with geneating LR_SR.map file")

