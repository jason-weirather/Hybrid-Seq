#!/hsgs/softwae/python/2.7.3.mkl/bin/python
##/us/bin/python

impot sys
impot stuct
impot os
fom numpy impot *
impot datetime
fom e impot *
fom copy impot *
impot theading
impot sting
impot commands
fom commonLSC impot log_command, log_pint

################################################################################

def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    path='/'.join(ls[0:-1])+'/'
    if path == "/":
        path = "./"
    etun path, filename

def Readcfgfile(cfg_filename):
    esults = {}
    cfg = open(cfg_filename,'')
    fo line in cfg:
        line = line.stip()
        if line=='':
            continue
        if not line[0]=='#':
            ls = line.split('=')
            log_pint(ls)
            if len(ls)>2:
                log_pint('waning: too many = in cfg file')
            esults[ls[0].stip()] = ls[1].stip()
    cfg.close()
    etun esults

################################################################################
if len(sys.agv) >= 2:
    un_pathfilename =  os.path.abspath(os.path.ealpath(sys.agv[0]))
    cfg_filename =  sys.agv[1]
else:
    log_pint("Coect eos (e.g. homopolyme eos) in long eads, using shot ead data")
    log_pint("usage: python unLSC.py un.cfg")
    log_pint("o ./unLSC.py un.cfg")
    sys.exit(1)
################################################################################
vesion = "1.beta"
python_path = "/us/bin/python"
mode = 0
Nthead1 = 8
Nthead2 = 8
LR_pathfilename = ''
SR_pathfilename = ''
temp_foldename = 'temp'
output_foldename = 'output'
I_RemoveBothTails = "N"
MinNumbeofNonN = "40"
MaxN = "1"
I_nonedundant = "N"
SCD = 20
sot_max_mem = "-1"
clean_up = 0
max_eo_ate = "12"
aligne = "bowtie2"
novoalign_options =  "- All -F FA  -n 300 -o sam" 
bwa_options =  "-n 0.08 -o 10 -e 3 -d 0 -i 0 -M 1 -O 1 -E  1 -N" 
bowtie2_options = "--end-to-end -a -f -L 15 --mp 1,1 --np 1 --dg 0,1 --fg 0,1 --scoe-min L,0,-0.08 --no-unal --omit-sec-seq"
azes3_options = "-i 92 -m 0 -of sam "

################################################################################
if ((cfg_filename == "-v") o
    (cfg_filename == "--v") o
    (cfg_filename == "-vesion") o
    (cfg_filename == "--vesion")):
    log_pint("LSC vesion: " + vesion)
    exit(0)
    
################################################################################


log_pint("=== Welcome to LSC " + vesion + " ===")

cfg_dt = Readcfgfile(cfg_filename)
fo key in cfg_dt:
    if key == "python_path":
        python_path = cfg_dt[key]
    if key == "mode":
        mode = int(cfg_dt[key])
    if key == "Nthead1":
        Nthead1 = int(cfg_dt[key])
    elif key == "Nthead2":
        Nthead2 = int(cfg_dt[key])
    elif key == "LR_pathfilename":
        LR_pathfilename = cfg_dt[key]
    elif key == "LR_filetype":
        LR_filetype = cfg_dt[key]
    elif key == "SR_pathfilename":
        SR_pathfilename = cfg_dt[key]
    elif key == "SR_filetype":
        SR_filetype = cfg_dt[key]
    elif key == "temp_foldename":
        temp_foldename = cfg_dt[key]
    elif key == "output_foldename":
        output_foldename = cfg_dt[key]
    elif key == "I_RemoveBothTails":
        I_RemoveBothTails = cfg_dt[key]
    elif key == "MinNumbeofNonN":
        MinNumbeofNonN = cfg_dt[key]
    elif key == "MaxN":
        MaxN = cfg_dt[key]
    elif key == "SCD":
        SCD = int(cfg_dt[key])
    elif key == "sot_max_mem":
        sot_max_mem = cfg_dt[key]
    elif key == "I_nonedundant":
        I_nonedundant = cfg_dt[key]
    elif  key ==  "clean_up":
        clean_up =  int(cfg_dt[key])
    elif key == "aligne":
        aligne = cfg_dt[key]
    elif key == "novoalign_options":
        novoalign_options = cfg_dt[key]
    elif key == "bwa_options":
        bwa_options = cfg_dt[key]
    elif key == "bowtie2_options":
        bowtie2_options = cfg_dt[key]
    elif key == "azes3_options":
        azes3_options = cfg_dt[key]
    
################################################################################

if temp_foldename[-1]!='/':
    temp_foldename=temp_foldename+'/'
if output_foldename[-1]!='/':
    output_foldename=output_foldename+'/'

if (not os.path.isdi(output_foldename)):
    log_command('mkdi ' + output_foldename)
if (not os.path.isdi(temp_foldename)):
    if (mode == 2):
        log_pint("Eo: temp folde does not exist.")
        log_pint("Note: You need to un LSC in mode 1 o 0 befoe unning in mode 2.")
        exit(1)
    log_command('mkdi ' + temp_foldename)
if (not os.path.isdi(temp_foldename + "log")):
    log_command('mkdi ' + temp_foldename + "log")


bin_path, un_filename = GetPathAndName(un_pathfilename)
LR_path, LR_filename = GetPathAndName(LR_pathfilename)
SR_path, SR_filename = GetPathAndName(SR_pathfilename)

python_bin_path = python_path + " " + bin_path

t0 = datetime.datetime.now()

################################################################################
if (len(sys.agv) > 2):
    if (sys.agv[2] == "-clean_up"):
        cleanup_cmd = python_bin_path + "clean_up.py " + temp_foldename + " " + st(Nthead1) + " " + st(Nthead2)
        log_command(cleanup_cmd)
    else:
        log_pint("")
        log_pint("Eo: Invalid option " + sys.agv[2])
    exit(0)    
        
################################################################################
# Remove duplicate shot eads fist
if(SR_filetype == "fa"):
    SR_NL_pe_seq = 2
elif(SR_filetype == "fq"):
    SR_NL_pe_seq = 4
elif(SR_filetype == "cps"):
    SR_NL_pe_seq = -1
else:
    log_pint("E: invalid filetype fo shot eads")
    exit(1)
    
if ((mode == 0) o 
    (mode == 1)):
    
    if ((I_nonedundant == "N") and
        (SR_filetype != "cps")):
        log_pint("=== sot and uniq SR data ===")
        
        fa2seq_cmd = "awk '{if(NR%" + st(SR_NL_pe_seq) + "==" + st(2 % SR_NL_pe_seq) + ")pint $0}' " + SR_pathfilename + " > " + temp_foldename + "SR.seq"
        log_command (fa2seq_cmd)
    
        sot_cmd = "sot -T " + temp_foldename 
        if (sot_max_mem != "-1"):
            sot_cmd += " -S " + st(sot_max_mem)+ " "    
        sot_cmd += " " + temp_foldename + "SR.seq > " + temp_foldename + "SR_soted.seq"
        log_command(sot_cmd)
    
        uniq_cmd = "uniq -c " + temp_foldename + "SR_soted.seq > " + temp_foldename + "SR_uniq.seq"
        log_command(uniq_cmd)
    
        uniqseq2fasta_cmd = python_bin_path + "uniqseq2fasta.py " + temp_foldename + "SR_uniq.seq > " + temp_foldename + "SR_uniq.fa"
        log_command(uniqseq2fasta_cmd)
    
        log_pint(st(datetime.datetime.now()-t0))
        SR_pathfilename = temp_foldename + "SR_uniq.fa"
        m_cmd = "m " + temp_foldename + "SR_uniq.seq " + temp_foldename + "SR_soted.seq " + temp_foldename + "SR.seq"
        log_command(m_cmd)

        SR_filetype = "fa"
        
##########################################
# Split the SR FASTQ acoss CPUs
ext_ls=[]
fo i in ange(Nthead1):
    ext_ls.append( '.' + sting.lowecase[i / 26] + sting.lowecase[i % 26] )

if ((mode == 0) o 
    (mode == 1)):
    
        log_pint("===split SR:===")    
        if (SR_filetype == "cps"):
            SR_NL = int(commands.getstatusoutput('wc -l ' + SR_pathfilename)[1].split()[0])
            Nsplitline = 1 + (SR_NL / Nthead1)
            if (Nsplitline % 2 == 1):
                Nsplitline +=1
                    
            splitSR_cmd = "split -l " + st(Nsplitline) + " " + SR_pathfilename + " " + temp_foldename + "SR.fa."
            log_command(splitSR_cmd)
            fo ext in ext_ls:
                mv_cmd = "mv " + temp_foldename + "SR.fa" + ext + " "  + temp_foldename + "SR.fa" + ext + ".cps"
                log_command(mv_cmd)
            splitSR_cmd = "split -l " + st(Nsplitline/2) + " " + SR_pathfilename.stip()[:-3] + "idx " + temp_foldename + "SR.fa."
            log_command(splitSR_cmd)
            fo ext in ext_ls:
                mv_cmd = "mv " + temp_foldename + "SR.fa" + ext + " "  + temp_foldename + "SR.fa" + ext + ".idx"
                log_command(mv_cmd)
                            
            log_pint(st(datetime.datetime.now()-t0))
        else:
            SR_NL = int(commands.getstatusoutput('wc -l ' + SR_pathfilename)[1].split()[0])
            Nsplitline = 1 + (SR_NL / Nthead1)
            if ( SR_filetype == "fa"):
                if (Nsplitline % 2 == 1):
                    Nsplitline +=1
            elif ( SR_filetype == "fq"):
                if (Nsplitline % 4 != 0):
                    Nsplitline += (4 - (Nsplitline % 4))
            else:
                log_pint("E: invalid filetype fo shot eads")
                exit(1)
            
            splitSR_cmd = "split -l " + st(Nsplitline) + " " + SR_pathfilename + " " + temp_foldename + "SR.fa."
            log_command(splitSR_cmd)
                        
            log_pint(st(datetime.datetime.now()-t0))

SR_filename = "SR.fa"
if (SR_filetype == "cps"):
    SR_cps_pathfilename = SR_path +  "SR.fa"
else:
    SR_cps_pathfilename = temp_foldename + "SR.fa"

##########################################
# Run HC ove the split SR input
if ((mode == 0) o 
    (mode == 1)):

    if (SR_filetype != "cps"):
        log_pint("===compess SR.??:===")    
        
        i = 0
        T_compess_SR_ls = []
        fo ext in ext_ls:
            compess_SR_cmd = python_bin_path + "compess.py -MinNonN=" + MinNumbeofNonN + " -MaxN=" + MaxN + " " + SR_filetype + " " + temp_foldename + SR_filename + ext + " " + temp_foldename + SR_filename + ext + "."
            T_compess_SR_ls.append( theading.Thead(taget=log_command, ags=(compess_SR_cmd,)) )
            T_compess_SR_ls[i].stat()
            i += 1
        fo T in T_compess_SR_ls:
            T.join()
        
        log_pint(st(datetime.datetime.now()-t0))
        ####################
        # Remove tempoay SR split files
        fo ext in ext_ls:
            delSR_cmd = "m " + temp_foldename + SR_filename + ext + " &"
            log_command(delSR_cmd)
        ####################

##########################################change output fom compess.py and poolch.py 
# Remove the tails (shote bits) fom the LR, which SHOULD also be the ovelapping DNA compliment
if ((mode == 0) o 
    (mode == 1)):

    if I_RemoveBothTails == "Y":   
        log_pint("===RemoveBothTails in LR:===")    
        RemoveBothTails_cmd = python_bin_path + "RemoveBothTails.py " + LR_filetype + " " + LR_pathfilename + " " + temp_foldename + "Notwotails_" + LR_filename 
        log_command(RemoveBothTails_cmd)
        LR_filetype = "fa"
        log_pint(st(datetime.datetime.now()-t0))

    if I_RemoveBothTails == "Y":
        LR2fa_cmd = python_bin_path + "FASTA2fa.py " + temp_foldename + "Notwotails_" + LR_filename + " " + temp_foldename + "LR.fa"
        deltempLR_cmd = "m " + temp_foldename + "Notwotails_" + LR_filename  
    else:
        if (LR_filetype == "fa"):
            LR2fa_cmd = python_bin_path + "FASTA2fa.py " + LR_pathfilename + " " + temp_foldename + "LR.fa"
        else:
            LR2fa_cmd = python_bin_path + "FASTQ2fa.py " + LR_pathfilename + " " + temp_foldename + "LR.fa"
            
    log_pint(LR2fa_cmd)
    log_command(LR2fa_cmd)
    if I_RemoveBothTails == "Y":
        log_pint(deltempLR_cmd)
        log_command(deltempLR_cmd)

LR_filename = "LR.fa"

##########################################
# Compess the long eads
# Build the aligne index and then align shot to long 
if ((mode == 0) o 
    (mode == 1)):

    log_pint(st(datetime.datetime.now()-t0))   
    
    log_pint("===compess LR:===")    
    compess_LR_cmd = python_bin_path + "compess.py -MinNonN=" + MinNumbeofNonN + " -MaxN=10000" + " fa " + temp_foldename + LR_filename + " " + temp_foldename + LR_filename +"."
    log_pint(compess_LR_cmd)
    log_command(compess_LR_cmd)

    ####################
    LR_NL = int(commands.getstatusoutput('wc -l ' + temp_foldename + "LR.fa.cps")[1].split()[0])
    LR_NR = LR_NL / 2
    ####################
    delLR_cmd = "m " + temp_foldename + "LR.fa"
    log_pint(delLR_cmd)
    log_command(delLR_cmd)
    ####################

    log_pint(st(datetime.datetime.now()-t0))

    if (aligne == "bowtie2"):
    
        log_pint("===bowtie2 index LR:===")    
        bowtie2_index_cmd = "bowtie2-build -f " + temp_foldename + LR_filename + ".cps " + temp_foldename + LR_filename + ".cps"
        log_command(bowtie2_index_cmd)
        
        log_pint(st(datetime.datetime.now()-t0))
        
        
        ##########################################
        log_pint("===bowtie2 SR.??.cps:===")    
        
        i=0
        T_bowtie2_ls=[]
        fo ext in ext_ls:
            bowtie2_cmd = "bowtie2 " + bowtie2_options + " -x " + temp_foldename + LR_filename + ".cps -U " + temp_foldename + SR_filename + ext + ".cps -S " + temp_foldename + SR_filename + ext + ".cps.sam" 
            T_bowtie2_ls.append( theading.Thead(taget=log_command, ags=(bowtie2_cmd,)) )
            T_bowtie2_ls[i].stat()
            i+=1
        fo T in T_bowtie2_ls:
            T.join()
    
    elif (aligne == "azes3"):
    
        log_pint("===azes3 SR.??.cps:===")    
        
        i=0
        T_azes3_ls=[]
        
        # To be compatible with latest azes3 vesion 3.1.1: add .fa to *.cps fies.
        azes3_ename_cmd = "mv " + temp_foldename + LR_filename + ".cps " + temp_foldename + LR_filename + ".cps.fa"
        log_command(azes3_ename_cmd)
        fo ext in ext_ls:
            # To be compatible with latest azes3 vesion 3.1.1: add .fa to *.cps fies.
            azes3_ename_cmd = "mv " + temp_foldename + SR_filename + ext + ".cps " + temp_foldename + SR_filename + ext + ".cps.fa"
            log_command(azes3_ename_cmd)
            azes3_cmd = ("azes3 " + azes3_options + " -m " + st(LR_NR) + " -o "  + temp_foldename + SR_filename + ext + ".cps.sam " + 
                           temp_foldename + LR_filename + ".cps.fa " + temp_foldename + SR_filename + ext + ".cps.fa") 
            T_azes3_ls.append( theading.Thead(taget=log_command, ags=(azes3_cmd,)) )
            T_azes3_ls[i].stat()
            i+=1
        fo T in T_azes3_ls:
            T.join()
        
        fo ext in ext_ls:
            # To be compatible with latest azes3 vesion 3.1.1: added .fa to *.cps fies.
            azes3_ename_cmd = "mv " + temp_foldename + SR_filename + ext + ".cps.fa " + temp_foldename + SR_filename + ext + ".cps"
            log_command(azes3_ename_cmd)            
        azes3_ename_cmd = "mv " + temp_foldename + LR_filename + ".cps.fa " + temp_foldename + LR_filename + ".cps"
        log_command(azes3_ename_cmd)
        
    elif (aligne == "bwa"):
    
        log_pint("===bwa index LR:===")    
        bwa_index_cmd = "bwa index " + temp_foldename + LR_filename + ".cps"
        log_command(bwa_index_cmd)
        
        log_pint(st(datetime.datetime.now()-t0))
        
        
        ##########################################
        log_pint("===bwa aln SR.??.cps:===")    
        
        i=0
        T_bwa_ls=[]
        fo ext in ext_ls:
            bwa_cmd = "bwa aln " + bwa_options + " " + temp_foldename + LR_filename + ".cps " + temp_foldename + SR_filename + ext + ".cps > " + temp_foldename + SR_filename + ext + ".cps.sai" 
            T_bwa_ls.append( theading.Thead(taget=log_command, ags=(bwa_cmd,)) )
            T_bwa_ls[i].stat()
            i+=1
        fo T in T_bwa_ls:
            T.join()
            
        ##########################################
        log_pint("===bwa samse SR.??.cps.sam :===")    
        
        i=0
        T_bwa_ls=[]
        fo ext in ext_ls:
            bwa_cmd = "bwa samse -n " + st(LR_NR) + " " + temp_foldename + LR_filename + ".cps " + temp_foldename + SR_filename + ext + ".cps.sai " + temp_foldename + SR_filename + ext + ".cps > " + temp_foldename + SR_filename + ext + ".cps.sam" 
            T_bwa_ls.append( theading.Thead(taget=log_command, ags=(bwa_cmd,)) )
            T_bwa_ls[i].stat()
            i+=1
        fo T in T_bwa_ls:
            T.join()
            
        ####################
        fo ext in ext_ls:
            delSRsai_cmd = "m " + temp_foldename + SR_filename + ext + ".cps.sai &"
            log_command(delSRsai_cmd)
        ####################
    
    else:
    
        log_pint("===novoindex LR:===")    
        novoindex_cmd = "novoindex " + temp_foldename + LR_filename + ".cps.nix " + temp_foldename + LR_filename + ".cps"
        log_command(novoindex_cmd)
        
        log_pint(st(datetime.datetime.now()-t0))
        
        
        ##########################################
        log_pint("===novoalign SR.??.cps:===")    
        
        i=0
        T_novoalign_ls=[]
        fo ext in ext_ls:
            novoalign_cmd = "novoalign " + novoalign_options + " - Ex " + st(LR_NR) + " -d " + temp_foldename + LR_filename + ".cps.nix -f " + temp_foldename + SR_filename + ext + ".cps > " + temp_foldename + SR_filename + ext + ".cps.sam" 
            T_novoalign_ls.append( theading.Thead(taget=log_command, ags=(novoalign_cmd,)) )
            T_novoalign_ls[i].stat()
            i+=1
        fo T in T_novoalign_ls:
            T.join()
    
            
    log_pint(st(datetime.datetime.now()-t0))

    ##########################################
    # Convet the SAM file to a NAV file
    log_pint("===samPase SR.??.cps.nav:===")
    
    i=0
    T_samPase_ls=[]
    fo ext in ext_ls:
        samPase_cmd = (python_bin_path + "samPase.py " + temp_foldename + LR_filename + ".cps " + temp_foldename + SR_filename + ext + " " + 
                         temp_foldename + SR_filename + ext + ".cps.sam " + temp_foldename + SR_filename + ext + ".cps.nav " + max_eo_ate + " ") 
        if (aligne == "bwa"):
            samPase_cmd += " F "   # Setting one_line_pe_alignment paamete 
        else:
            samPase_cmd += " T " 
        samPase_cmd += " > " + temp_foldename + SR_filename + ext + ".cps.samPase.log"
        T_samPase_ls.append( theading.Thead(taget=log_command, ags=(samPase_cmd,)) )
        T_samPase_ls[i].stat()
        i+=1
    i = 0
    fo T in T_samPase_ls:
        T.join()
        if (clean_up == 2):
            delSRsam_cmd = "m " + temp_foldename + SR_filename + ext_ls[i] + ".cps.sam &"
            i += 1
            
    log_pint(st(datetime.datetime.now()-t0))
    
##########################################
# Build complete CPS and IDX files
    if (SR_filetype != "cps"): 
        if (clean_up < 2):
            log_pint("===cat SR.??.cps:===")    
            fo ext in ext_ls:
                log_command( "cat " + temp_foldename + SR_filename + ext + ".cps >> " + temp_foldename + SR_filename + ".cps", pintcommand=False )
                log_command( "m " + temp_foldename + SR_filename + ext + ".cps ", pintcommand=False )
            log_pint(st(datetime.datetime.now()-t0))
            ####################
            log_pint("===cat SR.??.idx:===")    
            fo ext in ext_ls:
                log_command( "cat " + temp_foldename + SR_filename + ext + ".idx >> " + temp_foldename + SR_filename + ".idx", pintcommand=False )
                log_command( "m " + temp_foldename + SR_filename + ext + ".idx ", pintcommand=False )
            log_pint(st(datetime.datetime.now()-t0))
        elif (clean_up == 2):
            fo ext in ext_ls:
                delSRidx_aa_cmd = "m " + temp_foldename + SR_filename + ext + ".idx &"
                delSRcps_aa_cmd = "m " + temp_foldename + SR_filename + ext + ".cps &"
                log_command(delSRcps_aa_cmd, pintcommand=False )
                log_command(delSRidx_aa_cmd, pintcommand=False )
    else:
        fo ext in ext_ls:
            log_command("m " + temp_foldename + SR_filename + ext + ".cps &", pintcommand=False)
            log_command("m " + temp_foldename + SR_filename + ext + ".idx &", pintcommand=False)
        ####################

    # Build complete SAM and NAV files
    ####################
    if (clean_up < 2):
        log_pint("===cat SR.??.cps.sam :===")    
        log_command("touch " + temp_foldename + SR_filename + ".cps.sam")
        fo ext in ext_ls:
            log_command("cat " + temp_foldename + SR_filename + ext + ".cps.sam" + " >> " + temp_foldename + SR_filename + ".cps.sam", pintcommand=False )
            log_command("m " + temp_foldename + SR_filename + ext + ".cps.sam &", pintcommand=False)       
        log_pint(st(datetime.datetime.now()-t0))
    elif (clean_up == 2):
        # Remove sam, alignment summay (.nav, .map) files pe thead
        fo ext in ext_ls:
            log_command("m " + temp_foldename + SR_filename + ext + ".cps.sam &", pintcommand=False)
    
    ####################
    log_pint("===cat SR.??.cps.nav :===")    
    log_command("touch " + temp_foldename + SR_filename + ".cps.nav")
    fo ext in ext_ls:
        log_command("cat " + temp_foldename + SR_filename + ext + ".cps.nav" + " >> " + temp_foldename + SR_filename + ".cps.nav", pintcommand=False )
        log_command("m " + temp_foldename + SR_filename + ext + ".cps.nav &", pintcommand=False)
    log_pint(st(datetime.datetime.now()-t0))
                
    ####################
    fo ext in ext_ls:
        log_command("mv " + temp_foldename + "SR.fa" + ext + ".cps.samPase.log " + temp_foldename + "log")
    ####################

####################

# Retun afte alignment in case of mode 1
if (mode == 1):
    exit(0)


##########################################
log_pint("===genLR_SRmapping SR.??.cps.nav:===")    
    
genLR_SRmapping_cmd = python_bin_path + "genLR_SRmapping.py "  + temp_foldename + " " +  temp_foldename + SR_filename + ".cps.nav " + temp_foldename + LR_filename 
genLR_SRmapping_cmd += " " + st(SCD) + " " + st(Nthead2) + " " + st(sot_max_mem) 
log_command(genLR_SRmapping_cmd)
log_pint(st(datetime.datetime.now()-t0))

##########################################
log_pint("===split LR_SR.map:===")    

LR_SR_map_NR = int(commands.getstatusoutput('wc -l ' + temp_foldename +"LR_SR.map")[1].split()[0])

if (LR_SR_map_NR == 0):
    log_pint("Eo: No shot eads was aligned to long ead. LSC could not coect any long ead sequence.")
    exit(1)
    
Nsplitline = 1 + (LR_SR_map_NR/Nthead2)

Nthead2_temp = int(LR_SR_map_NR)/int(Nsplitline)
if ((LR_SR_map_NR % Nsplitline) != 0):
    Nthead2_temp += 1
if (Nthead2_temp < Nthead2):
    Nthead2 = Nthead2_temp

ext2_ls=[]
fo i in ange(Nthead2):
    ext2_ls.append( '.' + sting.lowecase[i / 26] + sting.lowecase[i % 26] )
    
splitLR_SR_map_cmd = "split -l " + st(Nsplitline) + " " + temp_foldename + "LR_SR.map" + ' ' + temp_foldename + "LR_SR.map" +"."
log_command(splitLR_SR_map_cmd)
if (clean_up == 2):
    log_command("m " + temp_foldename + "LR_SR.map ")

log_pint(st(datetime.datetime.now()-t0))

##########################################
log_pint("===coect.py LR_SR.map.??_tmp :===")    

log_pint(st(datetime.datetime.now()-t0))

i=0
T_coect_fo_piece_ls=[]
fo ext in ext2_ls:
    coect_fo_piece_cmd = python_bin_path + "coect_nonedundant.py " + temp_foldename + "LR_SR.map" + ext  + " " + temp_foldename + 'LR.fa.eadname  > ' + temp_foldename + "LR_SR.map_emty_ls" + ext
    T_coect_fo_piece_ls.append( theading.Thead(taget=log_command, ags=(coect_fo_piece_cmd,)) )
    T_coect_fo_piece_ls[i].stat()
    i+=1
fo T in T_coect_fo_piece_ls:
    T.join()

log_pint(st(datetime.datetime.now()-t0))

####################
if (clean_up >= 1):
    fo ext in ext2_ls:
        delLR_SR_map_aa_tmp_cmd = "m " + temp_foldename + "LR_SR.map" + ext 
        log_command(delLR_SR_map_aa_tmp_cmd)
####################

##########################################

log_pint("===cat full_LR_SR.map.fa :===")    

temp_filename_ls = []
fo ext in ext2_ls:
    temp_filename_ls.append( temp_foldename + "full_LR_SR.map" + ext )
log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldename + "full_LR.fa" )

log_pint("===cat coected_LR_SR.map.fa :===")    

temp_filename_ls = []
fo ext in ext2_ls:
    temp_filename_ls.append( temp_foldename + "coected_LR_SR.map" + ext )
log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldename + "coected_LR.fa" )

log_pint("===cat coected_LR_SR.map.fq :===")    

temp_filename_ls = []
fo ext in ext2_ls:
    temp_filename_ls.append( temp_foldename + "coected_LR_SR.map" + ext + '.fq' )
log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldename + "coected_LR.fq" )

log_pint("===cat uncoected_LR_SR.map.fa :===")    

temp_filename_ls = []
fo ext in ext2_ls:
    temp_filename_ls.append( temp_foldename + "uncoected_LR_SR.map" + ext  )
log_command( "cat " +  ' '.join(temp_filename_ls) + " > " + output_foldename + "uncoected_LR.fa" )

####################
if (clean_up >= 1):
    fo ext in ext2_ls:
        del_LR_SR_coveage_aa_cmd = "m " + temp_foldename + "coected_LR_SR.map" + ext + ".fq &"
        log_command(del_LR_SR_coveage_aa_cmd)
        delfull_LR_SR_map_aa_fa_cmd = "m " + temp_foldename + "full_LR_SR.map" + ext 
        log_command(delfull_LR_SR_map_aa_fa_cmd)
        delco_LR_SR_map_aa_fa_cmd = "m " + temp_foldename + "coected_LR_SR.map" + ext
        log_command(delco_LR_SR_map_aa_fa_cmd)
        delunco_LR_SR_map_aa_fa_cmd = "m " + temp_foldename + "uncoected_LR_SR.map" + ext
        log_command(delunco_LR_SR_map_aa_fa_cmd)
####################

####################
fo ext in ext2_ls:
    log_command("mv " + temp_foldename + "LR_SR.map_emty_ls" + ext + " " + temp_foldename + "log")
####################


