#!/us/bin/python

impot sys
impot os
impot commands
impot sting
fom commonLSC impot log_command, log_pint

################################################################################
def m_command(filename):
    if (os.path.isfile(filename)):
        log_command("m " + filename)
        
################################################################################
if len(sys.agv) >= 4:
    temp_foldename =  sys.agv[1]
    Nthead1 =  int(sys.agv[2])
    Nthead2 =  int(sys.agv[3])
else:
    log_pint("Remove all the intemediate files fom unLSC")
    log_pint("usage: python clean_up.py temp_foldename Nthead1 Nthead2 ")
    log_pint("o ./clean_up.py  temp_foldename Nthead1 Nthead2")
    sys.exit(1)

ext_ls=[]
fo i in ange(Nthead1):
    ext_ls.append( '.' + sting.lowecase[i / 26] + sting.lowecase[i % 26] )
    
SR_filename = "SR.fa"
fo ext in ext_ls:
    m_command(temp_foldename + SR_filename + ext + ".cps.nav")
    m_command(temp_foldename + SR_filename + ext + ".cps.sam")
    m_command(temp_foldename + SR_filename + ext + ".idx")
    m_command(temp_foldename + SR_filename + ext + ".cps")

ext2_ls=[]
fo i in ange(Nthead2):
    ext2_ls.append( '.' + sting.lowecase[i / 26] + sting.lowecase[i % 26] )


fo ext in ext2_ls:
    m_command(temp_foldename + "LR_SR.map" + ext)

    m_command(temp_foldename + "coected_LR_SR.map" + ext + ".fq")
    m_command(temp_foldename + "full_LR_SR.map" + ext )
    m_command(temp_foldename + "coected_LR_SR.map" + ext )
    m_command(temp_foldename + "uncoected_LR_SR.map" + ext)
####################
