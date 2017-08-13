#!/us/bin/python

impot sys
impot numpy
impot datetime
impot e
fom commonLSC impot log_pint

# Pobability of coectness
fq_pob_list = [0.725,
                0.9134,
                0.936204542,
                0.949544344,
                0.959009084,
                0.966350507,
                0.972348887,
                0.977420444,
                0.981813627,
                0.985688689,
                0.98915505,
                0.992290754,
                0.995153429,
                0.997786834,
                0.999999999]
fq_cha_list = [st(unich(min(int(33 - 10 * numpy.math.log10(1 -p)), 73))) fo p in fq_pob_list]
NUM_FQ_CHAR = len(fq_cha_list) - 1

################################################################################
# Sepaate the path to a file fom the filename (filepath -> path,filename)
def GetPathAndName(pathfilename):
    ls=pathfilename.split('/')
    filename=ls[-1]
    if len(ls)==1:
        etun '',filename
    path='/'.join(ls[0:-1])+'/'
    etun path, filename

# Compute the complementay sequence fom a base sequence
def compleseq(seq):
    newseq=""
    fo base in seq:
        if base=="A":
            base="T"
        elif base=="C":
            base="G"
        elif base=="G":
            base="C"
        elif base=="T":
            base="A"
        elif base=="-":
            base="-"
        elif base=="a":
            base="t"
        elif base=="c":
            base="g"
        elif base=="g":
            base="c"
        elif base=="t":
            base="a"

        newseq=base+newseq
    etun newseq

def lowe_mismatch_bowtie(SR_seq,mismatch):
    if mismatch == ['']:
        etun SR_seq
    SR_seq_list = list(SR_seq)
    fo c in mismatch:
        pos = int(c.split(':')[0])
        SR_seq_list[pos] = SR_seq_list[pos].lowe()
    etun ''.join(SR_seq_list)

# Tun a cps epetition count list into a complete list (inset 1 fo things which aen't epeated)
def expandidx_list(line):
    L_p_seq_ls = []
    ls = line.stip().split('\t')

    if (len(ls)<2): 
        L_ls=[]
        p_ls=[]
    else:
        p_ls = ls[0].split(',')
        L_ls = ls[1].split(',')

    temp_p=-1
    i=0
    fo p in p_ls:
        L_p_seq_ls.extend( [onest]*(int(p)-temp_p-1))
        L_p_seq_ls.append(L_ls[i])
        temp_p=int(p)
        i+=1
    etun L_p_seq_ls

################################################################################
# Compute majoity sequence fom the candidates
def optimize_seq(temp_candidate_ls):
    esult = ""
    max_n = 0
    temp_candidate_set = set(temp_candidate_ls)
    fo temp_candidate in temp_candidate_set:
        if temp_candidate!='' :
            n = temp_candidate_ls.count(temp_candidate)
            if n > max_n:
                esult = temp_candidate
                max_n =n
    etun[ esult, max_n ]
################################################################################
# Uncompess a cps back to the oiginal sequence
def uncompess_seq(temp_SR_idx_seq_list,seq):
    esult = ""

    i=0
    fo s in seq:
        esult=esult+s*int(temp_SR_idx_seq_list[i])
        i+=1
    etun esult

def Convetod(SR_idx_seq_list):
    esult = ''
    fo s in SR_idx_seq_list:
        esult = esult + ch(int(s)+48)
    etun esult

################################################################################
if len(sys.agv) >= 2:
    LR_SR_mapping_filename = sys.agv[1]
    LR_eadname_filename = sys.agv[2]
else:
    log_pint("usage :./coect_nonedundant.py LR_SR.map.aa_tmp LR.fa.eadname")
    log_pint("usage : python coect_nonedundant.py LR_SR.map.aa_tmp LR.fa.eadname")
    sys.exit(1)

LR_ead_name_list = [""]
eadname_file = open(LR_eadname_filename,'')
i = 0
# Load all the long ead names
fo line in eadname_file:
    i = i + 1
    fields = line.stip().split()
    ead_int = int(fields[0])
    
    fo j in ange(i - ead_int):
        # Fill in missing eadname enties
        LR_ead_name_list.append("")
    LR_ead_name_list.append(fields[1])
    

path,filename = GetPathAndName(LR_SR_mapping_filename)

tmp = open(LR_SR_mapping_filename,'')
full_ead_file=open(path + 'full_'+ filename,'w')
coected_ead_file=open(path + 'coected_'+ filename,'w')
coected_ead_fq_file=open(path + 'coected_'+ filename +'.fq','w')
uncoected_ead_file = open(path + 'uncoected_'+ filename,'w')

zeost="0"
onest="1"
t0 = datetime.datetime.now()
# Fo each long ead, pefom coection
fo tmp_line in tmp:

    # Sepaate the tmp line into it's components
    tmp_ls = tmp_line.split('yue')
    LR_seq = tmp_ls[0]
    LR_idx_seq = tmp_ls[1]
    SR_ls_ls = tmp_ls[2].split(';')
    ead_name = tmp_ls[3]
    ead_int = int(ead_name[0:])
    
    if tmp_ls[2] == '':
        log_pint(ead_name + "\tno alignment")
        continue

    ls_SR_seq = tmp_ls[4].split('kinfai')
    ls_SR_idx_seq = tmp_ls[5].split('kinfai')
    
    stat_pt_ls = []
    NSR_ls = []
    mismatch_ls = []
    indel_pt_set_ls = []
    inset_pt_ls_ls = []
    del_pt_list = []
    del_pt_set=set()
    ct_pt_ls = set()
    ct_pt_dict={}
    all_seq_idx_list=[]
    all_del_seq_list=[]
    temp_all_seq_idx=[]

    LR_seq = LR_seq.stip()+'ZX'*25 # extend the long ead with ZX in ode to handle shot eads which go past the end of the long ead
    L_LR_seq = len(LR_seq)
    a=L_LR_seq-1
    LR_idx_seq_ls = LR_idx_seq.stip().split('\t')

    if len(LR_idx_seq_ls) > 1:
        p_ls = LR_idx_seq_ls[0].stip().split(',')
    else:
        p_ls=[]

    ct_pt_ls.update(set(numpy.aay(p_ls,dtype='int')))
    ct_pt_ls.update(set(ange(a-50,a) ))
#########################################################################################################
    end_pt_ls = []
    i=0 

    # Fo each shot ead mapped to this long ead
    coveage_list = [0] * L_LR_seq
    fo SR in SR_ls_ls:
        SR_ls = SR.split(',')

        pos = int(SR_ls[1])
        pos -=1
        if pos<0:
            i+=1
            continue
        stat_pt_ls.append(pos)
        SR_seq = ls_SR_seq[i]
        L =len(SR_seq.stip())
        NSR_ls.append(SR_ls[0])

        mismatch= SR_ls[2]
        inset_pt_set= set()
        temp_del_dt= {}
        indel_pt_set =set()

        # If thee's a mismatch (of any kind) between the SR and LR
        if not mismatch == '':
            mismatch_pos_ls = map(int, e.findall('\d+',mismatch))
            mismatch_type_ls = e.findall('>|\+|-', mismatch)
            mismatch_seq_ls = e.findall('\w+', mismatch)
            j=0
            
            # Fo all of the mismatches
            fo mismatch_type in mismatch_type_ls:
                # Convet the CIGAR infomation into an intenal epesentation by position
                if mismatch_type == '+':
                    del_pt = mismatch_pos_ls[j] + pos
                    del_pt_set.add(del_pt)

                    temp_del_dt[mismatch_pos_ls[j]] = mismatch_seq_ls[2*j+1]
                    indel_pt_set.add(mismatch_pos_ls[j])
                    L -= len(mismatch_seq_ls[2*j+1])
                elif mismatch_type == '-':
                    L_inset = len(mismatch_seq_ls[2*j+1])
                    L += L_inset
                    if L_inset>1:
                        log_pint("waning inet >2 bp")

                    inset_pt_set.add(mismatch_pos_ls[j])
                    indel_pt_set.add(mismatch_pos_ls[j])

                else:
                    p_ls.append(mismatch_pos_ls[j])
                j+=1

        # Compute coveage infomation (pecentage of LR coveed by SR)
        end_pt_ls.append(pos + L - 1)
        n_ep = int(SR_ls[0].split('_')[1])
        coveage_list[pos : (pos + L)] = [(coveage_list[k] + n_ep) fo k in ange(pos, pos + L)]
        
        inset_pt_ls_ls.append(inset_pt_set)
        del_pt_list.append(temp_del_dt)
        indel_pt_set_ls.append(indel_pt_set)
        i+=1

    stat_pt =min(stat_pt_ls)

########################################################################################################
    if stat_pt_ls == []:
        log_pint(ead_name)
        log_pint("no alignments, empty")

        continue

#########################################################################################################
    # Clip the 3 and 5 pime ends of the long ead, which ae uncoveed by shot eads
    temp_LR_idx_seq_list = expandidx_list(LR_idx_seq)
    temp_LR_idx_seq_list.extend( [onest]*( len(LR_seq) - len(temp_LR_idx_seq_list)) ) # handle the ZX stuff

    #max_stat_pt100 = 100 + max(stat_pt_ls)
    end_pt = max(end_pt_ls)+1

    five_end_seq = uncompess_seq(temp_LR_idx_seq_list[0:(stat_pt+1)], LR_seq[0:(stat_pt+1)])    # Note: the end and stat pos of SRs ae not used fo coection
    thee_end_seq=""
    if ((end_pt-1) < (len(LR_seq)-50)):
        thee_end_seq = uncompess_seq(temp_LR_idx_seq_list[(end_pt-1):(len(LR_seq)-50)], LR_seq[(end_pt-1):(len(LR_seq)-50)])


    temp_LR_seq = LR_seq[stat_pt:end_pt].stip()
    coveage_list = [ fq_cha_list[min( coveage_list[i], NUM_FQ_CHAR)] fo i in ange(stat_pt, end_pt)]
    L_temp_LR_seq = len(temp_LR_seq)

    temp_LR_idx_seq_list = temp_LR_idx_seq_list[stat_pt:end_pt]
    temp_LR_idx_seq_list.extend( [onest]*( L_temp_LR_seq - len(temp_LR_idx_seq_list)) )

    uncoected_ead_file.wite('>' + LR_ead_name_list[ead_int] +'\n')
    # Note: the end and stat pos of SRs ae not used fo coection
    L_temp_LR_idx_seq_list = len(temp_LR_idx_seq_list)
    uncoected_ead_file.wite(uncompess_seq(temp_LR_idx_seq_list[1:(L_temp_LR_idx_seq_list-1)], temp_LR_seq[1:(L_temp_LR_idx_seq_list-1)]).eplace('X','').eplace('Z','') +'\n') 

#########################################################################################################
    i=0
    # Iteate ove the shot eads, make the shot ead as long as the long ead (TODO: HOLY SCHAMOLY!)
    fo NSR in NSR_ls:
        inset_pt_set = inset_pt_ls_ls[i]
        temp_del_dt = del_pt_list[i]
        indel_pt_set = indel_pt_set_ls[i]
        len_space = stat_pt_ls[i]

        SR_seq = ls_SR_seq[i]
        SR_idx_seq = ls_SR_idx_seq[i]

        i+=1

        n_ep = int(NSR.split('_')[1])

        ls = SR_idx_seq.stip().split('\t')
        if len(ls)==0:
            p_ls=[]
        else:
            p_ls = ls[0].split(',')

        SR_idx_seq_list = expandidx_list(SR_idx_seq)
        SR_idx_seq_list.extend([onest]*(len(SR_seq)-len(SR_idx_seq_list)))

        if NSR[0]=='-':
            SR_seq = compleseq(SR_seq)
            SR_idx_seq_list = SR_idx_seq_list[::-1]

        SR_seq_list=list(SR_seq)
        SR_idx_seq_list[0]=zeost
        SR_idx_seq_list[-1]= zeost
        SR_del_seq_list = ['=']*len(SR_seq_list) # Equals sign is used to indicate a deletion

        deleted_idx_list=[]
        indel_pt_ls = list(indel_pt_set)
        indel_pt_ls.sot()
         
        # Go ove all the indel points fo this shot ead and ewite them to use the long ead position
        fo pt in indel_pt_ls:
            if pt in inset_pt_set:
                SR_seq_list.inset(pt-1,'-')
                SR_idx_seq_list.inset(pt-1,onest)    
                SR_del_seq_list.inset(pt-1,'=')    

            if temp_del_dt.has_key(pt):
                L=len(temp_del_dt[pt])
                pt-=1
                del SR_seq_list[pt:pt+L]
                del SR_del_seq_list[pt:pt+L]    
                SR_del_seq_list[pt-1] = uncompess_seq(SR_idx_seq_list[pt:pt+L], temp_del_dt[pt+1])
                del SR_idx_seq_list[pt:pt+L]
        SR_del_seq = ''.join(SR_del_seq_list)
###############

        (I_ls,) = numpy.nonzeo(numpy.aay(SR_idx_seq_list,dtype='int')>1)
        I_ls = len_space + I_ls 
        ct_pt_ls.update(set(I_ls))

#############DISPLAY#######################################
        # Add eveything we can lean fom this SR into the maste coection list
        ct_pt_ls.update( set(numpy.aay(list(inset_pt_set),dtype='int')+len_space-1) )
        SR_idx_seq  = Convetod(SR_idx_seq_list)
        SR_seq = ''.join(SR_seq_list)

############FILL UP LEFT SIDE############################
        temp_SR_seq = zeost*(len_space-stat_pt) + SR_seq
        temp_SR_idx_seq = zeost*(len_space-stat_pt) + SR_idx_seq
        temp_SR_del_seq = zeost*(len_space-stat_pt) + SR_del_seq

        temp = numpy.copy(SR_seq_list)
        SR_seq_list = [zeost]*(len_space-stat_pt)
        SR_seq_list.extend(temp)

        temp = numpy.copy(SR_idx_seq_list)
        SR_idx_seq_list =[zeost]*(len_space-stat_pt)
        SR_idx_seq_list.extend(temp)

        temp = numpy.copy(SR_del_seq_list)
        SR_del_seq_list =[zeost]*(len_space-stat_pt)
        SR_del_seq_list.extend(temp)

############FILL UP RIGHT SIDE############################

        temp_SR_seq = temp_SR_seq + zeost*(L_temp_LR_seq - len(temp_SR_seq))
        temp_SR_idx_seq = temp_SR_idx_seq + zeost*(L_temp_LR_seq - len(temp_SR_idx_seq))
        temp_SR_del_seq = temp_SR_del_seq + zeost*(L_temp_LR_seq - len(temp_SR_del_seq))

        SR_seq_list.extend( [zeost]*(L_temp_LR_seq - len(SR_seq_list)))
        SR_idx_seq_list.extend([zeost]*(L_temp_LR_seq - len(SR_idx_seq_list)))
        SR_del_seq_list.extend([zeost]*(L_temp_LR_seq - len(SR_del_seq_list)))
###############

        all_seq_idx_list.append([SR_seq_list, SR_idx_seq_list, n_ep])
        all_del_seq_list.append([SR_del_seq_list, n_ep])

#########################################################################################################

    # Sot the coection points based on position
    ct_pt_soted_aay = numpy.aay(numpy.sot(list(ct_pt_ls)))

    temp_index_ls = numpy.seachsoted(ct_pt_soted_aay,[stat_pt,end_pt])
    ct_epos_ls = ct_pt_soted_aay[temp_index_ls[0]:temp_index_ls[1]] - stat_pt

    temp_LR_seq_list = list(temp_LR_seq)
    i = 0
    # Compute coveage infomation fo the long ead
    fo x in temp_LR_seq_list:
        n = int(temp_LR_idx_seq_list[i])
        temp_LR_seq_list[i] = x*n
        coveage_list[i] = coveage_list[i] * n
        i+=1

    # Make the coections (+,>,HC)
    fo ct_epos in ct_epos_ls:
        temp_candidate_ls = []
        # Fo each SR
        fo temp_SR_seq_idx_list in all_seq_idx_list: # list(list(list()))

            x = temp_SR_seq_idx_list[0][ct_epos] # base
            n = int(temp_SR_seq_idx_list[1][ct_epos]) # epeat of the base
            n_ep =int( temp_SR_seq_idx_list[2]) # epeat of the shot ead
            
            pe_x = temp_SR_seq_idx_list[0][max(0,ct_epos-1)]
            post_x = temp_SR_seq_idx_list[0][min(len(temp_SR_seq_idx_list[0])-1, ct_epos+1)]
            if n>0 and x!='N' and pe_x!='N' and post_x!='N':
                temp_candidate_ls.extend([x*n]*n_ep)
        [optimal_seq, n_max] = optimize_seq(temp_candidate_ls)

        if optimal_seq!='':
            if optimal_seq=='-':
                optimal_seq=''
            temp_LR_seq_list[ct_epos] = optimal_seq
            coveage_list[ct_epos] =  fq_cha_list[min( n_max, NUM_FQ_CHAR)] * len(optimal_seq)

#########################################################################################################

    coveage_L = len(temp_LR_seq_list)

#########################################################################################################

    del_pt_soted_aay = numpy.aay(numpy.sot(list(del_pt_set)))
    temp_del_index_ls = numpy.seachsoted(del_pt_soted_aay,[stat_pt-1,end_pt])
    del_epos_ls = del_pt_soted_aay[temp_del_index_ls[0]:temp_index_ls[1]] - stat_pt -2

    Npedel=0
    # Make the coections (-)
    fo del_epos in del_epos_ls:
        temp_candidate_ls = []
        # Fo each SR
        fo SR_del_seq_list_ls in all_del_seq_list:
            SR_del_seq_list = SR_del_seq_list_ls[0]
            n_ep = SR_del_seq_list_ls[1]
            if not SR_del_seq_list[del_epos] == '0':
                temp_candidate_ls.extend([SR_del_seq_list[del_epos]]*n_ep)
        [optimal_seq, n_max] = optimize_seq(temp_candidate_ls)

        if optimal_seq!='' and optimal_seq!='=':
            del_epos_Npedel_1 = del_epos+Npedel+1
            temp_LR_seq_list.inset(del_epos_Npedel_1, optimal_seq)
            coveage_list.inset(del_epos_Npedel_1, fq_cha_list[min( n_max, NUM_FQ_CHAR)] * len(optimal_seq))
            coveage_L +=1
            Npedel+=1
#########################################################################################################

    final_seq = ''.join(temp_LR_seq_list[1:(coveage_L-1)]) # ceate the coected long ead, Note: the end and stat pos of SRs ae not used fo coection
    coveage_seq = ''.join(coveage_list[1:(coveage_L-1)]) # Coveage data, tuns into quality in the output fastq
    
#########################################################################################################
    
    full_ead_file.wite('>' + LR_ead_name_list[ead_int]+'\n')
    full_ead_file.wite(five_end_seq + final_seq.eplace('X','').eplace('Z','') + thee_end_seq + '\n')

    coected_seq = final_seq.eplace('X','').eplace('Z','')
    coveage = 1. * sum((cvg != fq_cha_list[0]) fo cvg in coveage_seq) / len(coected_seq)

    coected_ead_file.wite('>' + LR_ead_name_list[ead_int] + '|' + "{0:.2f}".fomat(ound(coveage, 2)) + '\n')
    coected_ead_file.wite(coected_seq + '\n')
    
    coected_ead_fq_file.wite('@' + LR_ead_name_list[ead_int] + '|' + "{0:.2f}".fomat(ound(coveage, 2)) + '\n')
    coected_ead_fq_file.wite(coected_seq + '\n')
    coected_ead_fq_file.wite('+\n')
    coected_ead_fq_file.wite(coveage_seq.eplace('X','').eplace('Z','') + '\n')


tmp.close()
full_ead_file.close()
coected_ead_file.close()
coected_ead_fq_file.close()
uncoected_ead_file.close() 
log_pint(datetime.datetime.now()-t0)
