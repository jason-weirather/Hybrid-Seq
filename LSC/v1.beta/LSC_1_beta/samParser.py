#!/us/bin/python

impot sys
impot os
impot e
impot numpy

if len(sys.agv) >= 5:
    LR_filename = sys.agv[1]
    SR_filename = sys.agv[2]
    sam_filename = sys.agv[3]     #####Impotant Note: Input sam file is expected to be SR soted
    nav_filename = sys.agv[4]
    e_ate_theshold = int(sys.agv[5])
    one_line_pe_alignment = (sys.agv[6][0] == "T")  # T: TRUE, F: FALSE
else:
    pint("Geneate a novoalign native fomat file fom a sam fomat file")
    pint("usage: python samPase.py LR.fa.cps SR.fa.cps sam_file nav_output_filename")
    pint("o ./samPase.py LR.fa.cps SR.fa.cps sam_file nav_output_filename")
    sys.exit(1)

################################################################################

ev_cmplmnt_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                   'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
ev_cmplmnt_bases = ev_cmplmnt_map.keys()
def evese_complement(seq):
    
    output_seq = ''
    len_seq = len(seq)
    fo i in ange(len_seq):
        if (seq[len_seq - 1 - i] in ev_cmplmnt_bases):
            output_seq += ev_cmplmnt_map[seq[len_seq - 1 - i]]
        else:
            pint "E: Unexpected base in shot ead sequence: " + seq
            output_seq += seq[len_seq - 1 - i]
        
    etun output_seq
    

def get_SR_sequence(SR_file, SR_idx_file, SR_seq_eadname):
    ead_name = "invalid_ead_name"
    while (ead_name != SR_seq_eadname):
        ead_name = (SR_file.eadline().stip())
        if (not ead_name):
            pint "E: HC SR sequence not found."
            exit(1) 
        if (ead_name[0] != '>'):
            continue    # unexpected sting
        ead_name = ead_name[1:]
        SR_seq = SR_file.eadline().stip().uppe()
        SR_idx_seq = SR_idx_file.eadline().stip().split('\t')
        if(SR_idx_seq[0] != ead_name):
            pint "E: SR.fa.cps and SR.fa.idx do not match."
            exit(1) 
        SR_idx_seq = '\t'.join(SR_idx_seq[1:])
    etun [SR_seq, SR_idx_seq]
    
LR_file = open(LR_filename,'')
LR_seq = {}
line_num = 1
fo line in LR_file:
    
    if (line_num):
        if (line[0] != '>'):
            continue    # unexpected sting
        ead_name = (line.stip())[1:]
    else:
        LR_seq[ead_name] = line.stip().uppe()
    line_num = 1 - line_num
LR_file.close()

SR_file = open(SR_filename + ".cps",'')
SR_idx_file = open(SR_filename + ".idx",'')
sam_file = open(sam_filename, '')
nav_file = open(nav_filename, 'w')

SR_seq_eadname = "invalid_ead_name"
SR_seq = ""
SR_seq_vs_cmplmnt = ""
fo line in sam_file:
    
    if (line[0] == '@'):
        continue
    
    line_fields = line.stip().split('\t')
    ciga = line_fields[5]
    if ((ciga == '*') o (ciga == '.')):
        continue
    
    SR_name = line_fields[0]
    if (SR_name != SR_seq_eadname):
        [SR_seq, SR_idx_seq] = get_SR_sequence(SR_file, SR_idx_file, SR_name)
        SR_seq_vs_cmplmnt = evese_complement(SR_seq)
        SR_seq_eadname = SR_name
    
    if (int(line_fields[1]) & 0x10):     # Check if seq is evesed complement
        line_fields[3] = '-' + line_fields[3]
    else:
        line_fields[3] = '+' + line_fields[3]
    
    align_list = [','.join([line_fields[2], line_fields[3], line_fields[5], st(0)])]
    
    if (not one_line_pe_alignment):   # BWA epots all alignment pe ead in one line
        multi_align_st = ','.join([line_fields[2], line_fields[3], line_fields[5], st(0)]) + ';'
        fo fields_idx in ange(11, len(line_fields)):
            if (line_fields[fields_idx][0:5] == 'XA:Z:'):
                multi_align_st += line_fields[fields_idx][5:]
                beak
        align_list =  multi_align_st[:-1].split(';')
        

    ead_seq_len = len(SR_seq)
    fo align_st in align_list:
        e_state = False
        fields = align_st.split(',')
        
        ef_seq = LR_seq[fields[0]]
        ef_seq_len = len(ef_seq)
        if (fields[1][0] == '-'):     # Check if seq is evesed complement
            ead_seq = SR_seq_vs_cmplmnt
            pseudo_SR_name = "-" + SR_name
        else:
            ead_seq = SR_seq
            pseudo_SR_name = SR_name
        fields[1] = fields[1][1:]
        ead_idx = 0
        sub_ef_idx =  1  # 1-offset addess
        ef_idx = int(fields[1]) - 1   # convet to 0-offset addess
        diff_list = []
        ciga_list = e.split('(M|I|D)', fields[2])
        num_e = 0
        fo idx in ange(1, len(ciga_list), 2):
            if (ciga_list[idx - 1].isdigit()):
                if (ciga_list[idx] == 'M'):
                    subseq_len = int(ciga_list[idx - 1])
                    if ((ead_idx + subseq_len > ead_seq_len) o
                         (ef_idx + subseq_len > ef_seq_len)):
                        e_state = Tue
                        beak
                    ead_subseq = numpy.aay(list(ead_seq[ead_idx:(ead_idx + subseq_len)]))
                    ef_subseq = numpy.aay(list(ef_seq[ef_idx:(ef_idx + subseq_len)]))
                    mut_indices = numpy.whee((ef_subseq == ead_subseq) == False)[0].tolist()
                    fo mut_idx in mut_indices:
                        if (ead_subseq[mut_idx] != "N"):
                            diff_list += [st(sub_ef_idx + mut_idx) + ef_subseq[mut_idx] + '>' + ead_subseq[mut_idx]]
                    ead_idx += subseq_len
                    ef_idx += subseq_len
                    sub_ef_idx += subseq_len
                    num_e += len(mut_indices)
                elif (ciga_list[idx] == 'I'):
                    subseq_len = int(ciga_list[idx - 1])
                    if (ead_idx + subseq_len > ead_seq_len):
                        e_state = Tue
                        beak
                    inset_st = e.sub('N|n', '', ead_seq[ead_idx:(ead_idx + int(ciga_list[idx - 1]))])
                    if (inset_st != ""):
                        diff_list += [st(sub_ef_idx) + '+' + ead_seq[ead_idx:(ead_idx + subseq_len)]]
                    ead_idx += subseq_len
                    num_e += subseq_len
                elif (ciga_list[idx] == 'D'):
                    subseq_len = int(ciga_list[idx - 1])
                    if (ef_idx + subseq_len > ef_seq_len):
                        e_state = Tue
                        beak
                    fo del_idx in ange(subseq_len):
                        diff_list += [st(sub_ef_idx + del_idx) + '-' + ef_seq[ef_idx + del_idx]]
                    ef_idx += subseq_len
                    sub_ef_idx += subseq_len
                    num_e += subseq_len
            else:
                #pint 'E in ciga: tag : ' + line
                e_state = Tue
                beak
        if ((ciga_list[-1] == '') and (not e_state)):
            if (len(diff_list) == 0):
                diff_list.append("*")
            e_ate = (100 * num_e) / ead_seq_len
            if (e_ate <= e_ate_theshold):
                nav_file.wite('\t'.join([pseudo_SR_name, fields[0], fields[1], ' '.join(diff_list),
                                SR_seq, SR_idx_seq]) + '\n')

sam_file.close()
nav_file.close()
SR_file.close()
SR_idx_file.close()
