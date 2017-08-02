import pandas as pd
from tkinter import *
import os
import saveOutput

#Runs SeqProp for multiple sequences.
def multSeqProp(sequences, restriction_sites):
    sequences=sequences.split()
    index=1
    for seq in sequences:
        SeqProp(seq, restriction_sites,index)
        index=index+1

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)

def SeqProp(sequence, restriction_sites, index):
    # Add re_list+re_Sites?, save_path_file, 
    
    '''Calculates GC_content, Tm, start/stop indices, possible exons, re_sites, motifs, binding domains, 
    promoters/enahncers, etc'''
    
    # GC Content
    GC_num = 0
    for nuc in sequence:
        # Finding number of G's and C's
        if nuc is 'G' or nuc is 'C':
            GC_num += 1
    GC_content = str(round(float(GC_num)/float(len(sequence))*100, 2))+'%'
    
    # Approximate Melting Temperature based on Wallace Rule
        #only good for short oligos
    if len(sequence) <= 14:
        Tm = 2*(len(sequence)-GC_num) + 4*GC_num
    else:
        Tm= 64.9 +41*(GC_num-16.4)/len(sequence)
    Tm = str(round(Tm, 2))+ ' C'
    # Add one for long oligos
    
    # Reverse-Complement
    try:
        complement = ''.join([{'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a', '\n':'', ' ':''}[B] for B in list(sequence)])
    except:
        assert False, "Make sure only A, T, C, or G!"
    reverse_complement = complement[::-1]
    
    # Indices of possible start and stop codons
    if 'A' in sequence:
        start_idxs = list(find_all(sequence, 'ATG'))
        if len(start_idxs) == 0:
            start_idxs.append('None')
        # Stop Codons UAA UAG UGA
        UAA_idxs = list(find_all(sequence, 'TAA'))
        UAG_idxs = list(find_all(sequence, 'TAG'))
        UGA_idxs = list(find_all(sequence, 'TGA'))
        stop_idxs = UAA_idxs + UAG_idxs + UGA_idxs
        stop_idxs.sort()
        if len(stop_idxs) == 0:
            stop_idxs.append('None')
    else:
        start_idxs = list(find_all(sequence, 'atg'))
        if len(start_idxs) == 0:
            start_idxs.append('None')
        # Stop Codons UAA UAG UGA
        UAA_idxs = list(find_all(sequence, 'taa'))
        UAG_idxs = list(find_all(sequence, 'tag'))
        UGA_idxs = list(find_all(sequence, 'tga'))
        stop_idxs = UAA_idxs + UAG_idxs + UGA_idxs
        stop_idxs.sort()
        if len(stop_idxs) == 0:
            stop_idxs.append('None')
        
    # Indices of possible exons
    exons = []
    i = 0
    for start_codon in start_idxs:
        exons.append([start_codon, stop_idxs[i]])
    if len(exons) == 0:
        exons.append('None')
        
    # Restriction enzyme sites --> COPY AND PASTE FROM GEL.VIZ >> APPEND RE_NAMES TOO!
#    restriction_sites = pd.read_csv('ISTools/ISToolkit/data/restriction_sites2.csv', sep=',')
    restriction_enzymes = restriction_sites.values.T.tolist()
    restriction_enzymes = restriction_enzymes[0]
    restriction_sites  = restriction_sites.values.tolist()
    restriction_sites = restriction_sites[1:]
    restriction_seqs = []
    for site in restriction_sites:
        clean_site = [x for x in site if str(x) != 'nan']
        del clean_site[0]
        restriction_seqs.append(clean_site)
    del restriction_enzymes[0]
    # finding idxs
    cut_idxs = []
    for re in restriction_enzymes:
        re_index = restriction_enzymes.index(re)
        seq_index = 0
        re_idxs = []
        re_idxs.append(str(re))
        for cut_site in restriction_seqs[re_index]:
            site_idx = sequence.find(restriction_seqs[re_index][seq_index])
            if site_idx != -1:
                re_idxs.append(str(site_idx))
                while site_idx < len(sequence):
                    site_idx = sequence.find(restriction_seqs[re_index][seq_index], site_idx)
                    if site_idx == -1:
                        break
                    re_idxs.append(str(site_idx))
                    site_idx += len(restriction_seqs[re_index][seq_index])
            seq_index += 1
        if len(re_idxs) >= 2:
            re_idxs = ', '.join(re_idxs) + ' '
            cut_idxs.append(re_idxs)
    if len(cut_idxs) == 0:
        cut_idxs.append('None')
    cut_idxs.sort()
    
    # Cleaning NaNs
    restriction_seqs = []
    for site in restriction_sites:
        clean_site = [x for x in site if str(x) != 'nan']
        del clean_site[0]
        restriction_seqs.append(clean_site)
    del restriction_enzymes[0]    
    
    # Motifs
    
    # Binding Domains
    
    # Promoter and Enhancer regions
    
    ''' WRITING FILE '''
    lines = []
    
    lines.extend(('GC Content: ' + GC_content, ' '))
    lines.extend(('Melting Temperature (C): ' + Tm, ' '))
    lines.extend(('Reverse Complement: ' + reverse_complement, ' '))
    lines.extend(('Start Codon Indices: ' + ', '.join(str(start_idxs)), ' '))
    lines.extend(('Stop Codon Indices: ' + ', '.join(str(stop_idxs)), ' '))
    lines.extend(('Possible Exons :', ' ')) # Need to fix, need to loop?
    lines.extend(('Restriction Cut Sites :' + ' ' + ', '.join(cut_idxs), ' ')) # Need to add names of enzymes

    out_text=''
    for item in lines:
        out_text=out_text+("%s\n" % item)
        
    saveOutput.saveData(out_text,"Seqprop Analysis "+str(index))