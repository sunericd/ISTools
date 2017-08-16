import pandas as pd
from tkinter import *
import os
import mods.saveOutput as saveOutput

#Runs SeqProp for multiple sequences.
def multSeqProp(sequences, restriction_sites):
    sequences=sequences.split()
    index=1
    for seq in sequences:
        SeqProp(seq, restriction_sites,index)
        index += 1

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)

def SeqProp(sequence, restriction_sites, index):    
    '''Calculates GC_content, Tm, start/stop indices, possible exons, re_sites, motifs, binding domains, 
    promoters/enahncers, etc'''
    # GC Content
    GC_content, GC_num = getGC(sequence)
    # Melting Temp
    Tm = getTm(sequence, GC_num)
    # Reverse-Complement
    reverse_complement = getRevComp(sequence)
    # Indices of possible start and stop codons
    start_idxs, stop_idxs = getStartStop(sequence)
    # Indices of possible exons
    exons = getExons(sequence, start_idxs, stop_idxs)
    # Restriction enzyme sites 
    cut_idxs = getRe(sequence, restriction_sites)   
    
    # Motifs
    
    # Binding Domains
    
    # Promoter and Enhancer regions
    
    ''' WRITING FILE LINES '''
    out_text = getListFile(GC_content, Tm, reverse_complement, start_idxs, stop_idxs, exons, cut_idxs)
        
    saveOutput.saveData(out_text,"Seqprop Analysis "+str(index))


def getGC(sequence):
    GC_num = 0
    for nuc in sequence:
        if nuc.upper() not in 'ATCG':
            raise Exception('Please only input DNA sequences: ATCG!')
        # Finding number of G's and C's
        if nuc is 'G' or nuc is 'C':
            GC_num += 1
    GC_content = round(float(GC_num)/float(len(sequence))*100, 2)
    return (GC_content, GC_num)

def getTm(sequence, GC_num):
    # Approximate Melting Temperature based on Wallace Rule
    #only good for short oligos
    if len(sequence) <= 14:
        Tm = 2*(len(sequence)-GC_num) + 4*GC_num
    else:
        Tm= 64.9 +41*(GC_num-16.4)/len(sequence)
    Tm = round(Tm, 2)
    return (Tm)

def getRevComp(sequence):
    # gets reverse complement of sequence
    try:
        complement = ''.join([{'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a', '\n':'', ' ':''}[B] for B in list(sequence)])
    except:
        assert False, "Make sure only A, T, C, or G!"
    reverse_complement = complement[::-1]
    return (reverse_complement)

def getStartStop(sequence):
    if 'A' in sequence:
        start_idxs = list(find_all(sequence, 'ATG'))
        if len(start_idxs) == 0:
            start_idxs.append('None')
            print ("Warning: No start codons detected.")
        # Stop Codons UAA UAG UGA
        UAA_idxs = list(find_all(sequence, 'TAA'))
        UAG_idxs = list(find_all(sequence, 'TAG'))
        UGA_idxs = list(find_all(sequence, 'TGA'))
        stop_idxs = UAA_idxs + UAG_idxs + UGA_idxs
        stop_idxs.sort()
        if len(stop_idxs) == 0:
            stop_idxs.append('None')
            print ("Warning: No stop codons detected.")
    else:
        start_idxs = list(find_all(sequence, 'atg'))
        if len(start_idxs) == 0:
            start_idxs.append('None')
            print ("Warning: No start codons detected.")
        # Stop Codons UAA UAG UGA
        UAA_idxs = list(find_all(sequence, 'taa'))
        UAG_idxs = list(find_all(sequence, 'tag'))
        UGA_idxs = list(find_all(sequence, 'tga'))
        stop_idxs = UAA_idxs + UAG_idxs + UGA_idxs
        stop_idxs.sort()
        if len(stop_idxs) == 0:
            stop_idxs.append('None')
            print ("Warning: No stop codons detected.")
    return (list(start_idxs), list(stop_idxs))

def getExons(sequence, start_idxs, stop_idxs):
    exons = []
    for start_codon in start_idxs:
        for stop_codon in stop_idxs:
            if stop_codon - start_codon >= 3:
                exons.append([start_codon, stop_codon])
    if len(exons) == 0:
        exons.append('None')
        print ("Warning: No exons detected.")
    return(exons)

def getRe(sequence, restriction_sites):
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
        print ("Warning: No restriction sites detected.")
    cut_idxs.sort()
    # Cleaning NaNs
    # restriction_seqs = []
    # for site in restriction_sites:
    #     clean_site = [x for x in site if str(x) != 'nan']
    #     del clean_site[0]
    #     restriction_seqs.append(clean_site)
    # del restriction_enzymes[0] 
    return (cut_idxs)


def getListFile(gc, tm, rc, strts, stps, exns, res):

    lines = []

    lines.extend(('GC Content: ' + str(gc) + '%', ' '))
    lines.extend(('Melting Temperature (C): ' + str(tm), ' '))
    lines.extend(('Reverse Complement: ' + rc, ' '))
    lines.extend(('Start Codon Indices: ' + ''.join(str(strts)), ' '))
    lines.extend(('Stop Codon Indices: ' + ''.join(str(stps)), ' '))
    lines.extend(('Possible Exons :'+ ' ' + ''.join(str(exns)), ' '))
    lines.extend(('Restriction Cut Sites :' + ' ' + ''.join(str(res)), ' ')) # Need to add names of enzymes

    out_text=''
    for item in lines:
        out_text=out_text+("%s\n" % item)

    return (out_text)