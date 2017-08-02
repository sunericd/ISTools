import csv
from tkinter import *
import os
import saveOutput

def msm(filename, motif_length = 10, match_threshold = 3):
    '''
    filename = STRING (address of the tsv/csv file to be read; organized with 1. NAME, 2. SEQ)
    motif_length = INT (min length of the matching seqeunces/motifs)
    match_threshold = INT (number of genes needed for a matching sequence to be returned)
    '''
    
    # Type Errors
    try:
        motif_length  = int(motif_length)
    except:
        raise Exception('Motif Length should be an integer!')
    
    try:
        match_threshold  = int(match_threshold)
    except:
        raise Exception('Match Threshold should be an integer!')
    
    
    gene_list = []
    
    # Reading the file: Handles .tsv and .csv files
    with open(filename) as file:
        if '.tsv' in filename:
            for line in csv.reader(file, delimiter="\t"):
                gene_list.append(line)
        elif '.csv' in filename:
            for line in csv.reader(file, delimiter="\t"):
                gene_list.append(line)
        else:
            raise Exception('Please use a tsv or csv file as input!')
                
    sequences_list = []
    name_list = []
    
    if len(gene_list[1]) != 2:
        raise Exception('For each line (gene), the first entry should be the NAME and the second entry should be its SEQEUNCE.')

    # Collects the name of the gene (column 1) and the sequence (Column 2)
    for gene in gene_list:  
        gene_info = str(gene[0])
        name_list.append(gene_info)
        seq = gene[1]
        sequences_list.append(seq)
        
    length_list = []
    
    for sequence in sequences_list:
        length_list.append(len(sequence))
    max_length = max(length_list)
    
    x = 0
    y = motif_length
    
    sequence_dict = {}
    
    # Writing a dictionary for matches
    for n in range(max_length - y):
        for sequence in sequences_list:
            if y > len(sequence):
                continue
            elif not sequence[x:y] in sequence_dict:
                order = sequences_list.index(sequence)
                sequence_dict[sequence[x:y]] = []
                sequence_dict[sequence[x:y]].append(name_list[order])
            else:
                if name_list[order] not in sequence_dict[sequence[x:y]]:             
                    sequence_dict[sequence[x:y]].append(name_list[order])                
        x += 1
        y += 1
        
    new_seq_dict = {}
    # Cleaning for only matches that exceed the match_threshold
    for sequence,gene_name in sequence_dict.items():
        num_genes = len(gene_name)
        if num_genes >= match_threshold:
            new_seq_dict[sequence] = gene_name

    out_text=''
    for key, value in new_seq_dict.items():
        out_text=out_text+('%s:%s\n' % (key, value))

    saveOutput.saveData(out_text,"MSM")