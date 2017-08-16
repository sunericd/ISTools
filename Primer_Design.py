
# coding: utf-8

# In[411]:

import numpy as np


# In[499]:

#########################Select type of primer to design
######## 1 = primer pair
######## 2 = forward primer
######## 3 = reverse primer


####Use primer design function. inputs are primer_type, sequence, start index, end, primer restrictions
#####Input 0 for start and end indices for primer types 2 and 3
#primer_design(primer_type, sequence, started, ended, primer_length, temp_range, gc_range)


##########################################Test Sequences
######################for primer pair
#seq1 = 'ATTTCTGAGGAAAGAGGACTATACCCATTAGGAAACGAATTGCCCGAGTAGTCTCCTCTGCCGACTTAAACCAACCTTTTTCTATTTCTCTTTTCTTTTCTCCCTCTTTTTTCTCTGTACTAGCATCCAAAAGCAAGCATCCATCCGAGTCCCAGTCGCAATCTCACATCTCCAATTTAACGTATCCATTGCATTTCCTCATTCGGTTTAACTCCTCTGCATTTCTTTTCTGACCCATAGCATTTCTTACATTCCATTGCATCTCCCTTTTACTCTCGTTCAAGACACTGATTTGATACGCTTTCTGTACGATGGCCATATTGAAGGATACCATAATTAGATACGCTAATGCAAGGTATGCTACCGCTAGTGGCACTTCCACCGCCACTGCCGCCTCTGTCAGCGCTGCCTCATGTCCTAATTTGCCCTTGCTCTTGCA'
#started = 30
#ended = 411
######################for forward
#seqf = 'TGGATTTCTGAGGAAAGAGGACTATACCCA'
#started = 0
#ended = 0
######################for reverse
#seqr = 'TGCCTCATGTCCTAATTTGCCCTTGCTCTTGCA'
#started = 0
#ended = 0

#sequence = seq1
#primer_design(1, sequence, started, ended, primer_length, temp_range, gc_range)


# In[490]:

################################## Forward and Reverse Sequence Function
def split_seq(split_sequence, started, ended):
    #Error Handling
    if int(started) > int(ended):
        print("Error: start index is after than end index. Reorder indices so start is inputted before end index.")
    elif int(started) > len(split_sequence):
        print("Error: Start index is out of range of sequence input. Choose smaller start index.")
    elif int(ended) > len(split_sequence):
        print("Error: End index is out of range of sequence input. Choose smaller end index.")
    elif int(ended) == int(started):
        print("Error: Start and end indices are equal. No sequence is amplified.")
    else:
        fwd_seq = seq[0:started]
        rv_seq = seq[ended-2:len(seq)]
        return fwd_seq, rv_seq
    #function returns fwd_seq and rv_seq


# In[456]:

#################################Primer Find Function 
def primerfind(sequence, length = [20, 19, 21, 18, 22], temperature = [40, 60], gcontent = [.4, .6]):
    
    #error handling for sequence inputs that are too short 
    if len(sequence) < min(length):
        print("Error: Seqeunce input is shorter than primer. Input longer sequence for primer search.")
    
    else:
        #Initialize
        primer=[] #primer sequences
        gc=[] #gc content
        temp=[] #melting temmperature
        index=[] #index of the beginning of primer

        #identify start index
        start = len(sequence) - min(length) + 1

        #Flip through the sequence to find primers
        for i in range(start-1, 2, -1):

            #Separate into appropriate bp length sequences 
            for j in length: #flips through all possible ideal lengths 
                seq = sequence[i:(i + j-1)] #identify sequences
                GC_num = 0 #GC content 

                #Calculate GC% content and Tm
                #Approximate Melting Temperature based on Wallace Rule
                #only good for short oligos
                for nuc in seq:
                    if nuc not in "ATCGatcg":
                        raise Exception('Sequences must only contain A, a, T, t, C, c, or G, g!')
                    if nuc is 'G' or nuc is 'C' or nuc is 'g' or nuc is 'c': 
                        GC_num += 1
                        GC = round(float(GC_num)/float(len(seq)), 2)
                        Tm= 64.9 +41*(GC_num-16.4)/len(seq)
                        Tm=(round(Tm, 2))        

                #Filter out Tm and Content if not in range
                if Tm >= min(temperature) and Tm <= max(temperature) and GC >= min(gcontent) and GC <= max(gcontent):
                    primer.append(seq)
                    gc.append(GC)
                    temp.append(Tm)
                    index.append(i)
                    break
            
        return primer, gc, temp, index
########### End of primer find function


# In[445]:

#####################################Reverse Complement Function
def rv_comp(unrev, ended):
    
    rv_complement = []
    rv_gc = unrev[1] 
    rv_temp = unrev[2]
    rv_index = np.add(unrev[3], ended)
    
    for primer in unrev[0]:
        try:
            complement = ''.join([{'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a', '\n':'', ' ':''}[B] for B in list(primer)])
        except:
            assert False, "Make sure only A, T, C, or G!"
        rv_complement.append(complement[::-1])
    
    #reverse list in all elements
    rv_complement.reverse()    
    rv_gc = rv_gc[::-1]
    rv_temp = rv_temp[::-1]
    rv_index = rv_index[::-1]
    
    #create reverse final list
    rv_final = []
    rv_final.append(rv_complement)
    rv_final.append(rv_gc)
    rv_final.append(rv_temp)
    rv_final.append(rv_index)
    return rv_final
#########Returns reverse complement, gc, temp, and index of sequences    


# In[396]:

###########################################Forward and Reverse Match Function
#Tm should be +- 2 degrees celsius difference
#optimize for the shortest sequence amplified

def primer_pair(fwd_final, rv_final, started, ended):
    
    fwd_primers = fwd_final[0]
    fwd_gc = fwd_final[1]
    fwd_temp = fwd_final[2]
    fwd_index = fwd_final[3]
    
    rv_primers = rv_final[0]
    rv_gc = rv_final[1]
    rv_temp = rv_final[2]
    rv_index = rv_final[3]
    
    primer_pair = []

    for i in range(len(fwd_final[0])):
        for j in range(15): 
    #primer can only be paired with the reverse primer that is at most 15bp + amplified length. 
    #if nothing fits, it moves to the next primer and pairs sequences
        
            diff = abs(fwd_temp[i] - rv_temp[j])
            if diff<3:
                seq_length = rv_index[j] - fwd_index[i] 
                things = [fwd_primers[i], fwd_index[i], fwd_temp[i], rv_primers[j], rv_index[j], rv_temp[j], seq_length]                    
                primer_pair.append(things)
                break
                #if tm of fwd 20bp primer and first reverse primer checked is good, then it moves on the the second primer
                #this is so that not all sequences have similar fwd primers with only a few bp difference

        if len(primer_pair)>=10:
            #display the first 6 primer pair along with lengths and then stop. 
            break
            
    return primer_pair
#primer_pair outputs include fwd primer info, rv primer info, and length of sequence


# In[500]:

##############################################Primer Design Function
def primer_design(primer_type, sequence_input, started, ended, primer_length = [20, 19, 21, 18, 22], temp_range= [4, 60], gc_range = [.4, .6]):
    if primer_type == 2:
        forward = primerfind(sequence_input, primer_length, temp_range, gc_range)
        if len(forward[0]) == 0:
            print("Error: No forward primers found. Extend input sequence upstream for primer search.")
        else:
            print("Column 1: Forward primer sequence. "
                  "Column 2: GC content. "
                  "Column 3: Melting temprature. "
                  "Column 4: start index of primer.")
            print(forward)

    elif primer_type == 3:
        rev = primerfind(seqr, primer_length, temp_range, gc_range)
        reverse = rv_comp(rev, ended)
        if len(reverse[0]) == 0:
            print("Error: No forward primers found. Extend input sequence upstream for primer search.")
        else:
            print("Column 1: Reverse primer sequence. "
                  "Column 2: GC content. "
                  "Column 3: Melting temprature. "
                  "Column 4: start index of primer.")
            rv_index = reverse[3]

            rv_final = []
            rv_final.append(reverse[0])
            rv_final.append(reverse[1])
            rv_final.append(reverse[2])
            rv_final.append(rv_index)
            print(rv_final)
    
    elif primer_type == 1:
        splitted = split_seq(sequence_input, started, ended)
        fwd_seq = splitted[0]
        rv_seq = splitted[1]
        
        #Error handling if sequence unable to be split
        if len(fwd_seq) < 18:
            print("Error: Sequence cannot be split to identify forward primers. Extend sequence upstream of amplification and check start index.")
            if len(rv_seq) > 18:
                rev = primerfind(sequence_input, primer_length, temp_range, gc_range)
                reverse = rv_comp(rev, ended)
                if len(reverse[0]) == 0:
                    print("Error: No reverse primers found. Extend input sequence donstream for primer search.")
                else:
                    print("Column 1: Reverse primer sequence. "
                          "Column 2: GC content. "
                          "Column 3: Melting temprature. "
                          "Column 4: start index of primer.")
                    rv_index = reverse[3]
                    rv_final = []
                    rv_final.append(reverse[0])
                    rv_final.append(reverse[1])
                    rv_final.append(reverse[2])
                    rv_final.append(rv_index)
                    print(rv_final)

        elif len(rv_seq) < 18:
            print("Error: Sequence cannot be split to identify reverse sequences. Check reverse indices.")
            forward = primerfind(sequence_input, primer_length, temp_range, gc_range)
            print("Column 1: Forward primer sequence. "
                  "Column 2: GC content. "
                  "Column 3: Melting temprature. "
                  "Column 4: start index of primer.")
            print(forward)
        
        #Continue if sequence can be split properly
        else:
            #forward primer
            fwd_final = primerfind(fwd_seq, primer_length, temp_range, gc_range)
            fwd_primers = fwd_final[0]

            #Error handling if no forward sequences found. Still return reverse sequences
            if len(fwd_primers) == 0:
                print("Error: No forward primers found. Extend input sequence upstream for primer search.")
                reverse = rv_comp(primerfind(rv_seq, primer_length, temp_range, gc_range), ended)
                if len(reverse[0]) == 0:
                    print("Error: No reverse primers found. Extend input sequence downstream for primer search.")
                else:
                    rv_primers = reverse[0]
                    rv_gc = reverse[1]
                    rv_temp = reverse[2]
                    rv_index = reverse[3]

                    rv_final = []
                    rv_final.append(rv_primers)
                    rv_final.append(rv_gc)
                    rv_final.append(rv_temp)
                    rv_final.append(rv_index)
                    print("Column 1: Reverse primer sequence. "
                          "Column 2: GC content. "
                          "Column 3: Melting temprature. "
                          "Column 4: start index of primer.")
                    print(rv_final)
            
            #Coninue to find reverse sequences if forward ones are found. 
            else:
                #reverse primer
                rev = primerfind(rv_seq, primer_length, temp_range, gc_range)
                rv = rv_comp(rev, ended)
                rv_primers = rv[0]

                #Error handling if no reverse primer found. Still return forward sequence
                if len(rv_primers) == 0:
                    print("Error: No reverse primers found. Extend input sequence downstream for primer search.")
                    print("Column 1: Forward primer sequence. "
                          "Column 2: GC content. "
                          "Column 3: Melting temprature. "
                          "Column 4: start index of primer.")
                    print(fwd_final)
                
                #Coninue to primer pairing if forward and reverse primers found
                else:
                    #Define variables for primer_pair function
                    fwd_gc = fwd_final[1]
                    fwd_temp = fwd_final[2]
                    fwd_index = fwd_final[3]
                    rv_gc = rv[1]
                    rv_temp = rv[2]
                    rv_index = np.add(rv[3], ended)
                    rv_final = []
                    rv_final.append(rv_primers)
                    rv_final.append(rv_gc)
                    rv_final.append(rv_temp)
                    rv_final.append(rv_index)

                    #primer pair
                    amplified = primer_pair(fwd_final, rv_final, started, ended)

                    #Error Handling if no primer pair found. Still prints forward and reverse primer list
                    if len(amplified[0]) == 0:
                        print("Error: No primer pair found. Extend sequences upstream and downstream for further primer search.")
                        print("Forward and reverse sequences still displayed. Select primers at own discretion.")
                        print("Column 1: Forward primer sequence. "
                              "Column 2: GC content. "
                              "Column 3: Melting temprature. "
                              "Column 4: start index of primer.")
                        print(fwd_final)
                        print("Column 1: Reverse primer sequence. "
                              "Column 2: GC content. "
                              "Column 3: Melting temprature. "
                              "Column 4: start index of primer.")
                        print(rv_final)

                    else:
                        print("Column 1-3: Forward primer sequence, primer start index, melting temperature. " 
                              "Column 4-6: Reverse primer sequence, primer start index, melting temperature. "
                              "Column 7: Length of amplified sequence. Use Gel.Viz to visualize outcome on a gel.")
                        print(amplified)
    
    #Error Handling for invalid primer type input
    else:
        print("Error: Primer type is invalid. Select appropriate primer type.")

