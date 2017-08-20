import numpy as np
import os
import mods.saveOutput as saveOutput

################For finding gc content, melting temp, and reverse complement

#Approximate Melting Temperature based on Wallace Rule
#only good for short oligos
def gc_content(primer):
    GC_num = 0
    for nuc in primer:
        if nuc is 'G' or nuc is 'C' or nuc is 'g' or nuc is 'c': 
            GC_num += 1
            GC = round(float(GC_num)/float(len(primer)), 2)
    return GC

def melting_temp(primer):
    GC_num = 0
    for nuc in primer:
        if nuc is 'G' or nuc is 'C' or nuc is 'g' or nuc is 'c': 
            GC_num += 1
            GC = round(float(GC_num)/float(len(primer)), 2)
            Tm= 64.9 +41*(GC_num-16.4)/len(primer)
            Tm=(round(Tm, 2))
    return Tm

def rv_comp(unreversed):
    try:
        complement = ''.join([{'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a', '\n':'', ' ':''}[B] for B in list(unreversed)])
    except:
        assert False, "Make sure only A, T, C, or G!"
    rv_comped = complement[::-1]
    return rv_comped

###############################For storing data
class primer_final(object):
    def __init__(self, primer, gc, Tm, index):
        self.all = primer, gc, Tm, index
        self.primer = primer
        self.gc = gc
        self.Tm = Tm
        self.index = index
    def add_primer(self, new_primer):
        self.primer.append(new_primer)
    def add_gc(self, new_gc):
        self.gc.append(new_gc)
    def add_Tm(self, new_Tm):
        self.Tm.append(new_Tm)
    def add_index(self, new_index):
        self.index.append(new_index)

#################################Primer Find Function 
def primerfind(sequence, length = [20, 19, 21, 18, 22], temperature = [52, 65], gcontent = [.4, .6]):
    
    for nuc in sequence:
        if nuc not in "ATCGatcg":
            raise UserWarning('Sequences must only contain A, a, T, t, C, c, or G, g!')
    
    #error handling for sequence inputs that are too short 
    if len(sequence) < min(length):
        raise UserWarning("Error: Seqeunce input is shorter than primer. Input longer sequence for primer search.")
    
    else:
        #Initialize
        f_primer_info = primer_final([],[],[],[])

        #identify start index
        start = len(sequence) - min(length) + 1

        #Flip through the sequence to find primers
        for i in range(start-1, 2, -1):
            #Separate into appropriate bp length sequences 
            for j in length: #flips through all possible ideal lengths 
                seq = sequence[i:(i + j-1)] #identify sequences
                GC = gc_content(seq)
                Tm = melting_temp(seq)

                #Filter out Tm and Content if not in range
                if Tm >= min(temperature) and Tm <= max(temperature) and GC >= min(gcontent) and GC <= max(gcontent):
                    f_primer_info.add_primer(seq)
                    f_primer_info.add_gc(GC)
                    f_primer_info.add_Tm(Tm)
                    f_primer_info.add_index(i+1)
                    break
        return f_primer_info
########### End of primer find function

def rv_primerfind(sequence, length = [20, 19, 21, 18, 22], temperature = [52, 65], gcontent = [.4, .6]):
    
    rv = primerfind(sequence, length, temperature, gcontent)
    
    unreversed = rv.primer
    primers = []
    
    for primer in unreversed:
        primers.append(rv_comp(primer))
    
    primers.reverse
    gc = rv.gc
    Tm = rv.Tm
    index = rv.index
    
    primers.reverse()
    r_primer_info = primer_final(primers, gc[::-1], Tm[::-1], index[::-1])
    return r_primer_info

###############################For storing data
class primerpairs(object):
    def __init__(self, f_prim, f_Tm, f_ind, r_prim, r_Tm, r_ind, length):
        self.all = f_prim, f_Tm, f_ind, r_prim, r_Tm, r_ind, length
        self.f_prim = f_prim
        self.f_Tm = f_Tm
        self.f_ind = f_ind 
        self.r_prim = r_prim 
        self.r_Tm = r_Tm 
        self.r_ind = r_ind
        self.length = length
    def add_f_prim(self, new_f_prim):
        self.f_prim.append(new_f_prim)
    def add_r_prim(self, new_r_prim):
        self.r_prim.append(new_r_prim)
    def add_f_Tm(self, new_f_Tm):
        self.f_Tm.append(new_f_Tm)
    def add_r_Tm(self, new_r_Tm):
        self.r_Tm.append(new_r_Tm)
    def add_f_ind(self, new_f_ind):
        self.f_ind.append(new_f_ind)
    def add_r_ind(self, new_r_ind):
        self.r_ind.append(new_r_ind)
    def add_length(self, new_length):
        self.length.append(new_length)

################################## Forward and Reverse Sequence Function
def split_seq(split_sequence, started, ended):
    #Error Handling
    if int(started) > len(split_sequence):
        raise UserWarning("Error: Start index is out of range of sequence input. Choose smaller start index.")
    elif int(ended) > len(split_sequence):
        raise UserWarning("Error: End index is out of range of sequence input. Choose smaller end index.")
    elif int(ended) == int(started):
        raise UserWarning("Error: Start and end indices are equal. No sequence is amplified.")
    else:
        fwd_seq = split_sequence[0:started]
        rv_seq = split_sequence[ended-1:len(split_sequence)]
        return fwd_seq, rv_seq
    #function returns fwd_seq and rv_seq

###########################################Forward and Reverse Match Function
#Tm should be +- 2 degrees celsius difference
#optimize for the shortest sequence amplified

def primer_pair(f_primer_info, r_primer_info, ended):
    
    pair = primerpairs([], [], [], [], [], [], [])
    
    for i in range(len(f_primer_info.primer)-1):
        for j in range(len(r_primer_info.primer)-1):
            diff = abs(f_primer_info.Tm[i] - r_primer_info.Tm[j])
            if diff <= 3:
                pair.add_f_prim(f_primer_info.primer[i])
                pair.add_f_Tm(f_primer_info.Tm[i])
                pair.add_f_ind(f_primer_info.index[i])
                pair.add_r_prim(r_primer_info.primer[j])
                pair.add_r_Tm(r_primer_info.Tm[j])
                pair.add_r_ind(r_primer_info.index[j] + ended -2)
                pair.add_length(r_primer_info.index[j] - f_primer_info.index[i] + ended -2 + len(r_primer_info.primer[j]))
                break
            
        if len(pair.f_prim)>=10:
            #display the first 10 primer pair along with lengths and then stop. 
            break
            
    return pair
    
#primer_pair outputs include fwd primer info, rv primer info, and length of sequence

def user_info(what):
    if what == 'f':
        return ("List 1: Forward primer sequence. "
               "List 2: GC content. "
               "List 3: Melting temprature. "
               "List 4: start index of primer.")
    elif what == 'r':
        return ("List 1: Reverse primer sequence. "
               "List 2: GC content. "
               "List 3: Melting temprature. "
               "List 4: start index of primer.")

##############################################Primer Design Function
def primer_designer(primer_type, sequence_input, indices = [0, 0] , primer_length = [18, 22], Tm_range= [52, 65], GC_range = [.4, .6]):
    #primer_types should be 'r','fr', or 'f'. Also the default one should be 'fr'
    print(sequence_input, indices, primer_length, Tm_range,GC_range)
    temp_range = [float(Tm_range[0]), float(Tm_range[1])]   
    gc_range = [float(GC_range[0]), float(GC_range[1])] 
    primer_length = list(range(int(primer_length[0]), int(primer_length[1]) + 1))
    
    if primer_type == 'f':
        print('f')
        forward = primerfind(sequence_input, primer_length, temp_range, gc_range)
        if len(forward.primer) == 0:
            raise UserWarning("Error: No forward primers found. Extend input sequence upstream for primer search.")
        else:
            out_text = [user_info('f'), str(forward.primer), str(forward.gc), str(forward.Tm), str(forward.index)]
            saveOutput.saveData(out_text, "Primer Designer Analysis ")

    elif primer_type == 'r':
        print('r')
        rev = rv_primerfind(sequence_input, primer_length, temp_range, gc_range)
        if len(rev.primer) == 0:
            raise UserWarning("Error: No forward primers found. Extend input sequence upstream for primer search.")
        else:
            out_text = [user_info('r'), str(rev.primer), str(rev.gc), str(rev.Tm), str(rev.index)]
            saveOutput.saveData(out_text, "Primer Designer Analysis ")
    else:
        print('fr')
        started = min(int(indices[0]), int(indices[1]))
        ended = max(int(indices[0]), int(indices[1]))
        splitted = split_seq(sequence_input, started, ended)
        fwd_seq = splitted[0]
        rv_seq = splitted[1]
        
        #Error handling if sequence unable to be split
        #Still return reverse sequences if forward sequence unable to be split
        if len(fwd_seq) < min(primer_length):
            if len(rv_seq) > min(primer_length):
                rev = rv_primerfind(rv_seq, primer_length, temp_range, gc_range)
                if len(rev.primer) == 0:
                    raise UserWarning("Error: Region for forward primer search is too short. Check start index or extend sequence upstream of amplification. Error: No reverse primers found. Extend input sequence upstream for primer search.")
                else:
                    out_text = [user_info('r'), str(rev.primer), str(rev.gc), str(rev.Tm), str(rev.index)]
                    saveOutput.saveData(out_text, "Primer Designer Analysis ")
                    raise UserWarning("Error: Region for forward primer search is too short. Check start index or extend sequence upstream of amplification. Reverse primers still found and listed.")
                   
                    
        #Still return forward sequences if reverse sequence unable to be split
        elif len(rv_seq) < min(primer_length):
            raise UserWarning("Error: Region for reverse primer search is too short. Check end index or extend sequence downstream of amplification.")
            forward = primerfind(fwd_seq, primer_length, temp_range, gc_range)
            if len(forward.primer) == 0:
                raise UserWarning("Error: No forward primers found. Extend input sequence upstream for primer search.")
            else:
                out_text = [user_info('f'), str(forward.primer), str(forward.gc), str(forward.Tm), str(forward.index)]
                saveOutput.saveData(out_text, "Primer Designer Analysis ")
        
        #Continue if sequence can be split properly
        else:
            #forward primer
            fwd = primerfind(fwd_seq, primer_length, temp_range, gc_range)

            #Error handling if no forward sequences found. Still return reverse sequences
            if len(fwd.primer) == 0:
                rv = rv_primerfind(rv_seq, primer_length, temp_range, gc_range)
                raise UserWarning("Error: No forward primers found. Extend input sequence upstream for primer search.")
                if len(rv.primer) == 0:
                    raise UserWarning("Error: No reverse primers found. Extend input sequence downstream for primer search.")
                else:
                    out_text = [user_info('r'), str(rv.primer), str(rv.gc), str(rv.Tm), str(rv.index)]
                    saveOutput.saveData(out_text, "Primer Designer Analysis ")
            
            #Coninue to find reverse sequences if forward ones are found. 
            else:
                #reverse primer
                rv = rv_primerfind(rv_seq, primer_length, temp_range, gc_range)

                #Error handling if no reverse primer found. Still return forward sequence
                if len(rv.primer) == 0:
                    out_text = [user_info('f'), str(fwd.primer), str(fwd.gc), str(fwd.Tm), str(fwd.index)]
                    saveOutput.saveData(out_text, "Primer Designer Analysis ")
                    raise UserWarning("Error: No Reverse primers found. Extend input sequence downstream for primer search."
                                      "Forward primers still found. Listed below. ")
                
                #Coninue to primer pairing if forward and reverse primers found
                else:
                    #Find primer pairs
                    amplified = primer_pair(fwd, rv, ended)
                    
                    #Error Handling if no primer pair found. Still prints forward and reverse primer list
                    if len(amplified.f_prim) == 0:
                        raise UserWarning("Error: No primer pair found. Extend sequence upstream and downstream for further primer search."
                                          "Forward and reverse sequences still displayed. Select primers at own discretion.")
                        out_text = [user_info('f'), str(fwd.primer), str(fwd.gc), str(fwd.Tm), str(fwd.index), 
                                    user_info('r'), str(rv.primer), str(rv.gc), str(rv.Tm), str(rv.index)]
                        saveOutput.saveData(out_text, "Primer Designer Analysis ")

                    else:
                        info ="List 1-3: Forward primer sequence, primer start index, melting temperature. List 4-6: Reverse primer sequence, primer start index, melting temperature. List 7: Length of amplified sequence. Use Gel.Viz to visualize outcome on a gel."
                        out_text = [info, str(amplified.f_prim), str(amplified.f_ind), str(amplified.f_Tm), str(amplified.r_prim), 
                                    str(amplified.r_ind), str(amplified.r_Tm), str(amplified.length)]
                        saveOutput.saveData(out_text, "Primer Designer Analysis ")
    
    #Error Handling for invalid primer type input
    
######For testing
#seqf = 'TGGATTTCTGAGGAAAGAGGACTATACCCA'
#print(rv_comp(seqf))
#m = primerfind(seqf)

#new = []
#n = rv_comp(seqf)
#new.append(seqf)
#new.append(n)
#print(new)
#g = primerfind(seqf)
#print(g.all)

#v = rv_primerfind(seqf)
#print(v.all)

#seq1 = 'ATTTCTGAGGAAAGAGGACTATACCCATTAGGAAACGAATTGCCCGAGTAGTCTCCTCTGCCGACTTAAACCAACCTTTTTCTATTTCTCTTTTCTTTTCTCCCTCTTTTTTCTCTGTACTAGCATCCAAAAGCAAGCATCCATCCGAGTCCCAGTCGCAATCTCACATCTCCAATTTAACGTATCCATTGCATTTCCTCATTCGGTTTAACTCCTCTGCATTTCTTTTCTGACCCATAGCATTTCTTACATTCCATTGCATCTCCCTTTTACTCTCGTTCAAGACACTGATTTGATACGCTTTCTGTACGATGGCCATATTGAAGGATACCATAATTAGATACGCTAATGCAAGGTATGCTACCGCTAGTGGCACTTCCACCGCCACTGCCGCCTCTGTCAGCGCTGCCTCATGTCCTAATTTGCCCTTGCTCTTGCA'
#started = 30
#ended = 411
#seq = seq1[410:len(seq1)]

#print(primer_design('r', seq))

#seqr = 'TCATGTCCTAATTTGCCCTTGCTCTTGCA'
#seqf = 'ATTTCTGAGGAAAGAGGACTATACCCATTA'

#print(primer_design('f', seqf))
#print(primer_design('r', seqr))
#print(primer_design('fr', seq1))