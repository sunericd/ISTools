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
def primerfind(sequence, length, temperature, gcontent):
    
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

def rv_primerfind(sequence, length, temperature, gcontent):
    
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

def primer_pair(f_primer_info, r_primer_info, ended, product_range):
    
    pair = primerpairs([], [], [], [], [], [], [])
    
    for i in range(len(f_primer_info.primer)-1):
        for j in range(len(r_primer_info.primer)-1):
            diff = abs(f_primer_info.Tm[i] - r_primer_info.Tm[j])
            length = r_primer_info.index[j] - f_primer_info.index[i] + ended -2 + len(r_primer_info.primer[j])
            if diff <= 3 and max(product_range) > length and min(product_range) < length:
                pair.add_f_prim(f_primer_info.primer[i])
                pair.add_f_Tm(f_primer_info.Tm[i])
                pair.add_f_ind(f_primer_info.index[i])
                pair.add_r_prim(r_primer_info.primer[j])
                pair.add_r_Tm(r_primer_info.Tm[j])
                pair.add_r_ind(r_primer_info.index[j] + ended -2)
                pair.add_length(length)
                break
            
        if len(pair.f_prim)>=10:
            #display the first 10 primer pair along with lengths and then stop. 
            break
            
    return pair
    
#primer_pair outputs include fwd primer info, rv primer info, and length of sequence

#########For saving outputs
def getFwdFile(primer, gc, tm, index):

    lines = []
    
    lines.extend(('Forward Primers: ' + str(primer) + '%', ' '))
    lines.extend(('GC Content: ' + str(gc) + '%', ' '))
    lines.extend(('Melting Temperature (C): ' + str(tm), ' '))
    lines.extend(('Primer Index: ' + str(index), ' '))

    out_text=''
    for item in lines:
        out_text=out_text+("%s\n" % item)

    return (out_text)

def getRevFile(primer, gc, tm, index):

    lines = []
    
    lines.extend(('Reverse Primers: ' + str(primer) + '%', ' '))
    lines.extend(('GC Content: ' + str(gc) + '%', ' '))
    lines.extend(('Melting Temperature (C): ' + str(tm), ' '))
    lines.extend(('Primer Index: ' + str(index), ' '))

    out_text=''
    for item in lines:
        out_text=out_text+("%s\n" % item)

    return (out_text)

def getBothFile(f_primer, f_gc, f_tm, f_ind, r_primer, r_gc, r_tm, r_ind):
    
    lines = []
    
    lines.extend(('Forward Primers: ' + str(f_primer) + '%', ' '))
    lines.extend(('Forward GC Content: ' + str(f_gc) + '%', ' '))
    lines.extend(('Forward Melting Temperature (C): ' + str(f_tm), ' '))
    lines.extend(('Forward Primer Index: ' + str(f_ind), ' '))
    
    lines.extend(('Reverse Primers: ' + str(r_primer) + '%', ' '))
    lines.extend(('Reverse GC Content: ' + str(r_gc) + '%', ' '))
    lines.extend(('ReverseMelting Temperature (C): ' + str(r_tm), ' '))
    lines.extend(('Reverse Primer Index: ' + str(r_ind), ' '))

    out_text=''
    for item in lines:
        out_text=out_text+("%s\n" % item)

    return (out_text)


def getPairFile(fwd_primer, fwd_tm, fwd_ind, rv_primer, rv_tm, rv_ind, length):
    
    lines = []
    
    lines.extend(('Forward Primers: ' + str(fwd_primer) + '%', ' '))
    lines.extend(('Forward Melting Temperature (C): ' + str(fwd_tm), ' '))
    lines.extend(('Forward Primer Index: ' + str(fwd_ind), ' '))
    
    lines.extend(('Reverse Primers: ' + str(rv_primer) + '%', ' '))
    lines.extend(('Reverse Melting Temperature (C): ' + str(rv_tm), ' '))
    lines.extend(('Reverse Primer Index: ' + str(rv_ind), ' '))

    out_text=''
    for item in lines:
        out_text=out_text+("%s\n" % item)

    return (out_text)
   
##############################################Primer Design Function
def primer_designer(primer_type, sequence_input, indices = [0, 0] , product_range = [300, 400], primer_length = [18, 22], temp_range= [52, 65], gc_range = [40, 60]):
    #primer_types should be 'r','fr', or 'f'. Also the default one should be 'fr'
    #Start and end will be in one list like [start, end]
    #primer length will be like [18, 22], so maybe you can do primer_length=list(range(primer_length[0],primer_length[1])) or something like that?
    # temp range.. why's it 4 to 60?
    #The inputs will be lists of strings, so could you do a input=[int(i) for i in input] to make them all integers? I think it'd be easier to do it here than over there..
    
    temp_range = [float(temp_range[0]), float(temp_range[1])]   
    gc_range = [float(gc_range[0])/100, float(gc_range[1])/100] 
    primer_length = list(range(int(primer_length[0]), int(primer_length[1]) + 1))
    
    if primer_type == 'f':
    
        forward = primerfind(sequence_input, primer_length, temp_range, gc_range)
        if len(forward.primer) == 0:
            raise UserWarning("Error: No forward primers found. Extend input sequence upstream for primer search.")
        else:
            print(forward.all)
            out_text = getFwdFile(forward.primer, forward.gc, forward.Tm, forward.index)
            saveOutput.saveData(out_text, "Primer Designer Analysis ")

    elif primer_type == 'r':
        rev = rv_primerfind(sequence_input, primer_length, temp_range, gc_range)
        if len(rev.primer) == 0:
            raise UerWarning("Error: No reverse primers found. Extend input sequence downstreamstream for primer search.")
        else:
            print(rev.all)
            out_text = getRevFile(rev.primer, rev.gc, rev.Tm, rev.index)
            saveOutput.saveData(out_text, "Primer Designer Analysis ")
            
    else:
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
                    print(rev.all)
                    out_text = getRevFile(rev.primer, rev.gc, rev.Tm, rev.index)
                    saveOutput.saveData(out_text, "Primer Designer Analysis ")
                    raise UserWarning("Error: Region for forward primer search is too short. Check start index or extend sequence upstream of amplification. Reverse primers still found and listed.")
                   
                    
        #Still return forward sequences if reverse sequence unable to be split
        elif len(rv_seq) < min(primer_length):
            raise UserWarning("Error: Region for reverse primer search is too short. Check end index or extend sequence downstream of amplification.")
            forward = primerfind(fwd_seq, primer_length, temp_range, gc_range)
            if len(forward.primer) == 0:
                raise UserWarning("Error: No forward primers found. Extend input sequence upstream for primer search.")
            else:
                print(forward.all)
                out_text = getFwdFile(forward.primer, forward.gc, forward.Tm, forward.index)
                saveOutput.saveData(out_text, "Primer Designer Analysis ")

        #Continue if sequence can be split properly
        else:
            #forward primer
            forward = primerfind(fwd_seq, primer_length, temp_range, gc_range)

            #Error handling if no forward sequences found. Still return reverse sequences
            if len(forward.primer) == 0:
                rev = rv_primerfind(rv_seq, primer_length, temp_range, gc_range)
                raise UserWarning("Error: No forward primers found. Extend input sequence upstream for primer search.")
                if len(rev.primer) == 0:
                    raise UserWarning("Error: No reverse primers found. Extend input sequence downstream for primer search.")
                else:
                    print(rev.all)
                    out_text = getRevFile(rev.primer, rev.gc, rev.Tm, rev.index)
                    saveOutput.saveData(out_text, "Primer Designer Analysis ")
                   
            #Coninue to find reverse sequences if forward ones are found. 
            else:
                #reverse primer
                rev = rv_primerfind(rv_seq, primer_length, temp_range, gc_range)

                #Error handling if no reverse primer found. Still return forward sequence
                if len(rev.primer) == 0:
                    out_text = getFwdFile(forward.primer, forward.gc, forward.Tm, forward.index)
                    saveOutput.saveData(out_text, "Primer Designer Analysis ")
                    raise UserWarning("Error: No Reverse primers found. Extend input sequence downstream for primer search."
                                      "Forward primers still found. Listed below. ")
                    
                #Coninue to primer pairing if forward and reverse primers found
                else:
                    #Find primer pairs
                    amplified = primer_pair(forward, rev, ended, product_range)
                    
                    #Error Handling if no primer pair found. Still prints forward and reverse primer list
                    if len(amplified.f_prim) == 0:
                        print(rev.all)
                        print(forward.all)
                        out_text = getBothFile(forward.primer, forward.gc, forward.Tm, forward.index, rev.primer, rev.gc, rev.Tm, rev.index) 
                        saveOutput.saveData(out_text, "Primer Designer Analysis ")
                        raise UserWarning("Error: No primer pair found. Extend sequence upstream and downstream for further primer search."
                                          "Forward and reverse sequences still displayed. Select primers at own discretion.")               
                    else:
                        print(amplified.all)
                        info ="List 1-3: Forward primer sequence, primer start index, melting temperature. List 4-6: Reverse primer sequence, primer start index, melting temperature. List 7: Length of amplified sequence. Use Gel.Viz to visualize outcome on a gel."
                        out_text = getPairFile(amplified.f_prim, amplified.f_ind, amplified.f_Tm, amplified.r_prim, amplified.r_ind, amplified.r_Tm, amplified.length)
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

#seq1 = 'CTCCATTTTGTATTCTTCTTCACTGATTTAAACAATAAATTCAAATAGCTACAACAAAGAGAAATGAAAAGGAGAACCAGGACAAGAGACTTTAATTTATACACACATTTGCTATTGTCTAAAGTATCAAACATTTGAGGCAAAAAAACTGACAGTAGGTTTTTTTTTTTCATTTCACTTTTTGCCTGAAAGGCATGTTTTATGCATAACAGTAAATGCCACACGCTTTGGAAATATCTCAAAAGTTCTAAAAATTATGCATTTTTTCTTTTTAATTTTGTTTTATTACATATATAAGAGAGGGAAAGAATGTGTGGGCATATATACACACCATATAATATACATGTGGAGGTCAGAGTTCTCCACTTTTGTGGAGTTGGTTCTTTCCTTCTACCTTTAGGGTGTCCCGGGaattaaacttaggttgccaaatgcacatggcaagccttatcctattgagccatcttgctggccTTATTCATGCTTTAAGGAGAAGAAACTAAGGATCAAAGACCTTACACATTTCCCAAAGTTCAAACATCTACAACCTTTACAGTTGCACCCAGACTATCTAACTTTGTCTGTACAAAGCACTTACAATAATGATAATACATCAGTAGTAAAGGTGATAGCCATGGTCATAATGGTAATGATGATAAGGATGGGTGGGTAGTGGTGATGGATGGTTGGTGGTGGTGGTGGTGATGGTAGTGGTAGGGGGGTAATGGTGGTGGTAACGGCTCTGTTGATGCTAAATTGTTCATTGCCCCATTATATTCTAAGTTTCTGAAACTGAAAGATGACTTTTACAGATAAAAGAAGAATTAACACACTTGGGAAATAAAACATGATTCACAGAACAAGAGAAAACATGAACTAATTGTTACTTTAGAGACAACAATGGTCTCTAGAGTGACTATATCCCAGGAGATGATCCACACACACACACACATATATATAAACAGAACCCCCAATTTTTAGTCATATTATATCACCTTTCCAAATGACAATCTCCTTTTCTGGACTAGACTGGGCTCTTTAAGTAAATTTTGTGGTGTATACACATCTGCATACTTCTAGAATTTACTGAATATTGAAGC'
#started = 30
#ended = 411
#seq = seq1[410:len(seq1)]

#seq1='CTGGGGTGTCATTTCCCTACTGAAAGAACCAAGATTCTTGGTGAAATGAATGACTCAGCTCTGAAGAATGAAAAGAATGTGGTAAAGTCACAACATCTTGTCATACCAGATGGTAAGGAAGATATTAGAGACTAGTATATTAAGTTCAGAGACAAATTGAGAAGTCCAGAGACTGTCAATACCTGGGGCAGTAGCTTACTGGTTAAGTGTGCGCACTGCTGTTGCTGAGCACCTGAGTGTGGTTCTCAGCACATAGCTCAACTTATAACTCAAGGTCCAAGGGGGAATTTATTGCTGTCCCCTCCAAAGGTATTTGTGTTCATGTGAACATACCAACACAGACACATGTACAATTTAAAATAATGAAAATAAAATGTAAAACATAAAATTAAATTCTGATAAAAGCATAGTAATTGTATTGGGCTGAAGCACATTAAATATGTTCACACTAGTGAATTTATAAAGGGGGGGATATAAGAAACCAACTATTGTTTTTAAATTATAAAACATATCTTTTTATTATAACTGTAATAGCAAATACAGGAACTATTCTGTTGAATAGATAGCACAAAAATAAGCACCAAAATGAGTGTGATTAAAATTGGAGAGCTCTGAATAAAATTAAAACTAGTGTCAATATTCTGGCCATGATGAGTACAGCTGTTTTCAAAGACATTACTCTCTCCCCGGGAAGGAATTTTAGCAACAGGCTTCCACTTTGTGCTTCAAGTTAAACTTGCAAAACAAAGAAGCCAGGACCTGAATTTGAGAAATAATGTGTGTCATTTGTCTTTCTGTCCCTGAGTTGCCTCACTTGGCAAAACTTTGTTTCAACAAAACTCTTTAAAACAAGCATGTCAAGGCTAGCTCATGATGCATCATGGTGACTAGCTAACAATATGTAATAACCATCATGGCACTGAGAATGATGTTGCTGGTAGTTACTGTGGTGCCTATTGAGATGAAAATG'
#print(primer_design('r', seq))

''' seqr = 'TCATGTCCTAATTTGCCCTTGCTCTTGCA'
seqf = 'ATTTCTGAGGAAAGAGGACTATACCCATTA'
seqa = 'ATTTCTGAGGAAAGAGGACTATACCCATTAGGAAACGAATTGCCCGAGTAGTCTCCTCTGCCGACTTAAACCAACCTTTTTCTATTTCTCTTTTCTTTTCTCCCTCTTTTTTCTCTGTACTAGCATCCAAAAGCAAGCATCCATCCGAGTCCCAGTCGCAATCTCACATCTCCAATTTAACGTATCCATTGCATTTCCTCATTCGGTTTAACTCCTCTGCATTTCTTTTCTGACCCATAGCATTTCTTACATTCCATTGCATCTCCCTTTTACTCTCGTTCAAGACACTGATTTGATACGCTTTCTGTACGATGGCCATATTGAAGGATACCATAATTAGATACGCTAATGCAAGGTATGCTACCGCTAGTGGCACTTCCACCGCCACTGCCGCCTCTGTCAGCGCTGCCTCATGTCCTAATTTGCCCTTGCTCTTGCA'

print(primer_designer('fr', seqa, [30,411], [300, 500], [18,22] ))
print(primer_designer('f', seqf))
print(primer_designer('r', seqr))

 '''