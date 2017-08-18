import numpy as np

#Christina: Please change ALL "Exception"s or "print"s that you want to show to the user to "UserWarning"s!!!!
#Christina: What is this function's last product? Please save this as a string, like you see in PlasmidBUILDR. This is what the user will see. then import mods.saveOutput into this file and use it to set up the saving process for users.
#Don't make the overall function return stuff -- it won't go anywhere, we won't use it

################For finding gc content, melting temp, and reverse complement

#Approximate Melting Temperature based on Wallace Rule
#only good for short oligos
def gc_content(primer):
    GC_num = 0
    for nuc in primer:
        if nuc not in "ATCGatcg":
            raise Exception('Sequences must only contain A, T, C, or G!')
        if nuc is 'G' or nuc is 'C' or nuc is 'g' or nuc is 'c': 
            GC_num += 1
            GC = round(float(GC_num)/float(len(primer)), 2)
    return GC

def melting_temp(primer):
    GC_num = 0
    for nuc in primer: #Christina: This is the second time so far this exception occurs, so maybe do the check in the primer_design function before running any of these smaller functions on the code?
        if nuc not in "ATCGatcg":
            raise Exception('Sequences must only contain A, T, C, or G!')
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
        print("Error: Start index is out of range of sequence input. Choose smaller start index.")
    elif int(ended) > len(split_sequence):
        print("Error: End index is out of range of sequence input. Choose smaller end index.")
    elif int(ended) == int(started):
        print("Error: Start and end indices are equal. No sequence is amplified.")
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
        print("Column 1: Forward primer sequence. "
              "Column 2: GC content. "
              "Column 3: Melting temprature. "
              "Column 4: start index of primer.")
    elif what == 'r':
        print("Column 1: Reverse primer sequence. "
              "Column 2: GC content. "
              "Column 3: Melting temprature. "
              "Column 4: start index of primer.")

##############################################Primer Design Function
def primer_design(primer_type, sequence_input, indices = [0, 0] , primer_length = [18, 22], Tm_range= [52, 58], GC_range = [40, 60]):
    
    temp_range = [float(Tm_range[0]), float(Tm_range[1])]   
    gc_range = [float(GC_range[0])/float(100), float(GC_range[1])/float(100)] 
    primer_length = list(range(int(primer_length[0]), int(primer_length[1]) + 1))
    
    if primer_type == 'f':
    
        forward = primerfind(sequence_input, primer_length, temp_range, gc_range)
        if len(forward.primer) == 0:
            print("Error: No forward primers found. Extend input sequence upstream for primer search.")
        else:
            user_info('f') 
            return forward.all

    elif primer_type == 'r':
        rev = rv_primerfind(sequence_input, primer_length, temp_range, gc_range)
        if len(rev.primer) == 0:
            print("Error: No forward primers found. Extend input sequence upstream for primer search.")
        else:
            user_info('r')
            return rev.all
    
    else:
        started = min(int(indices[0]), int(indices[1]))
        ended = max(int(indices[0]), int(indices[1]))
        splitted = split_seq(sequence_input, started, ended)
        fwd_seq = splitted[0]
        rv_seq = splitted[1]
        print(fwd_seq)
        print(rv_seq)
        
        #Error handling if sequence unable to be split
        #Still return reverse sequences if forward sequence unable to be split
        if len(fwd_seq) < 18:
            print("Error: Cannot identify forward primers. Extend sequence upstream of amplification and check start index.")
            if len(rv_seq) > 18:
                rev = rv_primerfind(rv_seq, primer_length, temp_range, gc_range)
                if len(rev.primer) == 0:
                    print("Error: No forward primers found. Extend input sequence upstream for primer search.")
                else:
                    user_info('r')
                    return rev.all
                    
        #Still return forward sequences if reverse sequence unable to be split
        elif len(rv_seq) < 18:
            print("Error: Cannot identify reverse primers. Extend sequence downstream of amplification and check end index.")
            forward = primerfind(fwd_seq, primer_length, temp_range, gc_range)
            if len(forward.primer) == 0:
                print("Error: No forward primers found. Extend input sequence upstream for primer search.")
            else:
                user_info('f')
                return forward.all
        
        #Continue if sequence can be split properly
        else:
            #forward primer
            fwd = primerfind(fwd_seq, primer_length, temp_range, gc_range)

            #Error handling if no forward sequences found. Still return reverse sequences
            if len(fwd.primer) == 0:
                print("Error: No forward primers found. Extend input sequence upstream for primer search.")
                rv = rv_primerfind(rv_seq, primer_length, temp_range, gc_range)
                if len(rv.primer) == 0:
                    print("Error: No reverse primers found. Extend input sequence downstream for primer search.")
                else:
                    user_info('r')
                    return rv.all
            
            #Coninue to find reverse sequences if forward ones are found. 
            else:
                #reverse primer
                rv = rv_primerfind(rv_seq, primer_length, temp_range, gc_range)

                #Error handling if no reverse primer found. Still return forward sequence
                if len(rv.primer) == 0:
                    print("Error: No Reverse primers found. Extend input sequence downstream for primer search.")
                    print("Forward primers still found. Listed below. ")
                    user_info('f')
                    return fwd.all
                
                #Coninue to primer pairing if forward and reverse primers found
                else:
                    #Find primer pairs
                    amplified = primer_pair(fwd, rv, ended)
                    
                    #Error Handling if no primer pair found. Still prints forward and reverse primer list
                    if len(amplified.f_prim) == 0:
                        print("Error: No primer pair found. Extend sequence upstream and downstream for further primer search.")
                        print("Forward and reverse sequences still displayed. Select primers at own discretion.")
                        user_info('f')
                        print(fwd.all)
                        user_info('r')
                        return rv.all

                    else:
                        print("Column 1-3: Forward primer sequence, primer start index, melting temperature. " 
                              "Column 4-6: Reverse primer sequence, primer start index, melting temperature. "
                              "Column 7: Length of amplified sequence. Use Gel.Viz to visualize outcome on a gel.")
                        return amplified.all
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

seq1 = 'ATTTCTGAGGAAAGAGGACTATACCCATTAGGAAACGAATTGCCCGAGTAGTCTCCTCTGCCGACTTAAACCAACCTTTTTCTATTTCTCTTTTCTTTTCTCCCTCTTTTTTCTCTGTACTAGCATCCAAAAGCAAGCATCCATCCGAGTCCCAGTCGCAATCTCACATCTCCAATTTAACGTATCCATTGCATTTCCTCATTCGGTTTAACTCCTCTGCATTTCTTTTCTGACCCATAGCATTTCTTACATTCCATTGCATCTCCCTTTTACTCTCGTTCAAGACACTGATTTGATACGCTTTCTGTACGATGGCCATATTGAAGGATACCATAATTAGATACGCTAATGCAAGGTATGCTACCGCTAGTGGCACTTCCACCGCCACTGCCGCCTCTGTCAGCGCTGCCTCATGTCCTAATTTGCCCTTGCTCTTGCA'
#started = 30
#ended = 411
#seq = seq1[410:len(seq1)]

#print(primer_design('r', seq))

seqr = 'TCATGTCCTAATTTGCCCTTGCTCTTGCA'
seqf = 'ATTTCTGAGGAAAGAGGACTATACCCATTA'

print(primer_design('f', seqf))
print(primer_design('r', seqr))
print(primer_design('fr', seq1, [30, 411]))


#def prim_match(f_primer_info, r_primer_info):
 
    #if len(f_primer_info.primer) < len(r_primer_info.primer):
    #    first = f_primer_info
    #    second = r_primer_info
    #elif len(f_primer_info.primer) == len(r_primer_info.primer):
    #    first = f_primer_info
    #    second = r_primer_info
    #elif len(f_primer_info.primer) > len(r_primer_info.primer):
    #    second = f_primer_info
    #    first = r_primer_info
   # for i in range(len(f_primer_info.Tm)):
   #     
    #    r_prim = []
    #    r_Tm = []
    #    r_ind = []
    #    length = []
    #    for j in range(len(r_primer_info.Tm)-1): 
    #primer can only be paired with the reverse primer that is at most 15bp + amplified length. 
    #if nothing fits, it moves to the next primer and pairs sequences
        
     #       diff = abs(f_primer_info.Tm[i] - r_primer_info.Tm[j])
     #       if diff<3:
     #           pair.add_f_prim(f_primer_info.primer[i])
     #           pair.add_f_Tm(f_primer_info.Tm[i])
     #           pair.add_f_ind(f_primer_info.index[i])
     #           while len(r_prim) <= 3:
     #               r_prim.append(r_primer_info.primer[j])
     #               r_Tm.append(r_primer_info.Tm[j])
     #               r_ind.append(r_primer_info.index[j] + ended -2)
     #               length.append(r_primer_info.index[j] - f_primer_info.index[i] + ended -2 + len(r_primer_info.primer[j]))
     #           break 
                   
                    #if tm of fwd 20bp primer and first reverse primer checked is good, then it moves on the the second primer
                    #this is so that not all sequences have similar fwd primers with only a few bp difference                    

      #  pair.add_r_prim(r_prim)
      #  pair.add_r_Tm(r_Tm)
      #  pair.add_r_ind(r_ind)
      #  pair.add_length(length)

       # if len(pair.r_ind)>=10:
       #     #display the first 10 primer pair along with lengths and then stop. 
       #     break