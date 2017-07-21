#07/09/2017 Working on graphics
#07/16/2017 Can save protocol to a specific file path and file name.

import pandas as pd
import os
from graphics import *
from PIL import Image as NewImage
from tkinter import *

def plasmid_builder (plsmds, mrkrs,restriction_sites):

    mrkrs=mrkrs.split()
    
# Parsing biolabs restriction sites
    restriction_enzymes = restriction_sites.values.T.tolist()
    restriction_enzymes = restriction_enzymes[0]
    restriction_sites  = restriction_sites.values.tolist()
    restriction_sites = restriction_sites[1:]

    # Cleaning NaNs
    restriction_seqs = []
    for site in restriction_sites:
        clean_site = [x for x in site if str(x) != 'nan']
        del clean_site[0]
        restriction_seqs.append(clean_site)

    del restriction_enzymes[0]

    # Reading plasmids and markers
    #plsmds = []
    #for seq_record in SeqIO.parse(plasmid_fasta, "fasta"):
    #    plsmd = (repr(seq_record.seq))
    #    plsmds.append(plasmd)

    #mrkrs = []
    #for seq_record in SeqIO.parse(markers_fasta, "fasta"):
    #    mrkr = (repr(seq_record.seq))
    #    mrkrs.append(mrkr)

    plasmid = plsmds

    plasmid_re = []

    for seq in restriction_seqs:
        if seq[0] in plasmid:
            idx = restriction_seqs.index(seq)
            re_name = restriction_enzymes[idx]
            plasmid_re.append(re_name)
    index=1
    for marker in mrkrs:
        marker_re = []
        matched = []
        idx = 0
        for seq in restriction_seqs:
            if seq[0] in marker:
                idx = restriction_seqs.index(seq)
                re_name = restriction_enzymes[idx]
                marker_re.append(re_name)

        for re in marker_re:
            if re in plasmid_re:
                matched.append(re)

        
        protocol=[]
        if len(matched) > 0:
            # Add to PROTOCOL for (1) --- matched[0] gives plasmid name
            protocol.append ("Digest the recipient plasmid and gene with " + matched[0] + " following the manufacturer's instructions")
            protocol.append ("Run on gel electrophoresis to separate the fragments by size; use GelViz to visualize outcome of gel.")
            protocol.append ("Add the separated fragments into the following mixture:")
            protocol.append ("a. " + "Marker of sequence " + marker)
            protocol.append ("b. the recipient plasmid fragment")
            protocol.append ("c. DNA ligase")
            protocol.append ("e. follow the manufacturer's instruction for appropriate temperatures and incubation times")

        elif len(marker) > 400:
            primer1 = marker[0:30] + plasmid[idx:idx+30]
            primer2 = plasmid[idx+30:idx+60] + marker[len(marker)-30:len(marker)]

            protocol.append ("FOLLOWING GIBSON PROTOCOL")
            protocol.append ("Recipient Plasmid Instructions:")
            protocol.append ("a. Obtain" + plasmid_re[0] + " enzyme and recipient plasmid; follow manufacturer's instruction for digestion")
            protocol.append ("b. Conduct Gel electrophoresis to ensure cutting; check fragments with Gel.Viz")
            protocol.append ("Gene instructions:")
            protocol.append ("a. Obtain two primers:")
            protocol.append (str(primer1) + ' and ' + str(primer2))
            protocol.append ("b. Run a Two Step PCR with the appropriate annealing temperatures")
            protocol.append ("c. Conduct Gel electrophoresis to ensure extension; check fragments with Gel.Viz")
            protocol.append ("Final instructions:")
            protocol.append ("Incubate the amplified fragment and recipient plasmid with the following mixture:")
            protocol.append ("a.an exonuclease the chews back 5' ends of the fragment to create overhangs")
            protocol.append ("b. a polymerase to fill in gaps")
            protocol.append ("c. a DNA ligase tht seals the nicks of filled in gaps")
            protocol.append ("d. follow the manufacturer's instructions for appropriate buffers and incubation time")
        
        out_text=''
        for item in protocol:
            out_text=out_text+("%s\n" % item)

        outPg=Toplevel()
        outPg.title("Protocol "+str(index))

        outFrame=Frame(outPg)
        outFrame.pack(fill=X)

        output=Text(outFrame,wrap=WORD)
        output.delete('1.0',END)
        output.insert(END,out_text)
        output.configure(state=DISABLED)
        output.pack(fill=BOTH,expand=True,side=LEFT,padx=5,pady=5)

        yScroll=Scrollbar(outFrame)
        yScroll.pack(side=LEFT,fill=Y)
        yScroll.configure(command=output.yview)

        output.configure(yscrollcommand=yScroll.set)

        pFrame=Frame(outPg)
        pFrame.pack(fill=X,padx=5,pady=5)

        pLabel=Label(pFrame, text="File Path: ")
        pLabel.pack(fill=Y,padx=5,pady=5,side=LEFT)

        pathBox=Text(pFrame,height=1,state=DISABLED)
        pathBox.pack(fill=BOTH,side=LEFT,padx=5,pady=5)

        browseButton=Button(pFrame, text="Browse", command=lambda:getDir(pathBox))
        browseButton.pack(side=RIGHT,fill=Y,padx=5,pady=5)

        fFrame=Frame(outPg)
        fFrame.pack(fill=X,padx=5,pady=5)

        fLabel=Label(fFrame,text="File Name: ")
        fLabel.pack(fill=Y,padx=5,pady=5,side=LEFT)

        fileBox=Text(fFrame,height=1)
        fileBox.pack(fill=BOTH,side=LEFT,padx=5,pady=5)
        
        sButton=Button(fFrame, text="Save", command=lambda:writeFile(out_text,pathBox.get('1.0',END),fileBox.get('1.0',END)))
        sButton.pack(side=RIGHT,fill=Y,padx=5,pady=5)
       
        index=index+1

def writeFile(text,path,fname):
    name=path.strip()+fname.strip()
    try:
        os.chmod('*', 0o755)
        os.remove(name+".txt")
    except OSError:
        pass

    h=open(name.strip()+".txt","w")
    h.write(text)
    h.close()

def getDir(box):
    box.configure(state=NORMAL)
    box.delete('1.0',END)
    pathname=filedialog.askdirectory()
    pathname=pathname+'/'
    box.insert(END,pathname)
    box.configure(state=DISABLED)

#Drawing plasmid
def plasDraw(plasmid, mrkrs):
    total=len(plasmid)
    for i in mrkrs:
        total=total+len(i)
    
    scale=360/total
    padding=50
    
    master=Toplevel()
    master.title("Plasmid")
    w=Canvas(master,width=500,height=500)
    w.pack()

    x=w.winfo_width()
    y=w.winfo_height()

    plasEnd=scale*len(plasmid)

    w.create_arc(padding,padding,x-padding,y-padding,start=0,extent=plasEnd, style=ARC,width=5, outline="black")
    for i in mrkrs:
        w.create_arc(50,50,450,450,start=plasEnd,extent=len(i)*scale,style=ARC,width=5,outline="red")
        #Need to incorporate cuts in where the restriction sites are

#    else:
        # Primer protocol      
   
# Finding
    
# User is presented with a list of restriction enzymes and other reagents and selects what they have in the lab.
# All of the following will be in a for loop for each gene/marker.    
# 1. Restriction enzyme/ligation function:
#    a. Search for all restriction sites on plasmid that were selected by user.
#    b. Search for all restriction sites on genes/markers that were selected by user.
#    c. Search for matches between plasmids and genes.
#    d. If there are matches, print restriction enzyme names and protocol.
#            "Digest the #recipient plasmid and #gene with #restriction enzymes following the manufacturer's instructions"
#            "Run on gel electrophoresis to separate the fragments by size; use GelViz to visualize outcome of gel."
#            "Add the separated fragments into the following mixture:
#                a. the %gene frgment
#                b. the %recipient plasmid fragment
#                c. DNA ligase
#                e. follow the manufacturer's instruction for appropriate temperatures and incubation times 
#
#
# 2. Gibson method:
#   a. Pick first restriction site found on plasmid.
#   b. Copy 30 bps upstream and downstream of restriction site location.
#   c. Copy 20 bps on either end of gene/marker.
#   d. Combine to form 60 bps for each primer.
#   e. Check melting temperature of primers. Discard and repeat if temp < 40C or > 60C.
#   f. Check for primer dimer formation. Discard and repeat if 5 or more nucleotides form dimers.
#   g. Print primer sequences and protocol.
#             "Recipient Plasmid Instructions:
#               a. Obtain %restriction enzyme and %recipientplasmid; follow manufacturer's instruction for digestion
#               b. Conduct Gel electrophoresis to ensure cutting; check with following gel:"
#             "Gene instruction:
#               a. "Obtain primer with the following sequene %primer sequence
#               b. Run a Two Step PCR with the appropriate annealing temperatures
#             "Final instructions:
#              ubate the amplified fragment and recipient plasmid with the following mixture:
#                   a.an exonuclease the chews back 5' ends of the fragment to create overhangs
#                   b. a polymerase to fill in gaps
#                   c. a DNA ligase tht seals the nicks of filled in gaps
#                   d. follow the manufacturer's instructions for appropriate buffers and incubation time

# 3. Restriction Primers
#        Protocol:
#           a. find sequence of restriction enzyme of recipient plasmid
#           b. combine sequence with one end of the gene complemetary sequence
#           c. print primer sequence and protocol:
#             "Obtain primers with sequence: %primer_seq"
#             "
#             "Aliquot an appropriate amount of the %enzyme restriction enzyme into the mixture."
#             "Incubate the mixture for an appropriate duration of time."
#             "Remove the restriction enzymes."
#             "Add the primers into the mixture."