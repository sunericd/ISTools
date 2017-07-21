# 07/09/2017 Modified to take a list of sequences, change window sizes to accommodate more gels, and label columns.

import pandas as pd
from graphics import *
from math import *
from PIL import Image as NewImage


def gel_visualize(plasmid_seqs, re_list,restriction_sites):
    plasmid_seqs = plasmid_seqs.split()
    re_list = re_list.split()

    # ERROR HANDLING
    # 30 well maximum
    if len(plasmid_seqs) > 30:
        raise Exception(
            'You have exceeded the 30 well maximum. Also, check to ensure that the only spaces are between sequences'
            )
    for plasmid_seq in plasmid_seqs:
        # Max ladder is 3KB
        if len(plasmid_seq) > 3000:
            raise Exception(
                'Our gel cannot handle sequences larger than 3KB!'
                )
        # Only capital A,T,C,G
        for nuc in plasmid_seq:
            if nuc not in "ATCG":
                raise Exception(
                    'Sequences must only contain A, T, C, or G!'
                    )

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

    # Matching restriction enzyme names ERROR
    for re in re_list:
        if re not in restriction_enzymes:
            raise Exception(
                str(re) + ' is not in our list. Please check to make sure that the enzyme name matches the names in our Enzyme List.'
                )
    
    ##### Gel Visualization ######
    lengths_list = []
    max_lengths = []
    
    # Finding cut indicies
    for plasmid_seq in plasmid_seqs:
        
        cut_idxs = []
        
        for re in re_list:
            re_index = restriction_enzymes.index(re)
            seq_index = 0
            for cut_site in restriction_seqs[re_index]:
                site_idx = plasmid_seq.find(restriction_seqs[re_index][seq_index])
                seq_index += 1
                cut_idxs.append(site_idx)
    
        cut_idxs.sort()
        start_idx = 0

        cut_idxs.append(len(plasmid_seq))

        lengths = []
        
        # Collecting fragment lengths
        for idx in cut_idxs:
            frag = plasmid_seq[start_idx:idx]
            lengths.append(len(frag))
            start_idx = idx

        # Cleaning lengths
        lengths = [x for x in lengths if str(x) != '0']
        lengths_list.append(lengths)
        max_lengths.append(max(lengths))

    # Selecting appropriate ladder
    #big_ladder = []
    #small_ladder = []
            
    if max(max_lengths) > 1000:
        #print("big")
        bigDraw(lengths_list)

    else:
        #print("Small")
        smallDraw(lengths_list)

#Small Ladder
def smallDraw (lengths):
    numSeq=len(lengths)
    if numSeq>2:
        x=150*0.2*(numSeq+3)
    else:
        x=150
    win=GraphWin('Gel',x,300)

    xWin=win.getWidth()
    yWin=win.getHeight()


    standardX=150*0.3
    standardY=yWin*0.1
    lenBand=150*0.15

    scale=1.35

    listBText=['1.0 kb','0.8 kb', '0.6 kb', '0.5 kb', '0.4 kb', '0.3 kb', '0.2 kb','0.15 kb','0.1 kb','0.075 kb' , '0.05 kb']
    listBand=[float(1000),float(800),float(600),float(500),float(400),float(300),float(200),float(150),float(100),float(75),float(50)]
    listY=[]

    for i in listBand:
        listY.append(float(-1/0.098*(log10(i)-3.8494)))
    largest=max(listY)

    for i in range(len(listY)):
        listY[i]=float(listY[i])/largest*yWin

    for i in range(len(listY)):
        yVal=scale*(listY[i]-min(listY)+15)
        pt=Point(standardX,yVal)
        l=Line(Point(standardX+lenBand,yVal),pt)
        kblabel=Text(Point((standardX)/2,yVal),listBText[i])
        kblabel.setSize(7)

        l.draw(win)
        kblabel.draw(win)
    
    stdlabel=Text(Point(standardX+0.5*lenBand,yWin*0.02),'ladder')
    stdlabel.setSize(7)
    stdlabel.draw(win)
        
    for i in range(numSeq):
        testX=standardX+150*0.2*(i+1)
        testY=yWin*0.1

        numlabel=Text(Point(testX+0.5*lenBand,yWin*0.02),(i+1))
        numlabel.setSize(7)
        numlabel.draw(win)

        listTY=[]

        for j in lengths[i]:
            listTY.append(float(-1/0.098*(log10(j)-3.8494)))

        for j in range(len(listTY)):
            listTY[j]=float(listTY[j])/largest*yWin
        
        for j in range(len(listTY)):
            yVal=scale*(listTY[j]-min(listY)+15)
            pt=Point(testX,yVal)
            l=Line(Point(testX+lenBand,yVal),pt)

            l.draw(win)
        
        # Saving File
#        win.postscript(file="Gel.eps", colormode='color')
#        img = NewImage.open("Gel.eps")

#        PNG_Title = "Gel.png"
#        img.save(PNG_Title, "png")
        
    win.getMouse()
    win.close()

#Large Ladder
def bigDraw (lengths):
    numSeq=len(lengths)
    if numSeq>2:
        x=150*0.2*(numSeq+3)
    else:
        x=150
    win=GraphWin('Gel',x,250)

    xWin=win.getWidth()
    yWin=win.getHeight()

    standardX=150*0.3
    standardY=yWin*0.1
    lenBand=150*0.15

    scale=1.1

    listBText=['3.0 kb','2.0 kb', '1.5 kb', '1.0 kb', '0.8 kb', '0.6 kb', '0.5 kb','0.4 kb','0.3 kb','0.2 kb','0.15 kb' , '0.1 kb']
    listBand=[float(3000),float(2000),float(1500),float(1000),float(800),float(600),float(500),float(400),float(300),float(200),float(150), float(100)]
    listY=[]

    for i in listBand:
        listY.append(float(-1/0.098*(log10(i)-3.8494)))
    largest=max(listY)

    for i in range(len(listY)):
        listY[i]=float(listY[i])/largest*yWin

    for i in range(len(listY)):
        yVal=scale*(listY[i]-min(listY)+15)
        pt=Point(standardX,yVal)
        l=Line(Point(standardX+lenBand,yVal),pt)
        kblabel=Text(Point((standardX)/2,yVal),listBText[i])
        l.draw(win)
        kblabel.draw(win)
        kblabel.setSize(7)
    
    stdlabel=Text(Point(standardX+0.5*lenBand,yWin*0.02),'ladder')
    stdlabel.setSize(7)
    stdlabel.draw(win)

    for i in range(numSeq):
        testX=standardX+150*0.2*(i+1)
        testY=yWin*0.1

        numlabel=Text(Point(testX+0.5*lenBand,yWin*0.02),(i+1))
        numlabel.setSize(7)
        numlabel.draw(win)

        listTY=[]

        for j in lengths[i]:
            listTY.append(float(-1/0.098*(log10(j)-3.8494)))

        for j in range(len(listTY)):
            listTY[j]=float(listTY[j])/largest*yWin
        
        for j in range(len(listTY)):
            yVal=scale*(listTY[j]-min(listY)+15)
            pt=Point(testX,yVal)
            l=Line(Point(testX+lenBand,yVal),pt)

            l.draw(win)
  
    win.getMouse()
    win.close()