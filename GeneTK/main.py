#8/15/17 Added GUI for primer design. Changed formatting from pack to grid.
#07/13/2017 Tab view for tools. Additional changes to reflect structure from Exe Toolkit.
#07/08/2017 allows multiple files/sequences to be uploaded.

from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
from PIL import Image, ImageTk
from Bio import SeqIO
import mods.gel_visualizer as gv
import os
import mods.plasmid_builder as pb
import pandas as pd
import mods.seqprop as sp
import mods.MSM as msm

class Example(Frame):   
  
    def __init__(self, parent):
        Frame.__init__(self, parent)   
         
        self.parent = parent
        self.initUI()
                
    def initUI(self):
        self.parent.title('Integrated Sciences Gene Toolkit')
        self.style=Style()
        self.style.theme_use('default')  

        self.pack(fill=BOTH, expand=1)
        
        mlabel = Label(self, text='Integrated Sciences Gene Toolkit\n\nDo NOT close this window.', font=('Helvetica', 16),justify=CENTER,wraplength=350)
        mlabel.grid(row=0,column=0,pady=15,padx=10,sticky=S+N+E+W)

        #Navigation
        sButton=Button(self,text='Start',command=self.nav)
        sButton.grid(row=1,column=0,padx=15,pady=5)

        qButton=Button(self,text='Quit',command=self.quit)
        qButton.grid(row=2,column=0,padx=15,pady=5)

    # Needed for PyInstaller to read files: filepath = resource_path()

    def resource_path(self, relative_path):
        ''' Get absolute path to resource, works for dev and for PyInstaller '''
        try:
            # PyInstaller creates a temp folder and stores path in _MEIPASS
            base_path = sys._MEIPASS
        except Exception:
            base_path = os.path.abspath('.')

        return os.path.join(base_path, relative_path)

    def quit(self):
        quit()

    def nav(self):
        top=Toplevel()
        top.title('Integrated Sciences Gene Toolkit')

        reData= pd.read_csv(self.resource_path('data/restriction_sites2.csv'), sep=',') 

        #######################ALL Page Declarations###################

        ntbk=Notebook(top)
        mainPg=Frame(ntbk)
        gelPg=Frame(ntbk)
        plasPg=Frame(ntbk)
        credPg=Frame(ntbk)
        msmPg=Frame(ntbk)
        seqPg=Frame(ntbk)
        primPg=Frame(ntbk)

        #########################Front page#################################3
        
        #Insert personalized logo
        pic = PhotoImage(file=self.resource_path('data/logo.gif'))
        pic.zoom(10,10)
        mainPg_l1=Label(mainPg,image=pic)
        mainPg_l1.image=pic
        mainPg_l1.configure(image=pic)
        mainPg_l1.grid(row=0, column=0,padx=5, sticky=E+W+S+N)

        mainPg_l2=Label(mainPg,text='Welcome to the Integrated Sciences Gene Toolkit. Click on the tabs above to alternate between tools. \n\nSend any questions to support@integratedsciences.org\n\nVersion 1.0.0\n\nLast updated Aug 15, 2017.', font=('Helvetica', 12), wraplength=250)
        mainPg_l2.grid(row=0, column=1, padx=5, sticky=E+W+S+N)
        
        #!!!Add a label, maybe beginning instructions.

        ##############################Gel.Viz Page#################################
        #Gene sequence input
        gelPg_l1 = Label(gelPg, text='Sequence(s) (DNA bases only). Separate different sequences with spaces.')
        gelPg_l1.grid(row=0,column=0,padx=5,pady=5,sticky=W)
        
        gelPg_geneBox = Text(gelPg, height=10)     
        gelPg_geneBox.grid(row=1,column=0, sticky=E+W+S+N, columnspan=2,padx=5,pady=5)   
        
        #Restriction Enzyme list
        gelPg_l2 = Label(gelPg, text='Restriction Enzyme(s)')
        gelPg_l2.grid(row=2,column=0,padx=5,pady=5,sticky=W)

        allREs=StringVar()
        allREs.set(self.getREs(reData))

        gelPg_reList = Listbox(gelPg,selectmode='multiple',listvariable=allREs) 
        gelPg_reList.grid(row=3,column=0,sticky=W+S+N+E,columnspan=2,padx=5,pady=5)

        yScrollRE=Scrollbar(gelPg)
        yScrollRE.grid(row=3,column=1,sticky=E+N+S,pady=5)
        yScrollRE.configure(command=gelPg_reList.yview)

        gelPg_reList.configure(yscrollcommand=yScrollRE.set)
        
        #Navigation
        gelPg_clr=Button(gelPg,text='Clear',command=lambda:[f() for f in [gelPg_reList.selection_clear(0, END),gelPg_geneBox.delete('1.0', END)]])
        gelPg_clr.grid(row=4,column=0,sticky=W,padx=5,pady=5)
        
        gelPg_ent= Button(gelPg,text='Enter', command=lambda:self.runGV(gelPg_geneBox,gelPg_reList.curselection(),reData))
        gelPg_ent.grid(row=4,column=1,sticky=E,padx=5,pady=5)

        #Upload Gene Sequence File Button.
        gelPg_upload=Button(gelPg,text='Upload Gene Sequence(s)', command=lambda: self.seqUpload(gelPg_geneBox,True))
        gelPg_upload.grid(row=0,column=1,sticky=E,padx=5,pady=5)

        ####################################Plasmid BUILDR##################################

        #Plasmid sequence input        
        plasPg_l1 = Label(plasPg, text='One Recipient Plasmid Sequence (DNA bases only)')
        plasPg_l1.grid(row=0,column=0, padx=5, pady=5,sticky=W)
        
        plasPg_geneBox = Text(plasPg, height=10)
        plasPg_geneBox.grid(row=1,column=0,columnspan=2,padx=5,pady=5,sticky=E+W+S+N)
        
        #Restriction enzymes input        
        plasPg_l2 = Label(plasPg, text='Markers or Genes (Separate sequences with spaces)')
        plasPg_l2.grid(row=2,column=0, padx=5, pady=5,sticky=W)

        plasPg_mrkr = Text(plasPg, height=10)
        plasPg_mrkr.grid(row=3,column=0, columnspan=2, padx=5, pady=5,sticky=W+E+S+N)

        #Navigation
        plasPg_clr=Button(plasPg,text='Clear',command=lambda:[f() for f in [plasPg_mrkr.delete('1.0', END),plasPg_geneBox.delete('1.0', END)]])
        plasPg_clr.grid(row=4,column=0, padx=5, pady=5,sticky=W)
        
        plasPg_ent= Button(plasPg,text='Enter', command=lambda:self.runBUILDR(plasPg_geneBox,plasPg_mrkr,reData))
        plasPg_ent.grid(row=4,column=1, padx=5, pady=5,sticky=E)

        plasPg_upl1=Button(plasPg,text='Upload Recipient Plasmid Sequence', command=lambda: self.seqUpload(plasPg_geneBox,False))
        plasPg_upl1.grid(row=0,column=1, padx=5, pady=5,sticky=E)

        plasPg_upl2=Button(plasPg,text='Upload Marker/Gene Sequence(s)', command=lambda: self.seqUpload(plasPg_mrkr,True))
        plasPg_upl2.grid(row=2,column=1, padx=5, pady=5,sticky=E)

        ##############################MSM####################################

        #Inputs
        msmPg_l1 = Label(msmPg, text='Select a .csv or .tsv file containing the names and the sequences.')
        msmPg_l1.grid(row=0,column=0, padx=5, pady=5,sticky=W,columnspan=2)

        msmPg_l2 = Label(msmPg, text='Filename: ')
        msmPg_l2.grid(row=1,column=0, padx=5, pady=5,sticky=W)

        msmPg_browse=Button(msmPg, text='Browse', command=lambda: self.getDir(msmPg_box))
        msmPg_browse.grid(row=1,column=1, padx=5, pady=5,sticky=E)

        msmPg_box = Text(msmPg, height=1, state=DISABLED)
        msmPg_box.grid(row=2,column=0, columnspan=2, padx=5, pady=5,sticky=N+S+W)

        msmPg_l3 = Label(msmPg, text='Motif Length (Default is 10): ')
        msmPg_l3.grid(row=3,column=0, padx=5, pady=5,sticky=W)

        msmPg_motif = Text(msmPg, height=1)
        msmPg_motif.insert(END, '10')
        msmPg_motif.grid(row=4,column=0,padx=5, pady=5,sticky=W+N+S+E,columnspan=2)

        msmPg_l4 = Label(msmPg, text='Match Threshold (Default is 3): ')
        msmPg_l4.grid(row=5,column=0, padx=5, pady=5,sticky=W)

        msmPg_thresh = Text(msmPg, height=1)
        msmPg_thresh.insert(END, '3')        
        msmPg_thresh.grid(row=6,column=0, padx=5, pady=5,sticky=W+N+S+E,columnspan=2)

        #Navigation
        msmPg_clr=Button(msmPg,text='Clear',command=lambda:[f() for f in [self.delBox(msmPg_box),msmPg_motif.delete('1.0',END),msmPg_thresh.delete('1.0',END)]])
        msmPg_clr.grid(row=7,column=0, padx=5, pady=5,sticky=W)
        
        msmPg_ent= Button(msmPg,text='Enter', command=lambda:self.runMSM(msmPg_box,msmPg_motif,msmPg_thresh))
        msmPg_ent.grid(row=7,column=1, padx=5, pady=5,sticky=E)

        #############################SeqProp##################################

        seqPg_l1 = Label(seqPg, text='Sequence(s) (DNA bases only). Separate different sequences with spaces.')
        seqPg_l1.grid(row=0,column=0, padx=5, pady=5,sticky=W)
        
        seqPg_box = Text(seqPg, height=10)
        seqPg_box.grid(row=1,column=0, columnspan=2, padx=5, pady=5,sticky=E+W+N+S)

        #Navigation        
        seqPg_ul1 = Button(seqPg, text= 'Upload Sequence(s)', command = lambda: self.seqUpload(seqPg_box, True))
        seqPg_ul1.grid(row=0,column=1, padx=5, pady=5,sticky=E)    
        
        seqPg_clr=Button(seqPg,text='Clear',command=lambda:[f() for f in [seqPg_box.delete('1.0', END)]])
        seqPg_clr.grid(row=2,column=0, padx=5, pady=5,sticky=W)
        
        seqPg_ent= Button(seqPg,text='Enter', command=lambda:self.runSeqProp(seqPg_box,reData))
        seqPg_ent.grid(row=2,column=1, padx=5, pady=5,sticky=E)

        #############################PrimerDesign##################################
        primPg_l1 = Label(primPg, text='Sequence to amplify AND surrounding regions (DNA bases only)')
        primPg_l1.grid(row=0,column=0, padx=5, pady=5,sticky=W)
        
        primPg_box = Text(primPg, height=10)
        primPg_box.grid(row=1,column=0, columnspan=2, padx=5, pady=5,sticky=E+N+S+W)

        primPg_l2 = Label(primPg, text='Start and end indices of sequence to amplify (Ex. 196-291)')
        primPg_l2.grid(row=2,column=0, padx=5, pady=5,sticky=W)

        primPg_idx = Text(primPg, height=1)      
        primPg_idx.grid(row=3,column=0, columnspan=2, padx=5, pady=5,sticky=E+W+S+N)

        primPg_l3 = Label(primPg, text='GC Content (Default is 40-60%)')
        primPg_l3.grid(row=4,column=0, padx=5, pady=5,sticky=W)

        primPg_gc = Text(primPg, height=1)    
        primPg_gc.insert('1.0', '40-60')
        primPg_gc.grid(row=5,column=0, columnspan=2, padx=5, pady=5,sticky=E+N+S+W)

        primPg_l4 = Label(primPg, text='Melting Temperature (Default is 42-46C)')
        primPg_l4.grid(row=6,column=0, padx=5, pady=5,sticky=W)

        primPg_tm = Text(primPg, height=1)    
        primPg_tm.insert('1.0', '42-46')
        primPg_tm.grid(row=7,column=0, columnspan=2, padx=5, pady=5,sticky=W+S+E+N)

        primPg_l5 = Label(primPg, text='Primer length (Default is 20)')
        primPg_l5.grid(row=8,column=0, padx=5, pady=5,sticky=W)

        primPg_len = Text(primPg, height=1)    
        primPg_len.insert('1.0', '20')
        primPg_len.grid(row=9,column=0, columnspan=2, padx=5, pady=5,sticky=E+W+S+N)

        #Navigation        
        primPg_ul1 = Button(primPg, text= 'Upload Sequence', command = lambda: self.seqUpload(primPg_box, True))     
        primPg_ul1.grid(row=0,column=1, padx=5, pady=5,sticky=E)    
        
        primPg_clr=Button(primPg,text='Clear',command=lambda:[f() for f in [primPg_box.delete('1.0', END), primPg_idx.delete('1.0',END),primPg_gc.delete('1.0',END),primPg_tm.delete('1.0',END),primPg_len.delete('1.0',END)]])
        primPg_clr.grid(row=10,column=0, padx=5, pady=5,sticky=W)
        
        primPg_ent= Button(primPg,text='Enter', command=lambda:self.runPrim(primPg_box, primPg_idx, primPg_gc, primPg_tm, primPg_len))
        primPg_ent.grid(row=10,column=1, padx=5, pady=5,sticky=E)

        ##########################Credits#################################        
        cred = Text(credPg,wrap=WORD)
        cred.insert(END, 'The Integrated Sciences Group is dedicated to the development of free, open-source tools for a range of scientific research. We ship large toolkits with user inteferfaces in order to make the application of our tools as seamless as possible. We are currently working to improve our current tools and to build new ones!\n\ngel.Viz is a gel visualization tool that generates an image of a gel given a sequence and restriction enzymes. This tool is the perfect resource to plan and verify experiments involving gel electrophoresis.\n\nThe Plasmid BUILDR, or Benchmark Utility for Integrating Long DNA Reads, generates a protocol to build a desired plasmid given the sequences of the recepient plasmid and any insert genes or markers.\n\nDevelopers: Yein Christina Park, Eric Sun, and Yi Chen.\n\nVisit our website at integratedsciences.org.\n\nQuestions? Concerns? Email us at support@integratedsciences.org.\n\nMany thanks to Jimmy Thai and Siavash Zamirpour.\n\n\nCopyright 2017 Integrated Sciences\n\nPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \'Software\'), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\nThe above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \'AS IS\', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.')
        cred.grid(row=0,column=0, pady=5,sticky=E+W+S+N)
        
        yScrollCredits=Scrollbar(credPg)
        yScrollCredits.grid(row=0,column=1, pady=5,sticky=E+N+S)
        yScrollCredits.config(command=cred.yview)
        cred.config(state=DISABLED,yscrollcommand=yScrollCredits.set)

        #############################Compile notebook############

        ntbk.add(mainPg, compound=LEFT, text='Welcome',padding=5)
        ntbk.add(gelPg,text='gel.Viz', padding=5)
        ntbk.add(msmPg,text='Multiple Sequence Mapper',padding=5)
        ntbk.add(plasPg,text='Plasmid BUILDR', padding=5)
        ntbk.add(primPg,text='Primo',padding=5)
        ntbk.add(seqPg,text='SeqProp',padding=5)
        ntbk.add(credPg,text='Credits', padding=5)
        
        ntbk.pack(fill=BOTH)
        self.pack(fill=BOTH,expand=True)        

    def getREs(self,reData):
        re = reData.values.T.tolist()
        re = re[0]
        return re
    
    def seqUpload(self,tBox, multiple):
        ftypes=[('Fasta files','*.fasta'),('Text files','*.txt')]        
        fname = filedialog.askopenfilename(title='Select file', filetypes=ftypes)
        if fname != '':
            text=''
            if 'fasta' in fname:
                for dnaSeq in SeqIO.parse(fname,'fasta'):
                    temp=(dnaSeq.seq)
                    text=text+' '+temp
            else:
                temp= self.readFile(fname)
                text=text+' '+temp
            if not multiple:
                tBox.delete('1.0',END)
            tBox.insert(END, text)
            
    def delBox(self,box):
        box.configure(state=NORMAL)
        box.delete('1.0',END)
        box.configure(state=DISABLED)

    def readFile(self, filename):
        f = open(filename, 'r')
        text = f.read()
        return text  

    def getDir(self, tBox):
        ftypes=[('CSV files', '*.csv'),('TSV files', '*.tsv')]
        fname=filedialog.askopenfilename(title='Select file',filetypes=ftypes)
        tBox.configure(state=NORMAL)
        tBox.delete('1.0',END)
        tBox.insert(END,fname)
        tBox.configure(state=DISABLED)
    
    #Get info from text boxes to run gel.Viz.    
    def runGV(self, geneBox,reIndex,reData):
        seq=geneBox.get('1.0',END)
        seq = str(seq) 
        allre=self.getREs(reData)
        re=''
        for i in reIndex:
            re=re+str(allre[i])+' '
        gv.gel_visualize(seq,re,reData)

    
    #Get info from text boxes to run PlasBUILDR    
    def runBUILDR(self, geneBox,reBox,reData):
        seq=geneBox.get('1.0',END)
        re=reBox.get('1.0',END)
        seq=geneBox.get('1.0',END)
        seq = str(seq) 
        re=reBox.get('1.0',END)
        re=str(re)
        pb.plasmid_builder(seq, re,reData)

    #Get info from text boxes to run SeqProp    
    def runSeqProp(self, seqBox, reData):
        seq=seqBox.get('1.0', END)
        seq=str(seq)
        sp.multSeqProp(seq, reData)

    def runMSM(self,fileBox, motifBox, threshBox):
        fname=str(self.resource_path(fileBox.get('1.0',END)))
        fname=fname.strip()
        motif=str(motifBox.get('1.0',END))
        thresh=str(threshBox.get('1.0',END))
        ''' if motif!='' and thresh!='': '''
        motif=int(motif.strip())
        thresh=int(thresh.strip()) 
        print(motif, thresh)
        msm.msm(fname,motif_length=motif,match_threshold=thresh)
        ''' elif motif!='':
            motif=int(motif.strip())
            msm.msm(fname,motif_length=motif)
        elif thresh!='':
            thresh=int(thresh.strip()) 
            msm.msm(fname,match_threshold=thresh)
        else:
            msm.msm(fname) '''
            
    #Get info from text boxes to run Primer Design    
    def runPrim(self,seqBox, idx, gcBox, tmBox, lenBox):
        seq=seqBox.get('1.0',END)
        seq=str(seq)
        #seq_idx will be a string. First value is start, last value is end. Separated by '-'
        seq_idx=idx.get('1.0',END)
        seq_idx=str(seq_idx)
        gc=gcBox.get('1.0',END)
        gc=str(gc)
        tm=tmBox.get('1.0',END)
        tm=str(tm)
        length=lenBox.get('1.0',END)
        length=str(length)
        #pd.pd(seq, start, end, gc=gc, tm=tm, range=len_range)

#Runs the program
def main():
  
    root = Tk()
    root.geometry('350x300+300+300')
    
    app = Example(root)
    
    root.protocol('WM_DELETE_WINDOW', app.quit)

    app.mainloop()


if __name__ == '__main__':
    main() 

''' comingSoon=Label(msmPg,text='A tool for returning all exact DNA matches between a list of sequences in a specified read frame.\n\nComing Soon.',width=110)
comingSoon.pack(fill=BOTH,padx=10,pady=10) '''