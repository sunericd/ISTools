#07/08/2017 allows multiple files/sequences to be uploaded.
#07/13/2017 Tab view for tools. Additional changes to reflect structure from Exe Toolkit.

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

class Example(Frame):   
  
    def __init__(self, parent):
        Frame.__init__(self, parent)   
         
        self.parent = parent
        self.initUI()
                
    def initUI(self):
        self.parent.title("Integrated Sciences Toolkit")
        self.style=Style()
        self.style.theme_use("default")  

        self.pack(fill=BOTH, expand=1)
        
        mlabel = Label(self, text="Integrated Sciences Toolkit\n\nDo NOT close this window.", font=("Helvetica", 16), wraplength=350)
        mlabel.pack(pady=15)

        #Navigation
        sButton=Button(self,text="Start",command=self.nav)
        sButton.pack(padx=15,pady=5)

        qButton=Button(self,text="Quit",command=self.quit)
        qButton.pack(padx=15,pady=5)

    # Needed for PyInstaller to read files: filepath = resource_path()

    def resource_path(self, relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        try:
            # PyInstaller creates a temp folder and stores path in _MEIPASS
            base_path = sys._MEIPASS
        except Exception:
            base_path = os.path.abspath(".")

        return os.path.join(base_path, relative_path)

    def quit(self):
        quit()

    def nav(self):
        top=Toplevel()
        top.title("Integrated Sciences Toolkit")

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
        pic =PhotoImage(file=self.resource_path("data/logo.gif"))
        pic.zoom(10,10)
        mainPg_l1=Label(mainPg,image=pic)
        mainPg_l1.image=pic
        mainPg_l1.configure(image=pic)
        mainPg_l1.grid(row=0, column=0,padx=5, sticky=E+W+S+N)

        mainPg_l2=Label(mainPg,text="Welcome to the Integrated Sciences Toolkit. Click on the tabs above to alternate between tools. \n\nSend any questions to support@integratedsciences.org\n\nVersion 1.0.0\n\nLast updated July 27, 2017.", font=("Helvetica", 12), wraplength=250)
        mainPg_l2.grid(row=0, column=1, padx=5, sticky=E+W+S+N)
        
        #!!!Add a label, maybe beginning instructions.

        ##############################Gel.Viz Page#################################
        #Gene sequence input
        gelPg_f1 = Frame(gelPg)
        gelPg_f1.pack(fill=X)
        
        gelPg_l1 = Label(gelPg_f1, text="Sequence (DNA bases only)")
        gelPg_l1.pack(side=LEFT, padx=5, pady=5) 

        gelPg_f2=Frame(gelPg)
        gelPg_f2.pack(fill=X)
        
        gelPg_geneBox = Text(gelPg_f2, height=10)        
        gelPg_geneBox.pack(fill=X, padx=5)
        
        #Restriction Enzyme list
        gelPg_f3= Frame(gelPg)
        gelPg_f3.pack(fill=X)
        
        gelPg_l2 = Label(gelPg_f2, text="Restriction Enzyme(s)")
        gelPg_l2.pack(side=LEFT, padx=5, pady=5)        

        allREs=StringVar()
        allREs.set(self.getREs(reData))

        gelPg_reList = Listbox(gelPg_f3,selectmode="multiple",listvariable=allREs,width=110) 
        gelPg_reList.pack(padx=5, side=LEFT, fill=Y)

        yScrollRE=Scrollbar(gelPg_f3)
        yScrollRE.pack(side=LEFT,fill=Y)
        yScrollRE.configure(command=gelPg_reList.yview)

        gelPg_reList.configure(yscrollcommand=yScrollRE.set)
        
        #Navigation
        gelPg_clr=Button(gelPg,text="Clear",command=lambda:[f() for f in [gelPg_reList.selection_clear(0, END),gelPg_geneBox.delete('1.0', END)]])
        gelPg_clr.pack(side=RIGHT,padx=5,pady=5)
        
        gelPg_ent= Button(gelPg,text="Enter", command=lambda:self.runGV(gelPg_geneBox,gelPg_reList.curselection(),reData))
        gelPg_ent.pack(side=RIGHT, padx=5, pady=5)  

        #Upload Gene Sequence File Button.
        gelPg_upload=Button(gelPg_f1,text="Upload Gene Sequence(s)", command=lambda: self.fileUpload(gelPg_geneBox,True))
        gelPg_upload.pack(fill=Y,expand=True,side=RIGHT, padx=5, pady=5)  

        ####################################Plasmid BUILDR##################################

        #Plasmid sequence input
        plasPg_f1 = Frame(plasPg)
        plasPg_f1.pack(fill=X)
        
        plasPg_l1 = Label(plasPg_f1, text="One Recipient Plasmid Sequence (DNA bases only)")
        plasPg_l1.pack(side=LEFT, padx=5, pady=5)           
       
        plasPg_f2=Frame(plasPg)
        plasPg_f2.pack(fill=X)
        
        plasPg_geneBox = Text(plasPg_f2, height=10)
        plasPg_geneBox.pack(fill=X, padx=5)
        
        #Restriction enzymes input
        plasPg_f3 = Frame(plasPg)
        plasPg_f3.pack(fill=X)
        
        plasPg_l2 = Label(plasPg_f2, text="Markers or Genes (Separate sequences with spaces)")
        plasPg_l2.pack(side=LEFT, padx=5, pady=5)        

        plasPg_mrkr = Text(plasPg_f3,height=10)
        plasPg_mrkr.pack(fill=X, padx=5)

        #Navigation
        plasPg_clr=Button(plasPg,text="Clear",command=lambda:[f() for f in [plasPg_mrkr.delete('1.0', END),plasPg_geneBox.delete('1.0', END)
]])
        plasPg_clr.pack(side=RIGHT,padx=5,pady=5)
        
        plasPg_ent= Button(plasPg,text="Enter", command=lambda:self.runBUILDR(plasPg_geneBox,plasPg_mrkr,reData))
        plasPg_ent.pack(side=RIGHT)

        plasPg_upl1=Button(plasPg_f1,text="Upload Recipient Plasmid Sequence", command=lambda: self.fileUpload(plasPg_geneBox,False))
        plasPg_upl1.pack(fill=Y,expand=True,side=RIGHT, pady=5)  

        plasPg_upl2=Button(plasPg_f2,text="Upload Marker/Gene Sequence(s)", command=lambda: self.fileUpload(plasPg_mrkr,True))
        plasPg_upl2.pack(fill=Y,expand=True,side=RIGHT,pady=5)

        ##############################MSM####################################

        comingSoon=Label(msmPg,text="A tool for returning all exact DNA matches between a list of sequences in a specified read frame.\n\nComing Soon.",width=110)
        comingSoon.pack(fill=BOTH,padx=10,pady=10)

        #############################SeqProp##################################

        seqPg_f1 = Frame(seqPg)
        seqPg_f1.pack(fill=X)

        seqPg_l1 = Label(seqPg_f1, text="Sequence to analyze (DNA bases only)")
        seqPg_l1.pack(side=LEFT, padx=5, pady=5)      

        seqPg_f2=Frame(seqPg)
        seqPg_f2.pack(fill=X)
        
        seqPg_box = Text(seqPg_f2, height=10)
        seqPg_box.pack(fill=X, padx=5, pady=5)
        
        seqPg_ul1 = Button(seqPg_f1, text= "Upload Sequence", command = lambda: self.fileUpload(seqPg_box, True))     
        seqPg_ul1.pack(fill=Y, pady=5, padx=5, side=RIGHT)       
        
        seqPg_clr=Button(seqPg,text="Clear",command=lambda:[f() for f in [seqPg_box.delete('1.0', END)]])
        seqPg_clr.pack(side=RIGHT,padx=5,pady=5)
        
        seqPg_ent= Button(seqPg,text="Enter", command=lambda:self.runSeqProp(seqPg_box,reData))
        seqPg_ent.pack(side=RIGHT)

        #############################PrimerDesign##################################

        primPg_f1 = Frame(primPg)
        primPg_f1.pack(fill=X)

        primPg_l1 = Label(primPg_f1, text="Sequence to amplify (DNA bases only)")
        primPg_l1.pack(side=LEFT, padx=5, pady=5)      

        primPg_f2=Frame(primPg)
        primPg_f2.pack(fill=X)
        
        primPg_box = Text(primPg_f2, height=10)
        seqPg_box.pack(fill=X, padx=5, pady=5)
        
        primPg_ul1 = Button(primPg_f1, text= "Upload Sequence", command = lambda: self.fileUpload(primPg_box, True))     
        primPg_ul1.pack(fill=Y, pady=5, padx=5, side=RIGHT)       
        
        primPg_clr=Button(primPg,text="Clear",command=lambda:[f() for f in [primPg_box.delete('1.0', END)]])
        primPg_clr.pack(side=RIGHT,padx=5,pady=5)
        
        primPg_ent= Button(primPg,text="Enter", command=lambda:self.runPrim(primPg_box))
        primPg_ent.pack(side=RIGHT)

        ##########################Credits#################################
        credPg_f1 = Frame(credPg)
        credPg_f1.pack(fill=X)    
        
        cred = Text(credPg_f1,wrap=WORD)
        cred.insert(END, "The Integrated Sciences Group is dedicated to the development of free, open-source tools for a range of scientific research. We ship large toolkits with user inteferfaces in order to make the application of our tools as seamless as possible. We are currently working to improve our current tools and to build new ones!\n\ngel.Viz is a gel visualization tool that generates an image of a gel given a sequence and restriction enzymes. This tool is the perfect resource to plan and verify experiments involving gel electrophoresis.\n\nThe Plasmid BUILDR, or Benchmark Utility for Integrating Long DNA Reads, generates a protocol to build a desired plasmid given the sequences of the recepient plasmid and any insert genes or markers.\n\nDevelopers: Yein Christina Park and Eric Sun.\n\nVisit our website at integratedsciences.org.\n\nQuestions? Concerns? Email us at intsci.tools@gmail.com.\n\nMany thanks to Yi Chen and Siavash Zamirpour.\n\n\nCopyright 2017 Integrated Sciences\n\nPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\nThe above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.")
        cred.pack(fill=BOTH, padx=5,side=LEFT)
        
        yScrollCredits=Scrollbar(credPg_f1);
        yScrollCredits.pack(side=LEFT,fill=Y)
        yScrollCredits.config(command=cred.yview)
        cred.config(state=DISABLED,yscrollcommand=yScrollCredits.set)

        #############################Compile notebook############

        ntbk.add(mainPg, compound=LEFT, text="Welcome",padding=5)
        ntbk.add(gelPg,text="gel.Viz", padding=5)
        ntbk.add(msmPg,text="Multiple Sequence Mapper",padding=5)
        ntbk.add(plasPg,text="Plasmid BUILDR", padding=5)
        ntbk.add(primPg,text="PrimerDesign",padding=5)
        ntbk.add(seqPg,text="SeqProp",padding=5)
        ntbk.add(credPg,text="Credits", padding=5)
        
        ntbk.pack(fill=BOTH)
        self.pack(fill=BOTH,expand=True)        

    def getREs(self,reData):
        re = reData.values.T.tolist()
        re = re[0]
        return re
    
    def fileUpload(self,tBox, multiple):
        ftypes=[('Fasta files','*.fasta'),('Text files','*.txt')]        
        fname = filedialog.askopenfilename(title="Select file", filetypes=ftypes)
        if fname != '':
            text=""
            if 'fasta' in fname:
                for dnaSeq in SeqIO.parse(fname,"fasta"):
                    temp=(dnaSeq.seq)
                    text=text+" "+temp
            else:
                temp= self.readFile(fname)
                text=text+" "+temp
            if not multiple:
                tBox.delete('1.0',END)
            tBox.insert(END, text)
            

    def readFile(self, filename):
        f = open(filename, "r")
        text = f.read()
        return text  
        
    
    #Get info from text boxes to run gel.Viz.    
    def runGV(self, geneBox,reIndex,reData):
        seq=geneBox.get("1.0",END)
        seq = str(seq) 
        allre=self.getREs(reData)
        re=''
        for i in reIndex:
            re=re+str(allre[i])+' '
        gv.gel_visualize(seq,re,reData)

    
    #Get info from text boxes to run PlasBUILDR    
    def runBUILDR(self, geneBox,reBox,reData):
        seq=geneBox.get("1.0",END)
        re=reBox.get("1.0",END)
        seq=geneBox.get("1.0",END)
        seq = str(seq) 
        re=reBox.get("1.0",END)
        re=str(re)
        pb.plasmid_builder(seq, re,reData)

    #Get info from text boxes to run SeqProp    
    def runSeqProp(self, seqBox, reData):
        seq=seqBox.get("1.0", END)
        seq=str(seq)
        sp.multSeqProp(seq, reData)

    #Get info from text boxes to run PlasBUILDR    
    def runPrim(self,seqBox):
        print("hi")

#Runs the program
def main():
  
    root = Tk()
    root.geometry("350x300+300+300")
    
    app = Example(root)
    
    root.protocol("WM_DELETE_WINDOW", app.quit)

    app.mainloop()


if __name__ == '__main__':
    main() 