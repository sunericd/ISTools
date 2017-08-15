#8/15/17 Added/started GUI for ProtTK

from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
from PIL import Image, ImageTk
from Bio import SeqIO
import mods.protQuant as pq
import os
import pandas as pd

class Example(Frame):   
  
    def __init__(self, parent):
        Frame.__init__(self, parent)   
         
        self.parent = parent
        self.initUI()
                
    def initUI(self):
        self.parent.title('Integrated Sciences Protein Toolkit')
        self.style=Style()
        self.style.theme_use('default')  

        self.pack(fill=BOTH, expand=1)
        
        mlabel = Label(self, text='Integrated Sciences Protein Toolkit\n\nDo NOT close this window.', font=('Helvetica', 16),justify=CENTER,wraplength=350)
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
        top.title('Integrated Sciences Protein Toolkit')

        #######################ALL Page Declarations###################

        ntbk=Notebook(top)
        mainPg=Frame(ntbk)
        credPg=Frame(ntbk)
        quantPg=Frame(ntbk)

        #########################Front page#################################3
        
        #Insert personalized logo
        pic = PhotoImage(file=self.resource_path('data/logo.gif'))
        pic.zoom(10,10)
        mainPg_l1=Label(mainPg,image=pic)
        mainPg_l1.image=pic
        mainPg_l1.configure(image=pic)
        mainPg_l1.grid(row=0, column=0,padx=5, sticky=E+W+S+N)

        mainPg_l2=Label(mainPg,text='Welcome to the Integrated Sciences Protein Toolkit. Click on the tabs above to alternate between tools. \n\nSend any questions to support@integratedsciences.org\n\nVersion 1.0.0\n\nLast updated Aug 15, 2017.', font=('Helvetica', 12), wraplength=250)
        mainPg_l2.grid(row=0, column=1, padx=5, sticky=E+W+S+N)
        
        #!!!Add a label, maybe beginning instructions.

        ##############################ProtQuant Page##############################
        quantPg_l1 = Label(quantPg, text='Select a .csv file of the standards.')
        quantPg_l1.grid(row=0,column=0, padx=5, pady=5,sticky=W,columnspan=2)

        quantPg_l2 = Label(quantPg, text='Standard CSV File: ')
        quantPg_l2.grid(row=1,column=0, padx=5, pady=5,sticky=W)

        quantPg_browse=Button(quantPg, text='Browse', command=lambda: self.getDir(quantPg_std))
        quantPg_browse.grid(row=1,column=1, padx=5, pady=5,sticky=E)

        quantPg_std = Text(quantPg, height=1, state=DISABLED)
        quantPg_std.grid(row=2,column=0, columnspan=2, padx=5, pady=5,sticky=N+S+W)
        
        quantPg_l1 = Label(quantPg, text='Select a .csv file of the protein fluorescence data.')
        quantPg_l1.grid(row=3,column=0, padx=5, pady=5,sticky=W,columnspan=2)

        quantPg_l2 = Label(quantPg, text='Protein CSV File: ')
        quantPg_l2.grid(row=4,column=0, padx=5, pady=5,sticky=W)

        quantPg_browse=Button(quantPg, text='Browse', command=lambda: self.getDir(quantPg_data))
        quantPg_browse.grid(row=4,column=1, padx=5, pady=5,sticky=E)

        quantPg_data = Text(quantPg, height=1, state=DISABLED)
        quantPg_data.grid(row=5,column=0, columnspan=2, padx=5, pady=5,sticky=N+S+W)

        #Navigation
        quantPg_clr=Button(quantPg,text='Clear',command=lambda:[f() for f in [quantPg_std.delete('1.0', END),quantPg_data.delete('1.0', END)]])
        quantPg_clr.grid(row=6,column=0, padx=5, pady=5,sticky=W)
        
        quantPg_ent= Button(quantPg,text='Enter', command=lambda:self.runPQ(quantPg_std,quantPg_data))
        quantPg_ent.grid(row=6,column=1, padx=5, pady=5,sticky=E)

        ##########################Credits#################################        
        cred = Text(credPg,wrap=WORD)
        cred.insert(END, 'The Integrated Sciences Group is dedicated to the development of free, open-source tools for a range of scientific research. We ship large toolkits with user inteferfaces in order to make the application of our tools as seamless as possible. We are currently working to improve our current tools and to build new ones!\n\nDevelopers: Yein Christina Park, Eric Sun, and Yi Chen.\n\nVisit our website at integratedsciences.org.\n\nQuestions? Concerns? Email us at support@integratedsciences.org.\n\nMany thanks to Jimmy Thai and Siavash Zamirpour.\n\n\nCopyright 2017 Integrated Sciences\n\nPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \'Software\'), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\nThe above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \'AS IS\', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.')
        cred.grid(row=0,column=0, pady=5,sticky=E+W+S+N)
        
        yScrollCredits=Scrollbar(credPg)
        yScrollCredits.grid(row=0,column=1, pady=5,sticky=E+N+S)
        yScrollCredits.config(command=cred.yview)
        cred.config(state=DISABLED,yscrollcommand=yScrollCredits.set)

        #############################Compile notebook############

        ntbk.add(mainPg, compound=LEFT, text='Welcome',padding=5)
        ntbk.add(quantPg,text='ProtQuant', padding=5)
        ntbk.add(credPg,text='Credits', padding=5)        
        
        ntbk.pack(fill=BOTH)
        self.pack(fill=BOTH,expand=True)

    def getDir(self, tBox):
        ftypes=[('CSV files', '*.csv'),('TSV files','*.tsv')]
        fname=filedialog.askopenfilename(title='Select file',filetypes=ftypes)
        tBox.configure(state=NORMAL)
        tBox.delete('1.0',END)
        tBox.insert(END,fname)
        tBox.configure(state=DISABLED)
    
    #Get info from text boxes to run protQuant.    
    def runPQ(self,standards,proteins):
        stdfile=str(self.resource_path(standards.get('1.0',END)))
        stdfile=stdfile.strip()
        datafile=str(self.resource_path(proteins.get('1.0',END)))
        datafile=datafile.strip()
        pq.protQuant(stdfile,datafile)

#Runs the program
def main():
  
    root = Tk()
    root.geometry('350x300+300+300')
    
    app = Example(root)
    
    root.protocol('WM_DELETE_WINDOW', app.quit)

    app.mainloop()


if __name__ == '__main__':
    main() 

''' comingSoon=Label(quantPg,text='A tool for returning all exact DNA matches between a list of sequences in a specified read frame.\n\nComing Soon.',width=110)
comingSoon.pack(fill=BOTH,padx=10,pady=10) '''