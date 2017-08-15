from tkinter import *
import os

def saveData(out_text, appName):
    outPg=Toplevel()
    outPg.title(appName)

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
    saved=Toplevel()
    saved.title("Save confirm")
    saveLabel=Label(saved,text="Saved!",font=("Helvetica", 16))
    saveLabel.pack(fill=BOTH)

def getDir(box):
    box.configure(state=NORMAL)
    box.delete('1.0',END)
    pathname=filedialog.askdirectory()
    pathname=pathname+'/'
    box.insert(END,pathname)
    box.configure(state=DISABLED)