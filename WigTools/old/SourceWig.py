class header:
    def __init__(self, chrom, start, step, lines):
        self.lines=lines
        self.chrom=chrom
        self.start=start
        self.stop=len(lines)*step+start
        self.step=step
        self.currentline=0
        
    def getLine(self, pos):
        return self.lines[self.getLinePosition(pos)]

    def getLinePosition(self, pos):
        return (pos-self.start)/self.step
    
class SourceWig:
    def __init__(self, filename):
        # initializes the root member
        self.filename = filename
        self.currentHeader=None
        self.currentPosition=None
        self.filePipe=open(self.filename)
        self.nextline=self.filePipe.readline()
        if "metadata: " in self.nextline:
            self.metadata=self.nextline.split("metadata: ", 1)[1]
        else:
            self.metadata=self.filename
        self.is_file_open=True

    def walk(self):
        if not self.is_file_open: return None
        if self.currentHeader==None or self.currentPosition==self.currentHeader.stop:
            if not self.nextHeader(): return None
        result=self.currentHeader.getLine(self.currentPosition).strip()
        self.currentPosition+=self.currentHeader.step
        return result, self.currentHeader.chrom, self.currentPosition-self.currentHeader.step

    def nextHeader(self):
        while self.nextline != "":
            if "fixedStep" in self.nextline:
                self.parseHeaderAndLoadLines()
                return True
            else: 
                self.nextline=self.filePipe.readline()
        self.filePipe.close()
        self.is_file_open=False
        return False
        
    def parseHeaderAndLoadLines(self):
        chrom, start, step=map(lambda x: x.split("=")[1].strip().strip("chr"), self.nextline.split(" ")[1:4])
        lines=[]
        self.nextline=self.filePipe.readline()
        while "fixedStep" not in self.nextline and self.nextline !="":
            lines.append(self.nextline)
            self.nextline=self.filePipe.readline()
        self.currentHeader=header(chrom, int(start), int(step), lines)
        self.currentPosition=self.currentHeader.start

class SourceWigInMemory:
    def __init__(self, filename):
        # initializes the root member
        self.filename = filename
        self.currentHeader=None
        self.currentPosition=None
        f=open(self.filename)
        self.lines=f.read().split("\n")
        f.close()
        self.numLines=len(self.lines)
        self.currentLine=0
        self.nextline=self.lines[self.currentLine]
        if "metadata: " in self.nextline:
            self.metadata=self.nextline.split("metadata: ", 1)[1]
        else:
            self.metadata=self.filename
        self.currentHeaderLine=None

    def walk(self):
        if self.currentLine>self.numLines: return None
        if self.currentHeader==None or self.currentPosition==self.currentHeader.stop:
            if not self.nextHeader(): return None
        result=self.getLine().strip()
        #print self.currentPosition, self.currentHeader.stop
        toreturn=(result, self.currentHeader.chrom, self.currentPosition)
        self.currentPosition+=self.currentHeader.step
        return toreturn

    def getLine(self):
        offset=self.currentHeader.getLinePosition(self.currentPosition)
        return self.lines[self.currentHeaderLine+offset+1]

    def nextHeader(self):
        while self.nextline != "":
            if "fixedStep" in self.nextline:
                self.parseHeaderAndLoadLines()
                return True
            else: 
                self.getNextLine()
        return False
        
    def getNextLine(self):
        if self.currentLine+1<self.numLines:
            self.currentLine+=1
            self.nextline=self.lines[self.currentLine]
        else:
            self.nextline=""
    
    def parseHeaderAndLoadLines(self):
        chrom, start, step=map(lambda x: x.split("=")[1].strip().strip("chr"), self.nextline.split(" ")[1:4])
        self.currentHeaderLine=self.currentLine
        lines=[]
        self.getNextLine()
        i=0
        while "fixedStep" not in self.nextline and "track" not in self.nextline and self.nextline !="":
            #print self.nextline
            self.getNextLine()
            i+=1
        #print "i!!", i
        self.currentHeader=header(chrom, int(start), int(step), lines)
        self.currentHeader.stop=self.currentHeader.start+i*self.currentHeader.step
        self.currentPosition=self.currentHeader.start

class  SourceWigWithChrom:
    def __init__(self, filename, chrom):
        # initializes the root member
        self.filename = filename
        self.currentHeader=None
        self.currentPosition=None
        self.filePipe=open(self.filename)
        self.nextline=self.filePipe.readline()
        self.is_file_open=True
        self.chrom=chrom
        if "metadata: " in self.nextline:
            self.metadata=self.nextline.split("metadata: ", 1)[1]
        else:
            self.metadata=self.filename

    def walk(self):
        if not self.is_file_open: return None
        if self.currentHeader==None:
            if not self.nextHeader(): return None
        while self.currentHeader.chrom!=self.chrom:
            if not self.nextHeader(): return None
        if self.currentPosition==self.currentHeader.stop:
            if not self.nextHeader(): return None
        result=self.currentHeader.getLine(self.currentPosition).strip()
        self.currentPosition+=1
        return result, self.currentHeader.chrom, self.currentPosition

    def nextHeader(self):
        while self.nextline != "":
            if "fixedStep" in self.nextline:
                if self.parseHeaderAndLoadLines():
                    return True
            else: 
                self.nextline=self.filePipe.readline()
        self.filePipe.close()
        self.is_file_open=False
        return False
        
    def parseHeaderAndLoadLines(self):
        chrom, start, step=map(lambda x: x.split("=")[1].strip().strip("chr"), self.nextline.split(" ")[1:4])
        if chrom!=self.chrom:
            self.nextline=self.filePipe.readline()
            return False
        lines=[]
        self.nextline=self.filePipe.readline()
        while "fixedStep" not in self.nextline and self.nextline !="":
            if "track type" not in self.nextline:
                lines.append(self.nextline)
            self.nextline=self.filePipe.readline()
        self.currentHeader=header(chrom, int(start), int(step), lines)
        self.currentPosition=self.currentHeader.start
        return True

if __name__ == "__main__":
    #x=SourceWigWithChrom("/Users/cmelton/Documents/Classes/FromDropbox/BMI217/Project/Ch22Results/CoverageFiles/0ab8d063-62b4-4d47-82aa-e3351a60029d.coverage.subset.wig.txt", "22")
    #x=SourceWigWithChrom("/Users/cmelton/Documents/Lab/SnyderLab/LocalMutationRate/WigFiles/coverage.wig.corrected.txt.test", "1")
    x=SourceWigInMemory("/Users/cmelton/Documents/Lab/SnyderLab/LocalMutationRate/WigFiles/coverage.wig.corrected.txt.test")
    next=x.walk()
    while next != None:
        print next
        next=x.walk()