class header:
    def __init__(self, chrom, start, step, lines):
        self.lines=lines
        self.chrom=chrom
        self.start=start
        self.numLines=len(lines)
        self.stop=self.numLines*step+start
        self.step=step
        self.currentline=0
        
    def getLine(self, pos):
#         print "start", self.start
#         print "stop", self.stop
#         print (pos-self.start)/self.step
        index=(pos-self.start)/self.step
        if index<0 or index>=self.numLines:
            print "getLine index wrong!!!!!", index
            return None
        return self.lines[(pos-self.start)/self.step]

class OverlapWig:
    def __init__(self, filename):
        # initializes the root member
        self.filename = filename
        self.currentHeader=None
        self.nextHeader=None
        self.currentPosition=None
        self.filePipe=open(self.filename)
        self.nextline=self.filePipe.readline()
        self.is_file_open=True
        self.i=0
        if "metadata: " in self.nextline:
            self.metadata=self.nextline.split("metadata: ", 1)[1]
        else:
            self.metadata=self.filename

    def getWigResults(self, chrom, pos):
        # if file is done return None
        if not self.is_file_open: 
            #print "file closed"
            return None
        #else:
            #print "file is open"
        # if header doesn't exist, get next header
        if self.currentHeader==None:
            if not self.getNextHeader(): return None
#         else:
#             print "current Header exists"
        # while chrom isn't right get next header
        while chrom !=self.currentHeader.chrom:
#             print "chrom isn't right"
            #print "getting new header, old header had chrom:", self.currentHeader.chrom, "looking for:", chrom
            if not self.getNextHeader(): return None
#         if self.i==0: 
#             if "3prime_overlapping_ncrna_exonInfo.wig" in self.filename or "lincRNA_exonInfo" in self.filename:
#                 print self.filename, self.currentHeader.chrom
#         self.i+=1
        # if position is below header return None
 
        # if position is between current and next header, return None
        if self.nextHeader != None:
            if pos>=self.currentHeader.stop and pos < self.nextHeader.start:
#                print "between headers"
                return None
        while pos>=self.currentHeader.stop:
#            print "getting next header"
            if not self.getNextHeader(): return None
#         print "pos:", pos, "start:", self.currentHeader.start, "stop:", self.currentHeader.stop
        if pos<self.currentHeader.start:
            return None       
        # if position is in current header
        result=self.currentHeader.getLine(pos).strip()
        return result

    def getNextHeader(self):
        if self.currentHeader==None:
            self.nextLine=" "
            while "fixedStep" not in self.nextline and self.nextline!="": 
                self.nextline=self.filePipe.readline()
        if self.nextline != "":
            self.parseHeaderAndLoadLines()
            return True
#             else: 
#                 self.nextline=self.filePipe.readline()
        self.filePipe.close()
        self.is_file_open=False
        return False
        
    def parseHeaderAndLoadLines(self):
        lines=[]
        if "fixedStep" in self.nextline:
            chrom, start, step=map(lambda x: x.split("=")[1].strip().strip("chr"), self.nextline.split(" ")[1:4])
            self.nextline=self.filePipe.readline()
        else:
            chrom, start, step = self.nextHeader.chrom, self.nextHeader.start, self.nextHeader.step
            lines=self.currentHeader.lines[int(len(self.currentHeader.lines)/2):]
        i=len(lines)
        linelimit=1000000
        while "fixedStep" not in self.nextline and self.nextline !="" and i<linelimit:
            lines.append(self.nextline)
            self.nextline=self.filePipe.readline()
            i+=1
        self.currentHeader=header(chrom, int(start), int(step), lines)
        if "fixedStep" in self.nextline:
            chrom, start, step=map(lambda x: x.split("=")[1].strip().strip("chr"), self.nextline.split(" ")[1:4])
            self.nextHeader=header(chrom, int(start), int(step), [])
        else:
            self.nextHeader=header(chrom, int(int(start)+int(step)*linelimit/2), int(step), [])
        print "current header start", self.currentHeader.start
        self.currentPosition=self.currentHeader.start

if __name__ == "__main__":
    x=OverlapWig("/Users/cmelton/Documents/Classes/FromDropbox/BMI217/Project/Ch22Results/CoverageFiles/0ab8d063-62b4-4d47-82aa-e3351a60029d.coverage.subset.wig.txt")
    print x.getWigResults('22', 16050001)
    print x.getWigResults('22', 16050002)
    print x.getWigResults('22', 16050003)
    print x.getWigResults('22', 16050004)
    print x.getWigResults('22', 16050005)
    print x.getWigResults('22', 16050006)
    print x.getWigResults('22', 16050007)
    print x.getWigResults('22', 16050008)
    print x.getWigResults('22', 16050009)
    print x.getWigResults('22', 16050010)
    print x.getWigResults('22', 16050014)
    print x.getWigResults('22', 16050015)
    #print x.currentHeader.start, x.currentHeader.step, x.currentHeader.stop
    print "x.getWigResults('22', 224196428-1)"
    print x.getWigResults('22', 224196428-1)
    #print x.currentHeader.start, x.currentHeader.step, x.currentHeader.stop
    print "\nprint x.getWigResults('22', 224196428)"
    print x.getWigResults('22', 224196428)
    print "\nprint x.getWigResults('22', 224196429)"
    print x.getWigResults('22', 224196429)
    print "\nprint x.getWigResults('22', 224196428+1886-1)"
    print x.getWigResults('22', 224196428+1886-1)
    print "\nprint x.getWigResults('22', 224196428+1886)"
    print x.getWigResults('22', 224196428+1886)
    
    #print x.currentHeader.start, x.currentHeader.step, x.currentHeader.stop
#     print x.getWigResults('22', 16050002)
#     print x.getWigResults('22', 16050003)
#     print x.getWigResults('22', 16050004)
#     print x.getWigResults('22', 16050008)
#     print x.getWigResults('22', 16050009)
#     print x.getWigResults('22', 16050010)
#     print x.is_file_open
#     print x.getWigResults('22', 160500020)
#     print x.is_file_open
#     print x.getWigResults('22', 234196428)
#     print x.getWigResults('22', 224196528)
#     print x.is_file_open