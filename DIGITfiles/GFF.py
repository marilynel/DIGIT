'''
GFF.py contains two classes: GffLine for parsing individual lines from Blast or individual Query objects and
stringifying data in GFF format, and GffFile for reading JSON or tab input files and creating GFF outputfiles.

These classes follow gff-version 3.1.26 format and interface specifically with the JSON files from DIGIT as well as
Blast output in tab format.

GFF documentation:
http://useast.ensembl.org/info/website/upload/gff3.html

Blast documentation (tab):
https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
'''

from Query import *


class GffLine:
    def __init__(self):
        self.seqid = None
        self.source = "DIGIT"
        self.type = "nucleotide_match"
        self.start = 0
        self.end = 0
        self.score = 0
        self.strand = None
        self.phase = "."
        self.attributes = None

    def __checkStartAndEnd__(self):
        if self.start > self.end:
            temp = self.end
            self.end = self.start
            self.start = temp

    def __fromBlastTab__(self, line):
        self.seqid = line[1]
        self.source = "blast"
        self.type = "nucleotide_match"
        self.start = int(line[8])
        self.end = int(line[9])
        self.score = float(line[10])
        self.strand = self.__setStrand__(int(line[8]), int(line[9]))
        self.phase = "."
        self.attributes = f"ID={line[0]}"

        self.__checkStartAndEnd__()

    def __fromQuery__(self, query):
        self.seqid = query.chromosome
        self.source = "DIGIT"
        self.type = "nucleotide_match"
        self.start = query.sStart
        self.end = query.sEnd
        self.score = query.eValue
        self.strand = self.__setStrand__(query.sStart, query.sEnd)
        self.phase = "."
        self.attributes = f"ID={query.query};genome={query.genome}"

        self.__checkStartAndEnd__()

    def __setStrand__(self, num1, num2):
        if num1 > num2:
            return "-"
        else:
            return "+"

    def __stringGffLine__(self):
        return f"{self.seqid}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{self.attributes}\n"


class GffFile:
    def __init__(self, filename):
        self.listGffLines = []
        self.filename = filename

        if filename.find("tab") != -1:
            self.__readFromBlast__()

        if filename.find("json") != -1:
            self.__readFromJson__()

    def __readFromBlast__(self):
        with open(self.filename, "r") as bfile:
            for line in bfile:
                if line[0] != '#':
                    lineparts = line.split("\t")
                    gffLine = GffLine()
                    gffLine.__fromBlastTab__(lineparts)
                    self.listGffLines.append(gffLine)

    def __formatFromQueriesWorkingSet__(self, queriesWorkingSet):
        for query in queriesWorkingSet:
            gffLine = GffLine()
            gffLine.__fromQuery__(queriesWorkingSet[query])
            self.listGffLines.append(gffLine)

    def __readFromJson__(self):
        # This function needs to be updated
        pass
        # queriesWorkingSet = {}
        # createQueryStruct(self.filename, queriesWorkingSet)
        # self.__formatFromQueriesWorkingSet__(queriesWorkingSet)

    def __writeToGffFile__(self, filename):
        with open(filename, "w") as gfile:
            gfile.write(f"##gff-version 3.1.26\n")
            for gff in self.listGffLines:
                gfile.write(gff.__stringGffLine__())




