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
    def __init__(self, seqType):
        self.seqid = None
        self.source = "DIGIT"
        self.type = seqType
        self.start = 0
        self.end = 0
        self.score = 0
        self.strand = None
        self.phase = "."
        self.attributes = None


    def __checkStartAndEnd__(self):
        # GFF files require the lower position number to be the start, even if that is where the sequence ends. Used
        # internally.
        if self.start > self.end:
            temp = self.end
            self.end = self.start
            self.start = temp


    def __fromQuery__(self, query, seqType):
        # Converts Query object to GffLine object. Used in GffFile class.
        self.seqid = query.chromosome
        self.source = "DIGIT"
        self.type = seqType
        self.start = query.sStart
        self.end = query.sEnd
        self.score = query.eValue
        self.strand = self.__setStrand__(query.sStart, query.sEnd)
        self.phase = "."
        self.attributes = f"ID={query.query};genome={query.genome}"
        self.__checkStartAndEnd__()


    def __setStrand__(self, num1, num2):
        # Finds if allele is on plus or minus strand. Used internally.
        if num1 > num2:
            return "-"
        else:
            return "+"


    def __stringGffLine__(self):
        # Return GffLine in format compatible with GFF3 file formatting. Used in GffFile class.
        return f"{self.seqid}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t" + \
            f"{self.phase}\t{self.attributes}\n"


########################################################################################################################


class GffFile:
    def __init__(self):
        self.listGffLines = []


    def __formatFromQueriesWorkingSet__(self, workingSet, seqType):
        # Creates a GFF file from the selected working set of queries. Used in QueriesWorkingSet class.
        for query in workingSet:
            if workingSet[query]:
                gffLine = GffLine(seqType)
                gffLine.__fromQuery__(workingSet[query], seqType)
                self.listGffLines.append(gffLine)


    def __writeToGffFile__(self, filename):
        # Write data to GFF file. Used in QueriesWorkingSet class.
        with open(filename, "w") as gfile:
            gfile.write(f"##gff-version 3.1.26\n")
            for gff in self.listGffLines:
                gfile.write(gff.__stringGffLine__())




