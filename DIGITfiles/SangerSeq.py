import re

from Query import *


# Not currently in use

class BlastLine:
    def __init__(self, line):
        # # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

        # R96C12	chr5	100.000	417	0	0	1	417	132692185	132692601	0.0	771
        lineParts = line.split("\t")
        self.query = lineParts[0]
        self.chromosome = lineParts[1]
        self.perIdentity = float(lineParts[2])
        self.alignmentLength = int(lineParts[3])
        self.mismatches = int(lineParts[4])
        self.gapOpens = int(lineParts[5])
        self.qStart = int(lineParts[6])  # first bp of query
        self.qEnd = int(lineParts[7])  # last bp of query
        self.sStart = int(lineParts[8])  # first bp of query in chromosome
        self.sEnd = int(lineParts[9])  # last bp of query in chromosome
        self.eValue = lineParts[10]
        self.bitScore = float(lineParts[11])  # used to rank alleles


class SangerSeq:
    def __init__(self, query, order):
        self.query = query
        self.orderNum = order
        self.primerMatch = None
        self.sequence = None
        self.query = None
        self.genomeDict = {
            "A188v1": [],
            "B73v5": [],
            "W22v2": []
        }

    def __makePrimerMatch__(self, filename):
        # 28329_C02_R179H03_5-Ds-3_C02_006.seq
        self.orderNum = filename[:5]
        self.query = re.findall(r"R\d{1,4}[A-Z]*\d{1,4}", filename)
        if filename.find("DsGFP_3UTR") != -1:
            self.primerMatch = "DsGFP_3UTR"
        if filename.find("5-Ds-3") != -1:
            self.primerMatch = "5-Ds-3"

    def __addBlastLine__(self, genome, line):
        newBlastLine = BlastLine(line)
        self.genomeDict[genome].append(newBlastLine)

    def __numHitsPerGenome__(self, genome):
        return len(self.genomeDict[genome])

    def __bestHit__(self, genome):
        pass
