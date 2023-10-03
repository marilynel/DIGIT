import json

from Sequences import Sequences

'''
Plus strand:

-2000                               [12345678]                            +2000                
-------------------------------------------------------------------------------
                                    ^ s. start                          ^ s.end

first group: lower position indices
    [(s.start-2000), (s.start+7)]
second group: higher position indices
    [s.start, (s.start+2000)]
third group: wild type
    [(s.start-2000), (s.start+2000)]


Minus strand (numbered relative to plus strand):

-2000                              [87654321]                             +2000                
-------------------------------------------------------------------------------
                      ^ s.end              ^ s. start                  


first group: lower position indices
    [(s.start-2000), (s.start)]
second group: higher position indices
    [(s.start-7), (s.start+2000)]
third group: wild type
    [(s.start-2000), (s.start+2000)]
'''


class Query:
    def __init__(self, line, genome, num_hits):
        if line:
            items = line.split('\t')
        else:
            items = [i == 0 for i in range(12)]

        # Data from Blast output files
        self.query = items[0]  # allele / R number / query
        self.chromosome = items[1]  # chr location of allele instance
        self.perIdentity = items[2]
        self.alignmentLength = items[3]
        self.mismatches = items[4]
        self.gapOpens = items[5]
        self.qStart = int(items[6])  # first bp of query
        self.qEnd = int(items[7])  # last bp of query
        self.sStart = int(items[8])  # first bp of query in chromosome
        self.sEnd = int(items[9])  # last bp of query in chromosome
        self.eValue = items[10]
        self.bitScore = float(items[11])  # used to rank alleles
        self.genome = genome  # genome database where allele instance was found
        self.numHits = num_hits  # total number of instances for query in genome
        ### NEW: changed default diff to -1, do not delete this comment until I can check the new
        # outputs against old!!! 5/22 ###
        self.percentDiff = -1
        self.strand = 0  # plus or minus
        self.qStartStatus = False  # true if sStart = 1
        self.bitScoreStatus = False  # true if bitScore >= 80

        # Data determined based on strand direction
        self.primerNameLeft = None
        self.primerNameRight = None

        # Kept for personal ease of understanding
        self.sideMatchGFP3UTR = None  # right or left primer matches c GFP3UTR
        self.sideMatch3DsgG = None  # right or left primer matches c 3DsgG

        # Data from filterfasta output files, or calc'd from Blast output
        self.upperCoordinates = [-1, -1]  # around sStart to +2000 bp
        self.upperSequence = None
        self.lowerCoordinates = [-1, -1]  # -2000 bp to around sStart
        self.lowerSequence = None
        self.insertionSequence = None  # predicted DsGFP insertion sequence for allele
        self.wildtypeCoordinates = [-1, -1]  # -2000 bp to +2000 bp relative to sStart
        self.wildtypeSequence = None  # full allele without DsGFP insertion

        # Data from Primer3 output files
        self.primerSequenceLeft = "FAIL"
        self.primerSequenceRight = "FAIL"  # from primer3 output

        self.primerLeftProductSize = -1  # from generic primer3 output
        self.primerRightProductSize = -1  # from generic primer3 output
        self.primerPairProductSize = -1  # from WT verfication
        self.tmLeft = -1  # from Generic P3 output
        self.tmRight = -1  # from Generic P3 output
        self.tmPair = -1  # from Validate P3 output
        self.primerPenaltyLeft = -1  # from generic primer3 output
        self.primerPenaltyRight = -1  # from generic primer3 output
        self.primerPairPenalty = -1  # from WT verification

        # Set during comparison of alleles
        self.bestAlleleForGenome = False  # true if this allele has best bit score across all
        # instances in genome
        self.bestHitForAllele = False  # true if this instance in this genome is the best for all
        # versions of the allele

        self.bYearFamilies = None

    def __findBYearFamilies__(self, familyDict):
        if self.query in familyDict:
            self.bYearFamilies = familyDict[self.query]

    def __setValues__(self):
        self.__strandDirection__()
        self.__setqStartStatus__()
        self.__setBitScoreStatus__()
        self.__getUpperCoordinates__()
        self.__getLowerCoordinates__()
        self.__getWildtypeCoordinates__()

    def __makeBestGenomeString__(self):
        # For printing to best_queries_by_genome.csv
        return f"{self.bitScore},{self.numHits},{self.qStartStatus},"

    def __getUpperCoordinates__(self):
        if self.strand == 1:
            self.upperCoordinates = [self.sStart, self.sStart + 2000]
        elif self.strand == -1:
            self.upperCoordinates = [self.sStart - 7, self.sStart + 2000]
        else:
            self.upperCoordinates = [-1, -1]

    def __getLowerCoordinates__(self):
        if self.strand == 1:
            self.lowerCoordinates = [self.sStart - 2000, self.sStart + 7]
        elif self.strand == -1:
            self.lowerCoordinates = [self.sStart - 2000, self.sStart]
        else:
            self.lowerCoordinates = [-1, -1]

    def __getWildtypeCoordinates__(self):
        self.wildtypeCoordinates = [self.sStart - 2000, self.sStart + 2000]

    def __setBitScoreStatus__(self):
        if self.bitScore >= 80:
            self.bitScoreStatus = True

    def __strandDirection__(self):
        '''
        Determine if sequence is on the plus or minus strand and set corresponding variables
        '''
        if int(self.sStart < self.sEnd):
            # plus strand
            self.sideMatchGFP3UTR = "right"  # gfp3utr is a left primer, match with right
            self.sideMatch3DsgG = "left"  # 3dsgg is a right primer, match with left
            self.primerNameLeft = "b_" + self.query  # left primer matches iwth 3dsgg
            self.primerNameRight = "a_" + self.query  # right primer matches with gfp3utr
            self.strand = 1
        else:
            # minus strand
            self.sideMatchGFP3UTR = "left"  # gfp3utr is a right primer, match with left
            self.sideMatch3DsgG = "right"  # 3dsgg is a left primer, match with right
            self.primerNameLeft = "a_" + self.query  # left primer matches with gfp3utr
            self.primerNameRight = "b_" + self.query  # right primer matches witj 3dsgg
            self.strand = -1

    def __setqStartStatus__(self):
        if self.qStart == 1:
            self.qStartStatus = True

    def __buildInsertionSequence__(self):
        seqObj = Sequences()
        if self.strand == 1:
            self.insertionSequence = self.lowerSequence + seqObj.dsgfp + self.upperSequence
        elif self.strand == -1:
            self.insertionSequence = self.lowerSequence + seqObj.dsgfpRevComp + self.upperSequence
        else:
            self.insertionSequence = "FAILED TO MAKE INSERTION SEQUENCE"

    def __QueryFromJSON__(self, jsonObject):
        '''
        Annoying, but it works
        '''
        self.query = jsonObject["query"]
        self.chromosome = jsonObject["chromosome"]
        self.perIdentity = float(jsonObject["perIdentity"])
        self.alignmentLength = int(jsonObject["alignmentLength"])
        self.mismatches = int(jsonObject["mismatches"])
        self.gapOpens = int(jsonObject["gapOpens"])
        self.qStart = int(jsonObject["qStart"])
        self.qEnd = int(jsonObject["qEnd"])
        self.sStart = int(jsonObject["sStart"])
        self.sEnd = int(jsonObject["sEnd"])
        self.eValue = jsonObject["eValue"]
        self.bitScore = float(jsonObject["bitScore"])
        self.genome = jsonObject["genome"]
        self.numHits = int(jsonObject["numHits"])
        self.percentDiff = float(jsonObject["percentDiff"])
        self.strand = int(jsonObject["strand"])
        self.qStartStatus = jsonObject["qStartStatus"]
        self.bitScoreStatus = jsonObject["bitScoreStatus"]
        self.primerNameLeft = jsonObject["primerNameLeft"]
        self.primerNameRight = jsonObject["primerNameRight"]
        self.sideMatchGFP3UTR = jsonObject["sideMatchGFP3UTR"]
        self.sideMatch3DsgG = jsonObject["sideMatch3DsgG"]
        self.upperSequence = jsonObject["upperSequence"]
        self.lowerSequence = jsonObject["lowerSequence"]
        self.insertionSequence = jsonObject["insertionSequence"]
        self.wildtypeSequence = jsonObject["wildtypeSequence"]
        self.primerSequenceLeft = jsonObject["primerSequenceLeft"]
        self.primerSequenceRight = jsonObject["primerSequenceRight"]
        self.primerLeftProductSize = int(jsonObject["primerLeftProductSize"])
        self.primerRightProductSize = int(jsonObject["primerRightProductSize"])
        self.primerPairProductSize = int(jsonObject["primerPairProductSize"])
        self.upperCoordinates = jsonObject["upperCoordinates"]
        self.lowerCoordinates = jsonObject["lowerCoordinates"]
        self.wildtypeCoordinates = jsonObject["wildtypeCoordinates"]
        self.primerLeftProductSize = jsonObject["primerLeftProductSize"]
        self.primerRightProductSize = jsonObject["primerRightProductSize"]
        self.primerPairProductSize = jsonObject["primerPairProductSize"]
        self.tmLeft = jsonObject["tmLeft"]
        self.tmRight = jsonObject["tmRight"]
        self.tmPair = jsonObject["tmPair"]
        self.primerPenaltyLeft = jsonObject["primerPenaltyLeft"]
        self.primerPenaltyRight = jsonObject["primerPenaltyRight"]
        self.primerPairPenalty = jsonObject["primerPairPenalty"]
        self.bYearFamilies = jsonObject["bYearFamilies"]
        self.bestAlleleForGenome = jsonObject["bestAlleleForGenome"]
        self.bestHitForAllele = jsonObject["bestHitForAllele"]

    def __workingSetCsvLine__(self):
        return f"{self.query},{self.genome},{self.bYearFamilies},{self.chromosome}," + \
               f"{self.perIdentity},{self.alignmentLength},{self.mismatches},{self.gapOpens}," + \
               f"{self.qStart},{self.qEnd},{self.sStart},{self.sEnd},{self.eValue}," + \
               f"{self.bitScore},{self.numHits},{self.percentDiff},{self.strand}," + \
               f"{self.qStartStatus},{self.primerNameLeft},{self.primerSequenceLeft}," + \
               f"{self.primerLeftProductSize},{self.primerPenaltyLeft},{self.primerNameRight}," + \
               f"{self.primerSequenceRight},{self.primerRightProductSize}," + \
               f"{self.primerPenaltyRight},{self.primerPairProductSize},{self.primerPairPenalty},\n"

    def __validate__(self, selfObj, otherObj):
        if selfObj == otherObj:
            return True
        else:
            print(f"Warning in {self.query}: {selfObj} and {otherObj} are not equivalent.")
            return False

    def __updateQueryWithP3Output__(self, p3obj):
        from Primer3Object import Primer3Object
        if p3obj.task == "generic":
            if p3obj.primerNameLeft:
                self.primerNameLeft = p3obj.primerNameLeft
                self.primerSequenceLeft = p3obj.primerSequenceLeft
                self.primerLeftProductSize = p3obj.primerLeftProductSize
                self.primerPenaltyLeft = p3obj.primerPenaltyLeft
                self.tmLeft = p3obj.tmLeft

            if p3obj.primerNameRight:
                self.primerNameRight = p3obj.primerNameRight
                self.primerSequenceRight = p3obj.primerSequenceRight
                self.primerRightProductSize = p3obj.primerRightProductSize
                self.primerPenaltyRight = p3obj.primerPenaltyRight
                self.tmRight = p3obj.tmRight
        elif p3obj.task == "check_primers":
            self.primerPairProductSize = p3obj.primerPairProductSize
            self.primerPairPenalty = p3obj.primerPairPenalty
            self.tmPair = p3obj.tmPair
        else:
            print(f"Something went wrong with {self.primerNameRight}.")

    def __toJSON__(self):
        json0bj = json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
        return json0bj

    def __lookupPrint__(self):
        print()
        print(f"Allele\t\t\t{self.query}")
        print(f"Genome\t\t\t{self.genome}")
        print(f"Chromosome\t\t{self.chromosome}")
        print(f"Start\t\t\t{self.sStart}")
        print(f"End\t\t\t{self.sEnd}")
        print(f"B Year Families\t\t{self.bYearFamilies}")
        print(f"Left Primer Name\t{self.primerNameLeft}")
        print(f"Left Primer Sequence\t{self.primerSequenceLeft}")
        print(f"Right Primer Name\t{self.primerNameRight}")
        print(f"Right Primer Sequence\t{self.primerSequenceRight}")
        print()

    def __leftPrimerDataLine__(self):
        leftPrimerMatch, rightPrimerMatch = self.__primerSides__()
        return f"{self.primerNameLeft},{self.query},{self.genome},{self.primerSequenceLeft}," \
               f"{leftPrimerMatch},{self.primerLeftProductSize},{self.primerPairProductSize}," + \
               f"{self.tmLeft},{self.primerPenaltyLeft},{self.primerPairPenalty}," + \
               f"{self.bitScore},{self.numHits},{self.qStartStatus}\n"

    def __rightPrimerDataLine__(self):
        leftPrimerMatch, rightPrimerMatch = self.__primerSides__()
        return f"{self.primerNameRight},{self.query},{self.genome},{self.primerSequenceRight}," \
               f"{rightPrimerMatch},{self.primerRightProductSize},{self.primerPairProductSize}," + \
               f"{self.tmRight},{self.primerPenaltyRight},{self.primerPairPenalty}," + \
               f"{self.bitScore},{self.numHits},{self.qStartStatus}\n"

    def __primerSides__(self):
        leftPrimerMatch, rightPrimerMatch = "", ""
        if self.sideMatchGFP3UTR == "right" and self.sideMatch3DsgG == "left":
            leftPrimerMatch, rightPrimerMatch = "GFP3UTR", "3DsgG"
        elif self.sideMatchGFP3UTR == "left" and self.sideMatch3DsgG == "right":
            leftPrimerMatch, rightPrimerMatch = "3DsgG", "GFP3UTR"
        else:
            leftPrimerMatch, rightPrimerMatch = "FAIL", "FAIL"
        return leftPrimerMatch, rightPrimerMatch

    def __iter__(self):
        iters = dict((x, y) for x, y in Query.__dict__.items() if x[:2] != '__')
        iters.update(self.__dict__)
        for x, y in iters.items():
            yield x, y
