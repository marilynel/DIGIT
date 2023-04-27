import json
from Sequences import Sequences


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
        self.diff = 0  # TODO
        self.strand = 0  # plus or minus
        self.qStartStatus = False  # true if sStart = 1
        self.bitScoreStatus = False  # true if bitScore >= 80

        # Data determined based on strand direction
        self.primerNameLeft = None
        self.primerNameRight = None

        # TODO: determine if these are even needed
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

        self.primerLeftProductSize = -1  # from primer3 output
        self.primerRightProductSize = -1  # from primer3 output
        self.primerPairProductSize = -1  # TODO
        self.tm = -1  # TODO: where does this come from?
        self.primerPairPenalty = -1  # TODO: where should this come from? WT verification?
        self.primerPenaltyLeft = -1  # from primer3 output
        self.primerPenaltyRight = -1  # from primer3 output

        self.primerLeftExplain = None
        self.primerRightExplain = None
        self.primerPairExplain = None

        # Set during comparison of alleles
        self.bestAlleleForGenome = False  # true if this allele has best bit score across all
        # instances in genome
        self.bestHitForAllele = False  # true if this instance in this genome is the best for all
        # versions of the allele

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

    # TODO: write out explanation of coordinates
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
            self.primerNameLeft = self.query + "b"  # left primer matches iwth 3dsgg
            self.primerNameRight = self.query + "a"  # right primer matches with gfp3utr
            self.strand = 1
        else:
            # minus strand
            self.sideMatchGFP3UTR = "left"  # gfp3utr is a right primer, match with left
            self.sideMatch3DsgG = "right"  # 3dsgg is a left primer, match with right
            self.primerNameLeft = self.query + "a"  # left primer matches with gfp3utr
            self.primerNameRight = self.query + "b"  # right primer matches witj 3dsgg
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

    # TODO: is there a better way to do this?
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
        self.diff = float(jsonObject["diff"])
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
        self.tm = jsonObject["tm"]
        self.primerPairPenalty = jsonObject["primerPairPenalty"]
        self.primerPenaltyLeft = jsonObject["primerPenaltyLeft"]
        self.primerPenaltyRight = jsonObject["primerPenaltyRight"]
        self.bestAlleleForGenome = jsonObject["bestAlleleForGenome"]
        self.bestHitForAllele = jsonObject["bestHitForAllele"]
        self.primerLeftExplain = jsonObject["primerLeftExplain"]
        self.primerRightExplain = jsonObject["primerRightExplain"]
        self.primerPairExplain = jsonObject["primerPairExplain"]

    # TODO: get rid of this and replace with __updateQueryWithGenericP3Output__
    '''
    def __updateQueryPrimer3__(self, output_dict):
        from Primer3Object import Primer3Object
        # right primer
        if "SEQUENCE_PRIMER" in output_dict and output_dict["PRIMER_PAIR_NUM_RETURNED"] != "0":
            self.primerPenaltyRight = output_dict["PRIMER_RIGHT_0_PENALTY"]
            self.primerSequenceRight = output_dict["PRIMER_RIGHT_0_SEQUENCE"]
            #self.primerPositionRight = output_dict["PRIMER_RIGHT_0"].split(",")[0]
            self.primerRightProductSize = output_dict["PRIMER_PAIR_0_PRODUCT_SIZE"]
        #else:
        #    print(f"{self.query} has no right primer")
        if "SEQUENCE_PRIMER" in output_dict and output_dict["PRIMER_PAIR_NUM_RETURNED"] == "0":
            print(f"{self.query} has no right primer")
        # left primer
        if "SEQUENCE_PRIMER_REVCOMP" in output_dict and output_dict["PRIMER_PAIR_NUM_RETURNED"] 
        != "0":
            self.primerPenaltyLeft = output_dict["PRIMER_LEFT_0_PENALTY"]
            self.primerSequenceLeft = output_dict["PRIMER_LEFT_0_SEQUENCE"]
            #self.primerPositionLeft = output_dict["PRIMER_LEFT_0"].split(",")[0]
            self.primerLeftProductSize = output_dict["PRIMER_PAIR_0_PRODUCT_SIZE"]
        #else:
        #    print(f"{self.query} has no left primer")
        if "SEQUENCE_PRIMER_REVCOMP" in output_dict and output_dict["PRIMER_PAIR_NUM_RETURNED"] 
        == "0":
            print(f"{self.query} has no left primer")
    '''

    def __validate__(self, selfObj, otherObj):
        if selfObj == otherObj:
            return True
        else:
            print(f"Warning in {self.query}: {selfObj} and {otherObj} are not equivalent.")
            return False

    def __updateQueryWithGenericP3Output__(self, p3obj):
        from Primer3Object import Primer3Object
        if self.__validate__(self.insertionSequence, p3obj.insertionSequence):
            if p3obj.primerNameLeft:
                self.primerNameLeft = p3obj.primerNameLeft
                self.primerSequenceLeft = p3obj.primerSequenceLeft
                self.primerLeftProductSize = p3obj.primerLeftProductSize
                self.primerPenaltyLeft = p3obj.primerPenaltyLeft
                self.primerLeftExplain = p3obj.primerLeftExplain

            if p3obj.primerNameRight:
                self.primerNameRight = p3obj.primerNameRight
                self.primerSequenceRight = p3obj.primerSequenceRight
                self.primerRightProductSize = p3obj.primerRightProductSize
                self.primerPenaltyRight = p3obj.primerPenaltyRight
                self.primerRightExplain = p3obj.primerRightExplain

    def __updateQueryWithVerificationP3Output__(self):
        from Primer3Object import Primer3Object
        self.__validate__(self.wildtypeSequence, p3obj.wildtypeSequence)
        self.__validate__(self.primerNameLeft, p3obj.primerNameLeft)
        self.__validate__(self.primerSequenceLeft, p3obj.primerSequenceLeft)
        self.__validate__(self.primerLeftProductSize, p3obj.primerLeftProductSize)
        self.__validate__(self.primerPenaltyLeft, p3obj.primerPenaltyLeft)
        self.__validate__(self.primerLeftExplain, p3obj.primerLeftExplain)

        self.__validate__(self.primerNameRight, p3obj.primerNameRight)
        self.__validate__(self.primerSequenceRight, p3obj.primerSequenceRight)
        self.__validate__(self.primerRightProductSize, p3obj.primerRightProductSize)
        self.__validate__(self.primerPenaltyRight, p3obj.primerPenaltyRight)
        self.__validate__(self.primerRightExplain, p3obj.primerRightExplain)

    def __toJSON__(self):
        json0bj = json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
        return json0bj

    def __iter__(self):
        iters = dict((x, y) for x, y in Query.__dict__.items() if x[:2] != '__')
        iters.update(self.__dict__)
        for x, y in iters.items():
            yield x, y

