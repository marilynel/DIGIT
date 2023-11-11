'''
WARNING: I do not recommend changing this code AT ALL unless absolutely necessary. Please make sure you fully understand
how this class is used by DIGIT before augmenting it.
'''

'''
This function creates files structured to submit to Primer3 on CQLS. 

Written by: Marilyn Leary 2023
'''

from Sequences import Sequences


class Primer3Object:
    def __init__(self, task):

        # The following will be gathered from a given Query object in __initValsFromQueryObj__()
        # side refers to the side of the primer we are looking to find
        # Allele or query ID
        self.query = "query"

        # Allele name plus a or b
        self.primerNameLeft = None
        self.primerNameRight = None

        # Sequences of primers; changed in output functions
        self.primerSequenceLeft = "primerSequenceLeft"
        self.primerSequenceRight = "primerSequenceRight"

        self.strand = 0

        # task is generic or check_primers.
        # sequenceTemplate is DsGFP insertion sequence or wildtype sequence, depending on the task
        self.task = task
        self.insertionSequence = "insert"
        self.wildtypeSequence = "wildtype"

        # To set primers and associated data for generic input
        self.seqObj = Sequences()

        # From seqObj, depending on strand
        self.productSizeStrLeft = "productSizeStrLeft"
        self.productSizeStrRight = "productSizeStrRight"
        self.bobbyLeft = "bobbyLeft"
        self.bobbyRight = "bobbyRight"

        # Associated primer data from output
        self.primerPairNumReturned = 0
        self.primerLeftProductSize = 0
        self.primerRightProductSize = 0
        self.primerLeftExplain = "tbd"
        self.primerRightExplain = "tbd"
        self.primerPairExplain = "tbd"
        self.primerWarning = False

        self.primerPenaltyLeft = -1
        self.primerPenaltyRight = -1
        self.tmLeft = -1
        self.tmRight = -1

        # From WT verification
        self.primerPairProductSize = 0
        self.primerPairPenalty = -1
        self.tmPair = -1

        self.inputStr = ""

        self.oligoStr = (
            f"PRIMER_PICK_INTERNAL_OLIGO=0\n" +
            f"PRIMER_PICK_RIGHT_PRIMER=1\n" +
            f"PRIMER_MASK_FAILURE_RATE=0.1\n"
        )

        self.sizeStr = (
            f"PRIMER_OPT_SIZE=22\n" +
            f"PRIMER_MIN_SIZE=20\n" +
            f"PRIMER_MAX_SIZE=24\n" +
            f"PRIMER_OPT_TM=62.0\n" +
            f"PRIMER_MIN_TM=59.0\n" +
            f"PRIMER_MAX_TM=65.0\n" +
            f"PRIMER_TM_FORMULA=1\n" +
            f"PRIMER_PAIR_MAX_DIFF_TM=5.0\n"
        )

        self.lastStr = (
            f"PRIMER_EXPLAIN_FLAG=1\n" +
            f"PRIMER_MIN_GC=30.0\n" +
            f"PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1\n" +
            f"PRIMER_MAX_END_STABILITY=9.0\n" +
            f"PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=3\n" +
            f"PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3\n" +
            f"PRIMER_LIBERAL_BASE=1\n" +
            f"PRIMER_FIRST_BASE_INDEX=1\n" +
            f"PRIMER_MAX_TEMPLATE_MISPRIMING=12.00\n" +
            f"PRIMER_MAX_TEMPLATE_MISPRIMING_TH=47.00\n" +
            f"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00\n" +
            f"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=47.00\n" +
            f"=\n"
        )

    def __initValsFromQueryObj__(self, q):
        self.query = q.query
        self.primerNameLeft = q.primerNameLeft
        self.primerNameRight = q.primerNameRight
        self.strand = q.strand
        self.insertionSequence = q.insertionSequence
        self.wildtypeSequence = q.wildtypeSequence

        self.__setInputData__()

        if self.task == "check_primers":
            # Only needed for verification
            self.primerSequenceLeft = q.primerSequenceLeft
            self.primerSequenceRight = q.primerSequenceRight


    def __setInputData__(self):
        '''
        Set constant values for primer object using value of strand and task type.
        '''
        if self.task == "generic":
            # Searching for allele-specific primers to match 3DsgG and GFP3UTR (from DsGFP element)
            if self.strand == 1:
                # On the plus strand, 3DsgG will be used to find the left primer and GFP3UTR will be used to find the right primer
                self.primerSequenceRight = self.seqObj.dsgg3
                self.productSizeStrLeft = self.seqObj.dsgg3Size
                self.bobbyLeft = self.seqObj.dsgg3Hairpin

                self.primerSequenceLeft = self.seqObj.gfp3utr
                self.productSizeStrRight = self.seqObj.gfp3utrSize
                self.bobbyRight = self.seqObj.gfp3utrHairpin

            elif self.strand == -1:
                # On the minus strand, GFP3UTR will be used to find the left primer and 3DsgG will be used to find the right primer
                self.primerSequenceRight = self.seqObj.gfp3utr
                self.productSizeStrLeft = self.seqObj.gfp3utrSize
                self.bobbyLeft = self.seqObj.gfp3utrHairpin

                self.primerSequenceLeft = self.seqObj.dsgg3
                self.productSizeStrRight = self.seqObj.dsgg3Size
                self.bobbyRight = self.seqObj.dsgg3Hairpin

        elif self.task == "check_primers":
            # Verifying already gathered primers
            self.productSizeStrLeft = "300-1100"


    def __buildInputStrGeneric__(self):
        '''
        Compose begining of input string for Generic P3 input.
        '''
        self.inputStr = (
            # left primer
            f"SEQUENCE_ID={self.primerNameLeft}\n" +
            f"SEQUENCE_TEMPLATE={self.insertionSequence}\n" +
            f"SEQUENCE_PRIMER_REVCOMP={self.primerSequenceRight}\n" +
            f"PRIMER_MASK_KMERLIST_PATH=genomes/kmerLists/zea_mays\n" +
            f"PRIMER_TASK={self.task}\n" +
            f"PRIMER_PICK_LEFT_PRIMER=1\n" +
            f"{self.oligoStr}" +
            f"PRIMER_PRODUCT_SIZE_RANGE={self.productSizeStrLeft}\n" +
            f"{self.sizeStr}" +
            f"PRIMER_MAX_HAIRPIN_TH={self.bobbyLeft}\n" +
            f"{self.lastStr}" +

            # right primer
            f"SEQUENCE_ID={self.primerNameRight}\n" +
            f"SEQUENCE_TEMPLATE={self.insertionSequence}\n" +
            f"SEQUENCE_PRIMER={self.primerSequenceLeft}\n" +
            f"PRIMER_MASK_KMERLIST_PATH=genomes/kmerLists/zea_mays\n" +
            f"PRIMER_TASK={self.task}\n" +
            f"PRIMER_PICK_LEFT_PRIMER=1\n" +
            f"{self.oligoStr}" +
            f"PRIMER_PRODUCT_SIZE_RANGE={self.productSizeStrRight}\n" +
            f"{self.sizeStr}" +
            f"PRIMER_MAX_HAIRPIN_TH={self.bobbyRight}\n" +
            f"{self.lastStr}"
        )


    def __buildInputStrValidate__(self):
        '''
        Compose begining of input string for Validation P3 input.
        '''
        self.inputStr = (
            f"SEQUENCE_ID={self.query}\n" +
            f"SEQUENCE_TEMPLATE={self.wildtypeSequence}\n" +
            f"SEQUENCE_PRIMER={self.primerSequenceLeft}\n" +
            f"SEQUENCE_PRIMER_REVCOMP={self.primerSequenceRight}\n" +
            f"PRIMER_MASK_KMERLIST_PATH=genomes/kmerLists/zea_mays\n" +
            f"PRIMER_TASK={self.task}\n" +
            f"PRIMER_PICK_LEFT_PRIMER=1\n" +
            f"{self.oligoStr}" +
            f"PRIMER_PRODUCT_SIZE_RANGE={self.productSizeStrLeft}\n" +
            f"{self.sizeStr}" +
            f"{self.lastStr}"

        )

    def __initValsFromP3Output__(self, outData):
        '''
        Intialize values for a new P3 object from a dict parsed from the P3 output file.
        '''
        # outData will be one set of = to = data from p3 output, dict form
        self.query = outData["SEQUENCE_ID"][:-1]
        # The insertion sequence may be used to verify that the output data and query data match
        if "PRIMER_PAIR_NUM_RETURNED" in outData:
            self.primerPairNumReturned = int(outData["PRIMER_PAIR_NUM_RETURNED"])
        else:
            self.primerPairNumReturned = -1
        if self.primerPairNumReturned > 0:
            if self.task == "generic":
                if "SEQUENCE_PRIMER_REVCOMP" in outData:
                    # P3 output contains information on the left primer
                    self.primerNameLeft = outData["SEQUENCE_ID"]
                    self.primerSequenceLeft = outData["PRIMER_LEFT_0_SEQUENCE"]
                    self.primerPenaltyLeft = outData["PRIMER_LEFT_0_PENALTY"]
                    self.primerLeftProductSize = outData["PRIMER_PAIR_0_PRODUCT_SIZE"]
                    self.tmLeft = outData["PRIMER_PAIR_0_PRODUCT_TM"]

                if "SEQUENCE_PRIMER" in outData:
                    # P3 output contains information on the right primer
                    self.primerNameRight = outData["SEQUENCE_ID"]
                    self.primerSequenceRight = outData["PRIMER_RIGHT_0_SEQUENCE"]
                    self.primerPenaltyRight = outData["PRIMER_RIGHT_0_PENALTY"]
                    self.primerRightProductSize = outData["PRIMER_PAIR_0_PRODUCT_SIZE"]
                    self.tmRight = outData["PRIMER_PAIR_0_PRODUCT_TM"]

            if self.task == "check_primers":

                self.primerPairPenalty = outData["PRIMER_PAIR_0_PENALTY"]
                self.tmPair = outData["PRIMER_PAIR_0_PRODUCT_TM"]
                self.primerPairProductSize = outData["PRIMER_PAIR_0_PRODUCT_SIZE"]


    def __parseOutputStrValidate__(self, outData):
        self.query = outData["SEQUENCE_ID"]
        self.wildtypeSequence = outData["SEQUENCE_TEMPLATE"]
        self.primerPairNumReturned = int(outData["PRIMER_PAIR_NUM_RETURNED"])

        if self.primerPairNumReturned != 0:
            self.primerPairPenalty = outData["PRIMER_PAIR_0_PENALTY"]
            self.primerPenaltyLeft = outData["PRIMER_LEFT_0_PENALTY"]
            self.primerPenaltyRight = outData["PRIMER_RIGHT_0_PENALTY"]


    def __initializePrimerName__(self):
        if self.side == "left":
            self.primerName = self.q.primer_name_left
        if self.side == "right":
            self.primerName = self.q.primer_name_right

