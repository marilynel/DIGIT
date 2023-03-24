from Sequences import Sequences
from Query import Query


class Primer3Object:
    def __init__(self, task, side, q):
        # side refers to the side of the primer we are looking to find
        self.s = Sequences()
        self.productSizeStr = ""
        self.bobby = ""
        self.startingPrimer = ""
        self.task = task
        self.side = side
        # self.primerName = ""

        # self.__initializePrimerName__(q)
        self.__initializeVars__(q)
        self.inputStr = ""

        self.__buildInputStr__(q)

    def __buildInputStr__(self, q):
        if self.side == "left":
            self.inputStr = (f"SEQUENCE_ID={q.primerNameLeft}\n")
        elif self.side == "right":
            self.inputStr = (f"SEQUENCE_ID={q.primerNameRight}\n")
        else:
            print("Error in primer3 input init: side error")

        if self.task == "check_primers":
            self.inputStr += (
                    f"SEQUENCE_TEMPLATE={q.wildtypeSequence}\n" + \
                    f"SEQUENCE_PRIMER={q.primerSequenceLeft}" + \
                    f"SEQUENCE_PRIMER_REVCOMP={q.primerSequenceRight}\n"
                # self.productSizeStr = "300-1100"
            )

        elif self.task == "generic":
            self.inputStr += (
                f"SEQUENCE_TEMPLATE={q.insertionSequence}\n"
            )

            if self.side == "left":
                self.inputStr += (
                    f"SEQUENCE_PRIMER_REVCOMP="
                )
            elif self.side == "right":
                self.inputStr += (
                    f"SEQUENCE_PRIMER="
                )
            else:
                self.inputStr += "__error__"
                print(
                    f"Error in composing Primer3 initial input: side error at {q.query}"
                )

            '''if (q.strand == 1 and side == "left") or (q.strand == -1 and side == "right"):
                self.inputStr += f"{s.dsgg3}\n"
                productSizeStr = s.dsgg3Size
                bobby = s.dsgg3Hairpin
            elif (q.strand == -1 and side == "left") or (q.strand == 1 and side == "right"):
                self.inputStr += f"{s.gfp3utr}\n"
                productSizeStr = s.gfp3utrSize
                bobby = s.gfp3utrHairpin
            else:
                self.inputStr += "__error__"
                print(
                    f"Error in composing Primer3 initial input: side or strand error at {q.query}"
                )'''
            self.inputStr += self.startingPrimer

        else:
            self.inputStr += "__error__"
            print(
                f"Error in composing Primer3 inital input: task error at {q.query}"
            )

        self.inputStr += (
                f"\nPRIMER_MASK_KMERLIST_PATH=genomes/kmer_lists/zea_mays\n" + \
                f"PRIMER_TASK={self.task}\n" + \
                f"PRIMER_PICK_LEFT_PRIMER=1\n" + \
                f"PRIMER_PICK_INTERNAL_OLIGO=0\n" + \
                f"PRIMER_PICK_RIGHT_PRIMER=1\n" + \
                f"PRIMER_MASK_FAILURE_RATE=0.1\n" + \
                f"PRIMER_PRODUCT_SIZE_RANGE={self.productSizeStr}\n" + \
                f"PRIMER_OPT_SIZE=22\n" + \
                f"PRIMER_MIN_SIZE=20\n" + \
                f"PRIMER_MAX_SIZE=24\n" + \
                f"PRIMER_OPT_TM=62.0\n" + \
                f"PRIMER_MIN_TM=59.0\n" + \
                f"PRIMER_MAX_TM=65.0\n" + \
                f"PRIMER_TM_FORMULA=1\n" + \
                f"PRIMER_PAIR_MAX_DIFF_TM=5.0\n"
        )

        if self.task == "generic":
            self.inputStr += (
                f"PRIMER_MAX_HAIRPIN_TH={self.bobby}\n"
            )

        self.inputStr += (
                f"PRIMER_EXPLAIN_FLAG=1\n" + \
                f"PRIMER_MIN_GC=30.0\n" + \
                f"PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1\n" + \
                f"PRIMER_MAX_END_STABILITY=9.0\n" + \
                f"PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=3\n" + \
                f"PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3\n" + \
                f"PRIMER_LIBERAL_BASE=1\n" + \
                f"PRIMER_FIRST_BASE_INDEX=1\n" + \
                f"PRIMER_MAX_TEMPLATE_MISPRIMING=12.00\n" + \
                f"PRIMER_MAX_TEMPLATE_MISPRIMING_TH=47.00\n" + \
                f"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00\n" + \
                f"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=47.00\n" + \
                f"=\n"
        )

    def __initializeVars__(self, q):
        if (q.strand == 1 and self.side == "left") or (q.strand == -1 and self.side == "right"):
            self.startingPrimer = self.s.dsgg3
            self.productSizeStr = self.s.dsgg3Size
            self.bobby = self.s.dsgg3Hairpin
        elif (q.strand == -1 and self.side == "left") or (q.strand == 1 and self.side == "right"):
            self.startingPrimer = self.s.gfp3utr
            self.productSizeStr = self.s.gfp3utrSize
            self.bobby = self.s.gfp3utrHairpin
        else:
            print("Error in initializing vars in Primer3Object")

        if self.task == "check_primers":
            self.productSizeStr = "300-1100"

    def __initializePrimerName__(self, q):
        if self.side == "left":
            self.primerName = q.primer_name_left
        if self.side == "right":
            self.primerName = q.primer_name_right

    # def __parsePrimer3Output__(self):
