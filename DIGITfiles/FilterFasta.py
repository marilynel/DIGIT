'''
FilterFasta.py finds multiple genomic sequences in a genome at given coordinates.

Written by: Marilyn Leary 2023

    Input:
    filename        Path to fasta file of reference genome (many can be found at
                    https://www.ncbi.nlm.nih.gov/data-hub/genome/)
    coordinates     A dictionary containing the information to look up in reference genome. Format as follows:

                    coordinates = {
                        alleleName : [
                            chromosomeID, [startingPosition, endingPosition]
                        ],
                        "R01C02_wt" : [
                            "chr5", [66710080, 66714080]
                        ]
                    }

                    As many allele coordinate sets as needed may be listed in this format.

    Output:
    sequenceData    Dictionary containing sequences for each allele listed in coordinates. Format as follows:

                    sequenceData = {
                        alleleName : sequence01
                        "R01C02_wt" : "TCAAGGGCACCAGCCGTCCTGCTCACTACCATGTCTTGTGGGACGAGAACAACTTCACAGCCGACGCACTGCAGACCCTCA
                                        CAACAACCTTTGCTACACGTAAGCTAGCTGCTCACAAAAAGAGGTGTCAGTGTCAGTTCAGTTCCCTGAACACCGACCGT
                                        TAATATAATAGCTTGTCAAATTGCCGCTGCAGCTACGCGAGGTGCACGCGCTCTGTGTCCATTGGTAGGTTTGTCAAAGT
                                        TATCCATTGCTGAATCGATGCACGAGCAACAGATAAACCAAAAGATGCATATGCCTTCCTTGCAGTCCCGCCGGCGTACT
                                        GCTCACCTGGCCGCATTCCGCGCCCGGTTCTACATGGAGCCTGACAGCTCAGACAGCGGGTCGCTGGCGAGTGGCGCCT"
                    }

                    All of the alleles given in the coordinate set will be in this returned struct.

'''


# import sys

class FilterFasta:
    def __init__(self, coordinates):
        self.dictAllCoordinates = {}  # keys = chr, values = list of all coordinates needed on that chromosome
        self.dictMinCoordinates = {}  # keys = chr, values = lowest int from chr list in dictAllCoordinates
        self.dictMaxCoordinates = {}  # keys = chr, values = greatest int from chr list in dictAllCoordinates
        self.genomeStruct = {}  # keys = chr, values = genomic sequence between max and min from above
        self.listChr = []  # list of chromosomes data will need to be pulled from

        self.__setVals__(coordinates)
        self.__setMinMax__()

    def __setVals__(self, coordinates):
        for allele in coordinates:
            if coordinates[allele][0] not in self.listChr:
                self.listChr.append(coordinates[allele][0])

            if coordinates[allele][0] not in self.dictAllCoordinates:
                self.dictAllCoordinates[coordinates[allele][0]] = []
            self.dictAllCoordinates[coordinates[allele][0]].append(coordinates[allele][1][0])
            self.dictAllCoordinates[coordinates[allele][0]].append(coordinates[allele][1][1])

    def __setMinMax__(self):
        '''
        Finds minimum and maximum coordinates for each chromosome. Used to limit how much data is pulled from genome
        file.
        '''
        for chr in self.dictAllCoordinates:
            self.dictMinCoordinates[chr] = min(self.dictAllCoordinates[chr]) - 1
            self.dictMaxCoordinates[chr] = max(self.dictAllCoordinates[chr])


def filterFasta(filename, coordinates):
    ffobj = FilterFasta(coordinates)

    with open(filename, "r") as zfile:
        key, chrSeq = "", ""
        for line in zfile:
            if line[0] == ">":
                if chrSeq and (key in ffobj.listChr):
                    ffobj.genomeStruct[key] = chrSeq[ffobj.dictMinCoordinates[key]:
                                                     ffobj.dictMaxCoordinates[key]]
                chrSeq, key = "", ""
                key = line[1:].strip()
            else:
                chrSeq += line.strip()

    sequenceData = {}
    for c in coordinates:
        chr = coordinates[c][0]
        c1 = coordinates[c][1][0] - ffobj.dictMinCoordinates[chr] - 1
        c2 = coordinates[c][1][1] - ffobj.dictMinCoordinates[chr]
        seq = ffobj.genomeStruct[chr][c1:c2]
        sequenceData[c] = seq

    return sequenceData

