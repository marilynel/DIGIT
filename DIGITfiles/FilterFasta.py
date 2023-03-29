import sys

'''
filename = "DIGITfiles/Genomes/Zm-A188-REFERENCE-KSU-1.0.fa"

coordinates = {
    "allele1" : ["chr1", [500, 600]],
    "allele2" : ["chr1", [550, 650]]
}

listChr = ["chr1"]
'''


class FilterFasta:
    def __init__(self, coordinates):
        self.dictOfCoordinates = {}  # keys = allele, values = list where first val is the chromosome and second is a list of the coordinates
        self.dictAllCoordinates = {}  # keys = chromosome, values = list of all coordinates needed on that chromosome
        self.dictMinCoordinates = {}  # keys = chromosome, values =
        self.dictMaxCoordinates = {}
        self.genomeStruct = {}
        self.listChr = []

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
        Method finds minum and maximum coordinates for each chromosome. Used to limit how much data is pulled from genome file.
        '''
        for chr in self.dictAllCoordinates:
            self.dictMinCoordinates[chr] = min(self.dictAllCoordinates[chr])
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
                    charSeq, key = "", ""
                key = line[1:].strip()
            else:
                chrSeq += line.strip()

    sequenceData = {}
    for c in coordinates:
        chr = coordinates[c][0]
        c1, c2 = coordinates[c][1][0] - ffobj.dictMinCoordinates[chr], coordinates[c][1][1] - \
                 ffobj.dictMinCoordinates[chr]
        seq = ffobj.genomeStruct[chr][c1:c2]
        sequenceData[c] = seq

        # print(f"c: {c}")
        # print(f"c1: {c1}")
        # print(f"c2: {c2}")
        # print(f"chr: {chr}")
        # print(f"full section of genome: {genomeStruct[chr]}")
        # print(f"taget section of genome: {seq}")

        # print(f"sequenceData[c]: {sequenceData[c]}")
    return sequenceData


def main():
    filterFasta(filename, coordinates, listChr)


if __name__ == "__main__":
    main()