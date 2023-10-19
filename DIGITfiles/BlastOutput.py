import os
from Query import *


class BlastOutput:
    def __init__(self, groupID, groupType):
        self.hits = {}
        self.groupID = groupID
        self.outputDir = f"DIGIToutput/{groupType}/{groupID}"
        self.listQueries = []

    def __parseBlastOutputFiles__(self, genomeSpec):
        for subdir, dirs, files in os.walk(self.outputDir):
            for oneFile in files:
                filepath = os.path.join(subdir, oneFile)
                if filepath.find(".tab") != -1 and filepath.find(genomeSpec) != -1:
                    filename = filepath.split("/")[-1]
                    genome = filename.split("_")[0]
                    if genome not in self.hits:
                        self.hits[genome] = {}
                    self.__readBlastFile__(genome, filepath)
                    # allQueries[genome] = processBlastOutput(filepath,genome)

    def __readBlastFile__(self, genome, filename):
        # bYearFamilies = lookupFamilyAllele()
        # genomeDB is self.hits[genome]
        numHits = -1
        qName = ""
        with open(filename, "r+") as blastfile:
            for line in blastfile:
                newQuery = None
                if line[0] == '#':
                    blastData = line.split(' ')
                    if blastData[1] == "Query:":
                        qName = blastData[2].strip()
                    if blastData[2] == "hits":
                        numHits = int(blastData[1])
                    if numHits == 0:
                        newQuery = Query("", genome, numHits)
                        newQuery.query = qName
                        # newQuery.__findBYearFamilies__(bYearFamilies)
                        # if newQuery.query not in genomeDB:
                        #    genomeDB[newQuery.query] = []
                        # genomeDB[newQuery.query].append(newQuery)
                else:
                    newQuery = Query(line, genome, numHits)
                    newQuery.__setValues__()
                    # TODO: implement find b year family stuff
                    # newQuery.__findBYearFamilies__(bYearFamilies)
                    # if newQuery.query not in genomeDB:
                    #    genomeDB[newQuery.query] = []
                    # genomeDB[newQuery.query].append(newQuery)
                if newQuery:
                    if newQuery.query not in self.hits[genome]:
                        self.hits[genome][newQuery.query] = []
                    self.hits[genome][newQuery.query].append(newQuery)

    def __findBestForGenome__(self):

        listQueries = []
        for genome in self.hits:
            for allele in self.hits[genome]:
                if allele not in self.listQueries:
                    self.listQueries.append(allele)
                # getBestQuery(genome, allele)
                self.__setBestQuery__(genome, allele)
        return listQueries

    def __setBestQuery__(self, genome, query):
        bestQuery = self.hits[genome][query][0]
        secondBestQuery = None
        for i in range(1, len(self.hits[genome][query])):
            if self.hits[genome][query][i].bitScore > bestQuery.bitScore:
                secondBestQuery = bestQuery
                bestQuery = self.hits[genome][query][i]
            elif self.hits[genome][query][i].bitScore <= bestQuery.bitScore:
                if not secondBestQuery or self.hits[genome][query][
                    i].bitScore >= secondBestQuery.bitScore:
                    secondBestQuery = self.hits[genome][query][i]

        perDiff = -1
        if secondBestQuery:
            try:
                perDiff = abs((((bestQuery.bitScore - secondBestQuery.bitScore) / (
                            bestQuery.bitScore + secondBestQuery.bitScore)) / 2) * 100)
            except:
                perDiff = None

        for i in range(0, len(self.hits[genome][query])):
            if self.hits[genome][query][i] == bestQuery:
                self.hits[genome][query][i].bestAlleleForGenome = True
                self.hits[genome][query][i].percentDiff = perDiff

    # TODO: how do I find the working set genome for an allele?
    def __pickGenome__(self):
        bestDict = {}
        # for each genome in hits:
        for genome in self.hits:
            # for each allele in hits[genome]:
            for allele in self.hits[genome]:
                if allele not in bestDict:
                    bestDict[allele] = {}
                if genome not in bestDict[allele]:
                    bestDict[allele][genome] = -1
                # for each hit in hits[genome][allele]:
                for i in range(len(self.hits[genome][allele])):
                    # if .bestAlleleForGenome is True:
                    if self.hits[genome][allele][i].bestAlleleForGenome:
                        # need bit score for that genome
                        bestDict[allele][genome] = self.hits[genome][allele][i].bitScore

        # ultimately something like:
        # {
        #    allele1 : {genome1 : bitscore1, genome2: bitscore2, genome3 : bitscore3},
        #    allele2 : {genome1: bitscore4, genome2 : biscore5, genome3 : bitscore6},
        #    etc.
        # }

        # then:
        # for each allele
        for allele in bestDict:
            # find max(struct[allele].value()) to get max bit score

            maxBitScore = max(bestDict[allele].values())
            # if struct[allele][B73] == max:
            if bestDict[allele]["B73v5"] == maxBitScore:
                # for hit in hits[B73][allele]:
                for i in range(len(self.hits["B73v5"][allele])):
                    # if hits[genome][allele][hit].bitscore = max:
                    if self.hits["B73v5"][allele][i].bitScore == maxBitScore:
                        # hits[genome][allele][hit].bestHitForAllele = True
                        self.hits["B73v5"][allele][i].bestHitForAllele = True
            elif bestDict[allele]["W22v2"] == maxBitScore:
                # for hit in hits[B73][allele]:
                for i in range(len(self.hits["W22v2"][allele])):
                    # if hits[genome][allele][hit].bitscore = max:
                    if self.hits["W22v2"][allele][i].bitScore == maxBitScore:
                        # hits[genome][allele][hit].bestHitForAllele = True
                        self.hits["W22v2"][allele][i].bestHitForAllele = True
            elif bestDict[allele]["A188v1"] == maxBitScore:
                # for hit in hits[B73][allele]:
                for i in range(len(self.hits["A188v1"][allele])):
                    # if hits[genome][allele][hit].bitscore = max:
                    if self.hits["A188v1"][allele][i].bitScore == maxBitScore:
                        # hits[genome][allele][hit].bestHitForAllele = True
                        self.hits["A188v1"][allele][i].bestHitForAllele = True
            else:
                print("ok")

    def __allBlastOutputDataToJson__(self):
        with open(f"{self.outputDir}/JSONfiles/AllSangerHits_{self.groupID}.json", "w") as outfile:
            outfile.write(self.__toJSON__())
        return

    def __toJSON__(self):
        json0bj = json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
        return json0bj