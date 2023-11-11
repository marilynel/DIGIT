'''
BlastOutput is a class that parses and stores data generated from a Blast search. The Blast output must be in tab-
delimited format (-outfmt 7). Results are stored as Query objects in a struct with a group ID, path to an output
directory, a list of all the searched for queries, and a dictionary containing every hit found within the Blast
search(es) of interest.

For DIGIT, we are focusing on sequences Blasted against maize reference genomes A188v1, B73v5, and W22v2. This class is
used to parse the results of Blasting both initial flanking sequences and sanger sequencing results. The hits dictionary
within BlastOutput will be structered thusly within this framework:

    self.hits = {
        "A188v1" : {
            allele1 : [hit1_A, hit2_A, etc.],
            allele2 : [hit3_A, hit4_A, hit5_A, etc.],
            etc.
        },
        "B73v5" : {
            allele1 : [hit1_B, hit2_B, etc.],
            allele2 : [hit3_B, hit4_B, hit5_B, etc.],
            etc.
        },
        "W22v2" : {
            allele1 : [hit1_W, hit2_W, etc.],
            allele2 : [hit3_W, hit4_W, hit5_W, etc.],
            etc.
        }
    }
'''

import os

from Query import *


class BlastOutput:
    def __init__(self, groupID, groupType):
        self.hits = {}
        self.groupID = groupID
        self.outputDir = f"DIGIToutput/{groupType}/{groupID}"
        self.listQueries = []
        self.bYearFamilyDict = {}
        self.__lookupFamilyAllele__()


    def __fillBlastOutputObject__(self, b73only, primers):
        # Driver to populate BlastOutputObject and sort data to find "best" hit; writes data to JSON
        self.__parseBlastOutputFiles__(b73only, primers)
        self.__makeListOfQueries__()
        self.__pickGenome__()
        self.__allBlastOutputDataToJson__()


    def __goThroughFiles__(self, filepath):
        filename = filepath.split("/")[-1]
        genome = filename.split("_")[0]
        if genome not in self.hits:
            self.hits[genome] = {}
            self.__readBlastFile__(genome, filepath)


    def __parseBlastOutputFiles__(self, b73only, primers):
        # Finds relevant blast .tab files.
        for subdir, dirs, files in os.walk(self.outputDir):
            for oneFile in files:
                filepath = os.path.join(subdir, oneFile)
                if filepath.find(".tab") != -1:
                    # Assess query blast results (!primers) in all genomes (!b73only)
                    if not b73only and not primers:
                        if filepath.find("Primer") == -1:
                            self.__goThroughFiles__(filepath)

                    # Assess query blast results (!primers) just in B73 (b73only)
                    if b73only and not primers:
                        if filepath.find("Primer") == -1 and filepath.find("B73") != -1:
                            self.__goThroughFiles__(filepath)

                    # Assess primer blast results (primers) in all genomes (!b73only)
                    if not b73only and primers:
                        if filepath.find("Primer") != -1:
                            self.__goThroughFiles__(filepath)

                    # Assess primer blast results (primers) just in B73 (b73only)
                    if b73only and primers:
                        if filepath.find("Primer") != -1 and filepath.find("B73") != -1:
                            self.__goThroughFiles__(filepath)


    def __lookupFamilyAllele__(self):
        # Reads file and creates dictionary of family : allele associations for reference. There may be multiple
        # families that represent one allele.
        with open("DIGITfiles/YearB2022FlankingSequenceFamilyData", "r") as familyfile:
            familyfile.readline()
            for line in familyfile:
                allele, family = line.split("\t")
                if allele not in self.bYearFamilyDict:
                    self.bYearFamilyDict[allele] = []
                self.bYearFamilyDict[allele].append(int(family.strip()))


    def __readBlastFile__(self, genome, filename):
        # Parses blast .tab files for data relevant to project.
        numHits = -1
        qName = ""
        with open(filename, "r+") as blastfile:
            for line in blastfile:
                newQuery = None
                if line[0] == '#':
                    blastData = line.split(' ')
                    if blastData[1] == "Query:":
                        qName =  blastData[2].strip()
                    if blastData[2] == "hits":
                        numHits = int(blastData[1])
                    if numHits == 0:
                        newQuery = Query("", genome, numHits)
                        newQuery.query = qName
                        newQuery.__findBYearFamilies__(self.bYearFamilyDict)
                        numHits = -1
                else:
                    newQuery = Query(line, genome, numHits)
                    newQuery.__setValues__()
                    newQuery.__findBYearFamilies__(self.bYearFamilyDict)

                if newQuery:
                    if newQuery.query not in self.hits[genome]:
                        self.hits[genome][newQuery.query] = []
                    self.hits[genome][newQuery.query].append(newQuery)


    def __makeListOfQueries__(self):
        # List of query names as strings. Don't remember the purpose but I'm sure it's important.
        listQueries = []
        for genome in self.hits:
            for allele in self.hits[genome]:
                if allele not in self.listQueries:
                    self.listQueries.append(allele)
                self.__setBestQuery__(genome, allele)
        return listQueries


    def __setBestQuery__(self, genome, query):
        # Identifies the best query hit of a genome and calculates the percent difference in bit score between the best
        # and second best hit.
        bestQuery = self.hits[genome][query][0]
        secondBestQuery = None
        for i in range(1, len(self.hits[genome][query])):
            if self.hits[genome][query][i].bitScore > bestQuery.bitScore:
                secondBestQuery = bestQuery
                bestQuery = self.hits[genome][query][i]
            elif self.hits[genome][query][i].bitScore <= bestQuery.bitScore:
                if not secondBestQuery or self.hits[genome][query][i].bitScore >= secondBestQuery.bitScore:
                    secondBestQuery = self.hits[genome][query][i]

        perDiff = -1
        if secondBestQuery:
            try:
                perDiff = abs((((bestQuery.bitScore - secondBestQuery.bitScore)/(bestQuery.bitScore + secondBestQuery.bitScore)) / 2) * 100)
            except:
                perDiff = -1

        for i in range(0, len(self.hits[genome][query])):
            if self.hits[genome][query][i] == bestQuery:
                self.hits[genome][query][i].bestAlleleForGenome = True
                self.hits[genome][query][i].percentDiff = perDiff


    def __pickGenome__(self):
        # Finds the prefferred genome for an allele and identifies the hit that will be used in the working set later on
        # (bestHitForAllele). Best bit score is chosen first, then results are prioritied in genome order B73, W22, then
        # A188.
        bestDict = {}
        for genome in self.hits:
            for allele in self.hits[genome]:
                if allele not in bestDict:
                    bestDict[allele] = {}
                if genome not in bestDict[allele]:
                    bestDict[allele][genome] = -1
                for i in range(len(self.hits[genome][allele])):
                    if self.hits[genome][allele][i].bestAlleleForGenome:
                        bestDict[allele][genome] = self.hits[genome][allele][i].bitScore

        for allele in bestDict:
            maxBitScore = max(bestDict[allele].values())
            if bestDict[allele]["B73v5"] == maxBitScore:
                for i in range(len(self.hits["B73v5"][allele])):
                    if self.hits["B73v5"][allele][i].bitScore == maxBitScore:
                        self.hits["B73v5"][allele][i].bestHitForAllele = True
            elif bestDict[allele]["W22v2"] == maxBitScore:
                for i in range(len(self.hits["W22v2"][allele])):
                    if self.hits["W22v2"][allele][i].bitScore == maxBitScore:
                        self.hits["W22v2"][allele][i].bestHitForAllele = True
            elif bestDict[allele]["A188v1"] == maxBitScore:
                for i in range(len(self.hits["A188v1"][allele])):
                    if self.hits["A188v1"][allele][i].bitScore == maxBitScore:
                        self.hits["A188v1"][allele][i].bestHitForAllele = True
            else:
                print("ok")


    def __allBlastOutputDataToJson__(self):
        # Writes contents of BlastOutput object to JSON file.
        with open(f"{self.outputDir}/QueryData/JSONfiles/AllBlastHits_{self.groupID}.json", "w") as outfile:
            outfile.write(self.__toJSON__())
        return


    def __toJSON__(self):
        # For converting BlastOutput object to json object.
        json0bj = json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
        return json0bj