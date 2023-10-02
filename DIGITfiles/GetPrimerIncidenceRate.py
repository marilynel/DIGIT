from Utils import *
from Query import *
import sys


class Primer:
    def __init__(self, queryObj):
        # mr299b_R161H09_ogRefGen_B73v5
        self.primerID = queryObj.query
        primerNameParts = queryObj.query.split("_")
        self.primerName = primerNameParts[0] + "_" + primerNameParts[1]
        self.query = primerNameParts[1]
        self.selectedGenome = primerNameParts[3]
        self.chromosome = queryObj.chromosome
        self.genome = queryObj.genome
        self.numHits = queryObj.numHits
        self.sStart = queryObj.sStart
        self.sEnd = queryObj.sEnd
        self.primerFoundNearFlankSeq = False


def writePrimerNumHitsToCSV(dirname, flankseq, numHitsDir):
    with open(dirname + flankseq + "_PrimerOccuranceInGenome.csv", "w") as csvfile:
        csvfile.write(f"PrimerName,NumHits,Genome,\n")
        for genome in numHitsDir:
            for primer in numHitsDir[genome]:
                csvfile.write(f"{primer},{numHitsDir[genome][primer]},{genome},\n")


def findNumHitsPerPrimer():
    numHitsDir = {}
    for genome in allPrimers:
        if genome not in numHitsDir:
            numHitsDir[genome] = {}
        for primer in allPrimers[genome]:
            for hit in allPrimers[genome][primer]:
                parts = hit.query.split("_")
                primerName = parts[0] + "_" + parts[1]
                if primerName not in numHitsDir[genome]:
                    numHitsDir[genome][primerName] = 0
                numHitsDir[genome][primerName] = hit.numHits
    return numHitsDir


def matchChrAndGenome(primer, query):
    if primer.chromosome != query.chromosome:
        return False
    if primer.genome != query.genome:
        return False
    return True


def checkCoordinates(primer, query):
    if matchChrAndGenome(primer, query):
        if query.strand == 1:
            if (primer.sStart >= query.wildtypeCoordinates[0]) and (
                    primer.sEnd <= query.wildtypeCoordinates[1]):
                return True
        elif query.strand == -1:
            if (primer.sStart <= query.wildtypeCoordinates[0]) and (
                    primer.sEnd >= query.wildtypeCoordinates[1]):
                return True

    return False


def makePrimerDict(allPrimers):
    primerDict = {
        "A188v1": {},
        "B73v5": {},
        "W22v2": {}
    }
    for genome in allPrimers:
        for primer in allPrimers[genome]:
            for i in range(0, len(allPrimers[genome][primer])):
                newPrimer = Primer(allPrimers[genome][primer][i])
                if newPrimer.query not in primerDict[genome]:
                    primerDict[genome][newPrimer.query] = []
                primerDict[genome][newPrimer.query].append(newPrimer)
    return primerDict


def setPrimerFoundNearFlankSeq(primerDict, queriesWorkingSet):
    for genome in primerDict:
        for query in primerDict[genome]:
            for i in range(0, len(primerDict[genome][query])):
                primerDict[genome][query][i].primerFoundNearFlankSeq = checkCoordinates(
                    primerDict[genome][query][i], queriesWorkingSet[query])


def main():
    allPrimers = {}
    queriesWorkingSet = {}
    flankseq = sys.argv[1]
    genomeSpec = "PrimerBlast.tab"
    dirname = f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/PrimerBlast"

    # Locate blast output files for the primer sequences
    findBlastOutputFiles(dirname, allPrimers, genomeSpec)
    # Create dict of Primer objects; key=query (R num), val=list of Primer objects
    primerDict = makePrimerDict(allPrimers)

    # Create working queries struct of Query objects from flankseq's working set JSON
    workingQueriesFile = f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/JSONfiles/WorkingSet_{flankseq}.json"
    createQueryStruct(workingQueriesFile, queriesWorkingSet)

    # Check if primers are in the right spot
    setPrimerFoundNearFlankSeq(primerDict, queriesWorkingSet)

    with open(f"DIGIToutput/FlankingSequences/{flankseq}/PrimerData/TBD.csv", "w") as outfile:
        for genome in primerDict:
            for query in primerDict[genome]:
                for i in range(0, len(primerDict[genome][query])):
                    outfile.write(
                        f"{primerDict[genome][query][i].primerName},{primerDict[genome][query][i].numHits},{primerDict[genome][query][i].genome},{primerDict[genome][query][i].primerFoundNearFlankSeq}\n")


if __name__ == '__main__':
    main()