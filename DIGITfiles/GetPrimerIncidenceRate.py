'''
GetPrimerIncidenceRate.py creates a csv file out of number of times each primer occurs in each of the genomes, A188,
B73, and W22.

Written by: Marilyn Leary 2023

    Input:
    flankseq (sys.argv[1])  String representing the Flanking Sequence Set of interest. Should match with a subfolder of
                            DIGIToutput/FlankingSequences

    Output:
    Three Blast output files with path and naming convention:
        DIGIToutput/FlankingSequences/<flankseq>/PrimerData/<flankseq>_PrimerIncidenceRate.csv
'''

import sys

from Utils import *
from Query import *
from Primer import *
from BlastOutput import *


def main():
    allPrimers = {}
    flankseq = sys.argv[1]
    genomeSpec = "PrimerBlast.tab"
    dirname = f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/PrimerBlast"

    primerBlastResults = BlastOutput(flankseq, "FlankingSequences")

    # TODO: set this for b73 only
    primerBlastResults.__parseBlastOutputFiles__(False, True)

    workingQueriesJsonFile = f"DIGIToutput/FlankingSequences/{flankseq}/QueryData/JSONfiles/WorkingSet_{flankseq}.json"
    queriesWorkingSet = QueriesWorkingSet()
    queriesWorkingSet.__createQueryStructFromJson__(workingQueriesJsonFile)

    primerBlastResults.__allBlastOutputDataToJson__()

    primerDict = {}

    # TODO: move this function to QueriesWorkingSet class?
    for genome in primerBlastResults.hits:
        for query in primerBlastResults.hits[genome]:
            pNum, allele, somethingToIgnore, somethingElseToIgnore = query.split("_")
            primerName = pNum + "_" + allele
            if primerName not in primerDict:
                primerDict[primerName] = {}
            if genome not in primerDict[primerName]:
                primerDict[primerName][genome] = {"numHits":0, "anyHitInsideOgCoor": False}
            primerDict[primerName][genome]["numHits"] = primerBlastResults.hits[genome][query][0].numHits

            # Iterate through the hits for that genome
            for primer in primerBlastResults.hits[genome][query]:
                # IFF the chr in this hit is the same as the chr from our reference query; otherwise this does not need
                # to be assessed.
                if primer.chromosome == queriesWorkingSet.workingSet[allele].chromosome:
                    # Get min and max coordinate values for that particular hit
                    lowHitCoor = min(primer.sStart, primer.sEnd)
                    higHitCoor = max(primer.sStart, primer.sEnd)

                    lowHitQuer = queriesWorkingSet.workingSet[allele].wildtypeCoordinates[0]
                    higHitQuer = queriesWorkingSet.workingSet[allele].wildtypeCoordinates[1]

                    # Check if primer coordinates are inside of the wildtype sequence coordinates
                    # Purpose: make sure there is at least one instance of this primer inside the flanking sequence's
                    # wildtype sequence
                    if lowHitCoor >= lowHitQuer and higHitCoor <= higHitQuer:
                        primerDict[primerName][genome]["anyHitInsideOgCoor"] = True

    with open(f"DIGIToutput/FlankingSequences/{flankseq}/PrimerData/IncidenceRate.csv", "w") as outfile:
        outfile.write(f"PrimerID,Genome,NumHitsInGenome,PrimerOccursInPredictedWTSeq,Genome,NumHitsInGenome,PrimerOcc" +
                      f"ursInPredictedWTSeq,Genome,NumHitsInGenome,PrimerOccursInPredictedWTSeq,\n")
        for primer in primerDict:
            outfile.write(f"{primer},")
            for genome in primerDict[primer]:
                outfile.write(f"{genome},{primerDict[primer][genome]['numHits']}," +
                              f"{primerDict[primer][genome]['anyHitInsideOgCoor']},")
            outfile.write(f"\n")


if __name__ == '__main__':
    main()