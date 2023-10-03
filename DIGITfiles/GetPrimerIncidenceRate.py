import sys

from Utils import *
from Query import *
from Primer import *


def main():
    allPrimers = {}
    flankseq = sys.argv[1]
    genomeSpec = "PrimerBlast.tab"
    dirname = f"DIGIToutput/FlankingSequences/{flankseq}/BlastOutput/PrimerBlast"

    findBlastOutputFiles(dirname, allPrimers, genomeSpec)
    # Create dict of Primer objects; key=query (R num), val=list of Primer objects
    primerDict = makePrimerDict(allPrimers)

    workingQueriesJsonFile = f"DIGIToutput/FlankingSequences/" \
                             f"{flankseq}/QueryData/JSONfiles/WorkingSet_{flankseq}.json"
    queriesWorkingSet = QueriesWorkingSet()
    queriesWorkingSet.__createQueryStructFromJson__(workingQueriesJsonFile)

    # Check if primers are in the right spot
    setPrimerFoundNearFlankSeq(primerDict, queriesWorkingSet)
    writePrimerIncidenceDataToFile(primerDict, flankseq)


if __name__ == '__main__':
    main()