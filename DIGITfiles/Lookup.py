'''
This function allows the user to look up a specific flanking sequence or primer and it's
associated data.

Requires no input arguments. Results are printed to console.

Written by: Marilyn Leary 2023
'''

from Utils import *


def cleanupAlleleName(allele):
    allele = allele.upper()
    if allele[-1] == "A" or allele[-1] == "B":
        letter = allele[-1]
        letter.lower()
        allele = allele[:-1] + letter
    return allele


def lookupFlankingSequenceData(originalData):
    print("Please enter the allele you would like to look up:")
    allele = input()
    allele = cleanupAlleleName(allele)
    originalData.workingSet[allele].__lookupPrint__()


def lookupPrimerData(originalData):
    print("Please enter the primer name you would like to look up:")
    primer = input()
    pNum, allele = primer.split("_")
    pNum = pNum.lower()
    allele = cleanupAlleleName(allele)
    primer = pNum + "_" + allele
    if originalData.workingSet[allele].primerNameLeft == primer:
        print(f"\nLEFT PRIMER")
        originalData.workingSet[allele].__lookupPrimerPrintLeft__()
    elif originalData.workingSet[allele].primerNameRight == primer:
        print(f"\nRIGHT PRIMER")
        originalData.workingSet[allele].__lookupPrimerPrintRight__()
    else:
        print(f"Primer {primer} does not exist in dataset")


def main():
    print("Please select what you want to look up:")
    print()
    print("(1) Flanking sequence data")
    print("(2) Primer data")
    print("(3) Sanger sequence data (TBD)")
    print()
    sel = input()
    print()

    originalData = QueriesWorkingSet()
    originalData.__buildCompleteWorkingSet__("")

    if sel == "1":
        lookupFlankingSequenceData(originalData)
    elif sel == "2":
        lookupPrimerData(originalData)
    else:
        print("Goodbye.")
        exit()

    # TODO: Sanger lookup


if __name__ == "__main__":
    main()
