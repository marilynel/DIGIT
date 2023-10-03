from Utils import *


def main():
    # 1. Load entire query working set
    #       a.B73 only or all datasets?
    #       b. allQueries?
    originalData = buildFullWorkingSet(False)

    print("Please enter the allele you would like to look up:")
    allele = input()
    allele = allele.upper()
    if allele[-1] == "A" or allele[-1] == "B":
        letter = allele[-1]
        letter.lower()
        allele = allele[:-1] + letter
    for q in originalData:
        if originalData[q].query == allele:
            originalData[q].__lookupPrint__()

    # 2. Load sanger data
    #       a. B73 only or all datasets?
    # 3. prompt user for input
    # 4. Print


if __name__ == "__main__":
    main()
