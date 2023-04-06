import os
import subprocess
import time


# from DIGITfiles.Query import Query
# from DIGITfiles.Sequences import Sequences

def okGo(dirname, fileSubstrings):
    ok = False
    listDirs = os.listdir(dirname)
    for fileSubstring in fileSubstrings:
        for dir in listDirs:
            if dir.find(fileSubstring) != -1:
                ok = True
                break
        else:
            print(
                "You do not have the correct files available to complete this task. Check to make sure you have completed the previous steps in order to proceed. See README for more details.")
            exit()
    return ok


def runSequencesFromBlast():
    listDirs = os.listdir("DIGITfiles/BlastOutput")

    print(
        f"\nSelect a directory to work with. Make sure the folder you need is in the DIGITfiles/BlastOutput/ folder.\n")

    for i in range(0, len(listDirs)):
        print(f"({i + 1}) {listDirs[i]}")

    print()
    sel = input()

    try:
        if okGo("DIGITfiles/BlastOutput/" + listDirs[int(sel) - 1], ["A188", "B73", "W22"]):
            print(f"\nFinding genomic sequences for alleles in {listDirs[int(sel) - 1]}.")
            now = time.time()
            subprocess.run(
                [
                    "SGE_Batch",
                    "-c",
                    f"python3 DIGITfiles/SequencesFromBlast.py {listDirs[int(sel) - 1]}",
                    "-q",
                    "bpp",
                    "-P",
                    "8",
                    "-r",
                    f"sge.seqFromBlast_{listDirs[int(sel) - 1]}_{now}"
                ]
            )

            print(
                "Return later to continue to continue the primer making process with building the predicted insertion sequences. SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running time may vary wildly.")

    except:
        print(f"\nPlease don't be weird, just select a real directory. Exiting.")


def runGetPrimersAndVerify():
    listDirs = os.listdir("DIGIToutput/")

    print(
        f"\nSelect a directory to work with. Make sure you have completed Part 1 for the dataset you are looking for.\n")

    for i in range(0, len(listDirs)):
        print(f"({i + 1}) {listDirs[i]}")
    print()

    sel = input()

    try:
        if okGo("DIGIToutput/" + listDirs[int(sel) - 1], ["WorkingSet", "Primer3Output"]):
            # print(listFiles[int(sel) - 1])
            print(
                f"\nBuilding DsGFP insertion sequences and creating Primer3 input for alleles in {listDirs[int(sel) - 1]}.")
            now = time.time()
            subprocess.run(
                [
                    "SGE_Batch",
                    "-c",
                    f"python3 DIGITfiles/GetPrimers.py {listDirs[int(sel) - 1]}",
                    "-q",
                    "bpp",
                    "-P",
                    "8",
                    "-r",
                    f"sge.getPrimers_{listDirs[int(sel) - 1]}_{now}"
                ]
            )

            print(
                "Return later to continue to continue the primer making process with building the predicted insertion sequences. SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running time may vary wildly.")

    except:
        print(f"\nPlease don't be weird, just select a real directory. Exiting.")


def main():
    print(
        f"Welcome to DIGIT, the premiere tool for predicting DsGFP insertion" +
        f" sequences in maize and building primers for those sequences. " +
        f"Please select from the following menu options:\n\n" +
        f"(1) Parse Blast results and find primers\n" +
        f"(2) Verify primers from part 1\n"
    )
    sel = input()
    if sel == "1":
        runSequencesFromBlast()
    elif sel == "2":
        runGetPrimersAndVerify()
    else:
        print("Bye")


if __name__ == "__main__":
    main()