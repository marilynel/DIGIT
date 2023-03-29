import os
import subprocess


# from DIGITfiles.Query import Query
# from DIGITfiles.Sequences import Sequences

def runSequencesFromBlast():
    listDirs = os.listdir("DIGITfiles/BlastOutput")
    for i in range(0, len(listDirs)):
        print(f"({i + 1}) {listDirs[i]}")

    print(
        "Select a directory to work with. Make sure the folder you need is in the DIGITfiles/BlastOutput/ folder.")
    sel = input()

    try:
        print(f"Finding genomic sequences for alleles in {listDirs[int(sel) - 1]}.")
        subprocess.run(
            [
                "python3",
                "DIGITfiles/SequencesFromBlast.py",
                listDirs[int(sel) - 1]
            ]
        )

        print(
            "Return later to continue to continue the primer making process with building the predicted insertion sequences. SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running time may vary wildly.")

    except:
        print("Please don't be weird, just select a real directory. Exiting.")


def runBuildPredictedSequences():
    print(
        "Select the dataset to work with. Make sure you completed the first part, otherwise this section of the program will not work.")

    listFiles = os.listdir("DIGITfiles/")
    for i in range(0, len(listFiles)):
        if listFiles[i].find("json") != -1 and listFiles[i].find("AllBlastData") != -1:
            print(f"({i + 1}) {listFiles[i]}")

    sel = input()

    try:
        print(
            f"Building DsGFP insertion sequences and creating Primer3 input for alleles in {listFiles[int(sel) - 1]}.")
        subprocess.run(
            [
                "python3",
                "DIGITfiles/BuildPredictedSequences.py",
                listFiles[int(sel) - 1]
            ]
        )

        print(
            "Return later to continue to continue the primer making process with building the predicted insertion sequences. SGE may take a while to run. To check the status of your jobs, enter 'qstat' in the command line. Running time may vary wildly.")

    except:
        print("Please don't be weird, just select a real directory. Exiting.")


def main():
    print(
        f"Welcome to DIGIT, the premiere tool for predicting DsGFP insertion" +
        f" sequences in maize and building primers for those sequences. " +
        f"Enter 1 if you have new sequences you would like to find primers " +
        f"for. Enter 2 if your earlier jobs are complete and you are ready to " +
        f"start building your predicted DsGFP sequences."
    )
    sel = input()
    if sel == "1":
        runSequencesFromBlast()
    elif sel == "2":
        runBuildPredictedSequences()
    else:
        print("Bye")


if __name__ == "__main__":
    main()