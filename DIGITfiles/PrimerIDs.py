from Query import Query


class PrimerIDs:
    def __init__(self):
        self.listOfPrimers = []

        with open("DIGITfiles/compListPreviousPrimers", "r") as infile:
            for line in infile:
                self.listOfPrimers.append(line.strip())

        self.numPrimers = len(self.listOfPrimers) + 1

    def __setPrimerIDs__(self, queriesWorkingSet):
        for query in queriesWorkingSet:
            leftNameExists, rightNameExists = False, False
            for primerName in self.listOfPrimers:
                pts = primerName.split("_")
                noNum = pts[0][-1] + "_" + pts[1]
                # is left primer name in there? If so, rename left primer
                if noNum == queriesWorkingSet[query].primerNameLeft:
                    queriesWorkingSet[query].primerNameLeft = primerName
                    leftNameExists = True
                # is right primer in there? if so, rename right primer
                if noNum == queriesWorkingSet[query].primerNameRight:
                    queriesWorkingSet[query].primerNameRight = primerName
                    rightNameExists = True

            if not leftNameExists:
                queriesWorkingSet[query].primerNameLeft = (
                        "mr" + str(self.numPrimers) + queriesWorkingSet[query].primerNameLeft
                )
                self.listOfPrimers.append(queriesWorkingSet[query].primerNameLeft)
                self.numPrimers += 1
            if not rightNameExists:
                queriesWorkingSet[query].primerNameRight = (
                        "mr" + str(self.numPrimers) + queriesWorkingSet[query].primerNameRight
                )
                self.listOfPrimers.append(queriesWorkingSet[query].primerNameRight)
                self.numPrimers += 1
        self.__rewriteRecordFile__()

    def __rewriteRecordFile__(self):
        with open("DIGITfiles/compListPreviousPrimers", "w") as outfile:
            for primer in self.listOfPrimers:
                outfile.write(primer + "\n")