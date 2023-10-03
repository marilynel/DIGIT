# from Query import Query

class PrimerIDs:
    def __init__(self):
        self.listOfPrimers = []
        self.index = 0
        with open("DIGITfiles/compListPreviousPrimers", "r") as infile:
            for line in infile:
                self.listOfPrimers.append(line.strip())

        self.numPrimers = len(self.listOfPrimers) + 1

    def __iter__(self):
        return self

    def __next__(self):
        if self.index < len(self.listOfPrimers):
            result = self.listOfPrimers[self.index]
            self.index += 1
            return result
        else:
            raise StopIteration

    def __rewriteRecordFile__(self):
        with open("DIGITfiles/compListPreviousPrimers", "w") as outfile:
            for primer in self.listOfPrimers:
                outfile.write(primer + "\n")