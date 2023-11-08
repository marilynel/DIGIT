# This class is not currently in use. May be useful in future versions of DIGIT.

from Query import *


class Primer:
    def __init__(self, queryObj):
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

