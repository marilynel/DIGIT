import sys

'''

0.  get everything we need for a specific genome in one place

    {
        query_upper = [chrX, [c1,c2]],
        query_lower = [chrN, [c3,c4]],
        etc.
    }

1.  open that giant ass file
2.  save parts of file to a struct maybe?
    {
        chr1:sequence,
        chr2:sequence,
        etc.
    }

3.  use the coordinates to index in on that struct
4.  return sequence

'''

coordinates = {
    "test1" : ["chr1", [500, 600]],
    "test2" : ["chr1", [550, 650]]
}



def main():
    filename = "Zm-A188-REFERENCE-KSU-1.0.fa.gz"
    genomeStruct = {}
    with open(filename, "r") as zfile:
        key = ""
        for line in zfile:
            if line[0] == ">":
                key = line[1:]
            else:
                genomeStruct[key] = line

    # will need to accept the coordinates as in above

    sequenceData = {}
    for c in coordinates:
        c1, c2 = coordinates[1][0], coordinates[1][1]     # these should be ints
        chr = coordinates[0]

        seq = genomeStruct[chr][c1:c2]

        sequenceData[c] = seq


if __name__ == "__main__":
    main()