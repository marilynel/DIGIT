# Reference genomes in .fa format cannot be used with Blast directly. This script transforms the reference genome into
# a format that can be used with Blast. If a different reference genome needs to be used, simply copy the below format
# with the preferred reference genome and comment out (#) the current lines.
#
# Format:
#   SGE_Batch -c 'makeblastdb -in <path/to/reference/genome> -out <path/to/output/directory> -dbtype nucl -title <title>
#       -parse_seqids' -q bpp -P 8 -r sge.makeblastdb_<database>

SGE_Batch -c 'makeblastdb -in DIGITfiles/Genomes/A188v1/Zm-A188-REFERENCE-KSU-1.0.fa -out DIGITfiles/Genomes/A188v1/A188v1 -dbtype nucl -title A188v1 -parse_seqids' -q bpp -P 8 -r sge.makeblastdb_A188v1
SGE_Batch -c 'makeblastdb -in DIGITfiles/Genomes/B73v5/Zm-B73-REFERENCE-NAM-5.0.fa -out DIGITfiles/Genomes/B73v5/B73v5 -dbtype nucl -title B73v5 -parse_seqids' -q bpp -P 8 -r sge.makeblastdb_B73v5
SGE_Batch -c 'makeblastdb -in DIGITfiles/Genomes/W22v2/Zm-W22-REFERENCE-NRGENE-2.0.fa -out DIGITfiles/Genomes/W22v2/W22v2 -dbtype nucl -title A188v1 -parse_seqids' -q bpp -P 8 -r sge.makeblastdb_W22v2
