Find 1.1.0
This project has been made by Andrey SOBOLEVSKY, Franck TROUILLEZ and Tristan
SMEESTERS.
The purpose of this code is to align a query sequence with all the proteins
from the database. It uses the Smith-Watermann Algorithm to get scores and to
get the alignements in an output file (result.txt by default).

Usage : ./launch [OPTIONS]
Usage
-q      name of query file (required)
-d      name of data file (uniprot_sprot.fasta)
-o      name of output file (result.txt)
-n      number of results showed (10)
-m      scoring matrix used for Smith-Waterman(blosum62)
-gpo    gap penality opening (11)
-gpe    gap penality expansion (1)
