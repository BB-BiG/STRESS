# STRSS
Short Tandem Repeat Seed Selection

To compile the source code STRSS.c following command can be used

$gcc STRSS.c -DMAX=2900000000 -o STRSS

The option -DMAX indicates maximum length(in decimal) of the DNA sequence for STR mining. We have used it up to 2900000000 nt. 

To run the executable file 'STRSS' following command can be used

$STRSS input.txt -l 12 -n 1 -m 6 -p 0 > output.txt

where l stands for cutoff length(here it is 12), n(here it is 1) and m(here it is 6) are minimum and maximum motif length and 0 value for p stands for without partial and the value 1 is to mine STRs with partial motif.

File 'input.txt' contains the DNA sequence from which STRs are mined and stored in the file 'output.txt'.  

