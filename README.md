# STRESS
Short Tandem Repeat Extraction by Seed Selection

To compile the source code STRESS.c following command can be used

$gcc STRESS.c -DMAX=2900000000 -o STRESS

The option -DMAX indicates maximum length(in decimal) of the DNA sequence for STR mining. We have used it up to 2900000000 nt. 

To run the executable file 'STRESS' following command can be used

$STRESS input.txt -l 12 -n 1 -m 6 -p 0 > output.txt

where l stands for cutoff length(here it is 12), n(here it is 1) and m(here it is 6) are minimum and maximum motif length and 0 value for p stands for without partial where as the value 1 is to mine STRs with partial motif.

File 'input.txt' contains the DNA sequence from which STRs are to be mined and the output STRs will be stored in the file 'output.txt'.  

