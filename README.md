MakeGencCAGETSS
===============

Makes TSS clusters from Gencode TSS and CAGE peaks

This script takes as input:
- a stranded gencode TSS file in gff format (already filtered for TSS quality and with no header)
- a stranded cage peak file in bed format (already filtered for quality and with no header)

This script provides as output in the directory where it is executed the following files (indexed by the inputs' names):
- a "genccage" (i.e. gencode+cage) tss file in gff 2 format with information about the initial elements it is made of and in particular the coordinates of the initial cage peak from which it is made
- some statistics in standard output and some log information in standard error 

It is important to note that this script requires **bedtools**


