#!/bin/python
import os
import sys

#Usage: runPathoPairedRefseq.py <text_file_with_reads_filenames>
#Note: The text file needs to have pairs of reads, an output foldername, and a samplename for them all on a single line, separated by tabs. Full paths.

pairedReads = sys.argv[1]

with open(pairedReads,"r") as f:
    for line in f:
        R1,R2,folder,sample = line.strip().split()
        os.system("qsub -P pathoscope -cwd -V -pe omp 32 -l h_rt=48:00:00 -e " + folder + "/errorLog" + sample + ".txt -o " + folder + "/outLog" + sample + ".txt -N mapping" + sample + "  /restricted/projectnb/pathoscope/code/wrapperForPathoScope/scratchPairedRefseq.sh " + R1 + " " + R2 + " " + folder + " " + sample)

