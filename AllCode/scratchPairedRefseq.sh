#!/bin/bash

#This script should run Pathoscope2 on any paired-end reads given to it. It will make an output directory for it.
#If this script gets moved or copied, it will need adjustments for the new filepaths

#Usage: python scratchPairedRefseq.sh <R1.fastq> <R2.fastq> <outfolder_full_path> <SampleName>

mkdir $3

python2 /restricted/projectnb/pathoscope/code/PathoScope/pathoscope/pathoscope2.py  MAP -1 $1 -2 $2 -targetIndexPrefixes bacteria_0,bacteria_1,bacteria_2,bacteria_3,bacteria_4,bacteria_5 -outDir $TMPDIR -outAlign $4.sam -numThreads 32 -targetAlignParams "--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 10 --score-min L,100,1.28"


python2 /restricted/projectnb/pathoscope/code/PathoScope/pathoscope/pathoscope2.py ID --noUpdatedAlignFile -outDir $3 -alignFile $TMPDIR/$4.sam

rm -f $TMPDIR/$4.sam
