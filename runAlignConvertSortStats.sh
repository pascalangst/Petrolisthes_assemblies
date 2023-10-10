#!/bin/bash

#You must provide the following. Note variable DBDIR does not need a "/" at the end.
# sh runAlignConvertSortStats.sh 28 sequence_1.fastq sequence_2.fastq /work/GIF/remkv6/files genome.fa


PROC=16
R1_FQ="$1"
R2_FQ="$2"
DBDIR="$3"
GENOME="$4"

#Trim with trimgalore!
module load Trim_Galore/0.6.5-foss-2018b-Python-3.6.6

trim_galore --paired ${R1_FQ} ${R2_FQ} --cores 4


#Align with Hisat2!
#If aligning more than one set of RNAseq, build the database outside the script

module purge
module load HISAT2/2.2.1-foss-2018b

#hisat2-build ${GENOME} ${GENOME%.*}
hisat2 -p ${PROC} -x ${GENOME%.*} -1 $(basename $(echo ${R1_FQ%%.*})_val_1.fq.gz) -2 $(basename $(echo ${R2_FQ%%.*})_val_2.fq.gz) -S $(basename $(echo ${R1_FQ%.*}).sam)


#Convert sam to bam, and sort!
module purge
module load SAMtools/1.7-goolf-1.7.20
samtools view --threads ${PROC} -b -o $(basename $(echo ${R1_FQ%.*}).bam) $(basename $(echo ${R1_FQ%.*}).sam)
samtools sort -o $(basename $(echo ${R1_FQ%.*})_sorted.bam) -T $(basename $(echo ${R1_FQ%.*})_temp) --threads ${PROC} $(basename $(echo ${R1_FQ%.*}).bam)
rm $(basename $(echo ${R1_FQ%.*}).sam)