#!/bin/bash

#SBATCH --job-name=HIFIASM_Cinq					#Job name
#SBATCH --cpus-per-task=8	                  	#Number of cores reserved
#SBATCH --mem-per-cpu=32G              			#Memory reserved per core.
												#Total memory reserved: 16GB
#SBATCH --time=168:00:00	        			#Maximum time the job will run
#SBATCH --qos=1week          					#The job queue (time based)

#This is the stdout file
#SBATCH --output=/scicore/home/ebertd/dexter0000/crabs/logs/out_HIFIASM_Cinq

#This is the stderr file
#SBATCH --error=/scicore/home/ebertd/dexter0000/crabs/logs/err_HIFIASM_Cinq

#Things remember when running on sciCore:
#This job runs from the current working directory
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

################################################################################
#Load modules and set variables
################################################################################

#Load conda environment
conda activate eukaryotic_genome_assembly

#load hifiasm
module load hifiasm

#Fasta of reads to be assembled
READS_FASTA_1="/scicore/home/ebertd/GROUP/Jonathon_genomeproj/PacBioData/PB745_01_Cinc1_PorcelainCrab_HiFi/r64069_20220516_224343/F1/m64069_220521_034620.hifi_reads.fasta.gz"
READS_FASTA_2="/scicore/home/ebertd/GROUP/Jonathon_genomeproj/PacBioData/PB745_01_Cinc1_PorcelainCrab_HiFi/r64069_20220516_224343/G1/m64069_220522_144333.hifi_reads.fasta.gz"
READS_FASTA_3="/scicore/home/ebertd/GROUP/Jonathon_genomeproj/PacBioData/PB745_01_Cinc1_PorcelainCrab_HiFi/r64069_20220420_183557/F1/outputs/m64069_220427_224156.hifi_reads.fasta.gz"

#Prefix for assembly output files
ASSEMBLY="Cinc_HIFIASM"

################################################################################
Assemble genome
################################################################################

#Navigate to project directory
cd /scicore/home/ebertd/dexter0000/crabs

#Create output directory
mkdir Cinc/HIFIASM

#This runs amazingly fast (~1 day)
hifiasm -f0 --primary -t 8 -o Cinc/HIFIASM/"$ASSEMBLY".hifiasm "$READS_FASTA_1" "$READS_FASTA_2" "$READS_FASTA_3"