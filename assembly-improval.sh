# Code for the removal of non-target contigs

## first let's get an overview of the initial assembly to afterward remove non-target contigs

## 1) PB745_01_Cinc1_PorcelainCrab_HiFi/ 

## #HiFi read location
## /scicore/home/ebertd/GROUP/Jonathon_genomeproj/PacBioData/PB745_01_Cinc1_PorcelainCrab_HiFi/r64069_20220420_183557/F1/outputs/m64069_220427_224156.hifi_reads.fasta.gz 
## /scicore/home/ebertd/GROUP/Jonathon_genomeproj/PacBioData/PB745_01_Cinc1_PorcelainCrab_HiFi/r64069_20220516_224343/F1/m64069_220521_034620.hifi_reads.fasta.gz 
## /scicore/home/ebertd/GROUP/Jonathon_genomeproj/PacBioData/PB745_01_Cinc1_PorcelainCrab_HiFi/r64069_20220516_224343/G1/m64069_220522_144333.hifi_reads.fasta.gz 

## #Assembly location
## /scicore/home/ebertd/GROUP/HIFIASM_assemblies/Cinc/HIFIASM/Cinc_HIFIASM.hifiasm.p_ctg.fa

# BUSCO analysis
busco -i Cinc_HIFIASM.hifiasm.p_ctg.fa -l /scicore/home/ebertd/angpas00/crab/hifiasm-assembly-Eric/Mani/HIFIASM/busco_downloads/lineages/arthropoda_odb10 -o Cinc_busco -m genome --offline -f --cpu 8

# prepare files for blobtools analysis

# blast all contigs against ncbi and uniprot
blastn -db /home/pascal/bioinformatics/ncbi/nt \
       -query Cinc_HIFIASM.hifiasm.p_ctg.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -out blast.out \
       -num_threads 20
       
diamond blastx \
      --query Cinc_HIFIASM.hifiasm.p_ctg.fa \
      --db uniprot/reference_proteomes.dmnd \
      --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
      --sensitive \
      --max-target-seqs 1 \
      --evalue 1e-25 \
      --threads 42 \
      -F 15 -b4 -c1 \
      > diamond-scicore.all.out    

# mapping reads to assembly
minimap2 -ax map-hifi -t 54 Cinc_HIFIASM.hifiasm.p_ctg.fa m64069_220427_224156.hifi_reads.fasta.gz m64069_220521_034620.hifi_reads.fasta.gz m64069_220522_144333.hifi_reads.fasta.gz > Cinc_HIFIASM.hifiasm.p_ctg.sam

samtools view -Sb -@ 10 Cinc_HIFIASM.hifiasm.p_ctg.sam | samtools sort -@ 10 -T xyz -o Cinc_HIFIASM.hifiasm.p_ctg.bam
rm Cinc_HIFIASM.hifiasm.p_ctg.sam
samtools index -b Cinc_HIFIASM.hifiasm.p_ctg.bam

# mapping RNA-seq reads to assembly
find /scicore/home/ebertd/GROUP/Jonathon_genomeproj/RNAseqData/Cinctipes_Data/ -name *fastq.gz > RNA-seq_files.list
cat RNA-seq_files.list | tr '\n' ','

STAR --runThreadN 20 --runMode genomeGenerate --genomeDir . --genomeFastaFiles Cinc_HIFIASM.hifiasm.p_ctg.fa
STAR --runThreadN 20 --genomeDir . --readFilesIn RNA-seq_files.list --readFilesCommand zcat

samtools view -Sb -@ 10 Aligned.out.sam | samtools sort -@ 10 -T xyz -o Aligned.out.bam
rm Aligned.out.sam
samtools index -b Aligned.out.bam

# mapping 10x reads to assembly
bwa-mem2 index Cinc_HIFIASM.hifiasm.p_ctg.fa
bwa-mem2 mem -t 55 -M Cinc_HIFIASM.hifiasm.p_ctg.fa /scicore/home/ebertd/angpas00/crab/Crab2_S1_L006_R1_001.fastq.gz /scicore/home/ebertd/angpas00/crab/Crab2_S1_L006_R2_001.fastq.gz > Cinc_HIFIASM.hifiasm.p_ctg.sam

samtools view -Sb -@ 10 Cinc_HIFIASM.hifiasm.p_ctg.sam | samtools sort -@ 10 -T xyz -o Cinc_HIFIASM.hifiasm.p_ctg.bam
rm Cinc_HIFIASM.hifiasm.p_ctg.sam
samtools index -b Cinc_HIFIASM.hifiasm.p_ctg.bam


# run blobtools
blobtools nodesdb --nodes ~/bioinformatics/taxdump/nodes.dmp --names ~/bioinformatics/taxdump/names.dmp

blobtools create -i Cinc_HIFIASM.hifiasm.p_ctg.fa -b Cinc_HIFIASM.hifiasm.p_ctg.bam -b Cinc_HIFIASM.hifiasm.p_ctg_10x.bam \
 -t blast.out -t diamond-scicore.all.out -o Cinc_HIFIASM.hifiasm.p_ctg.blobtools

blobtools view -i Cinc_HIFIASM.hifiasm.p_ctg.blobtools.blobDB.json
blobtools plot -i Cinc_HIFIASM.hifiasm.p_ctg.blobtools.blobDB.json

# blobtools in covplot mode (using RNA-seq as second source)   
blobtools create -i Cinc_HIFIASM.hifiasm.p_ctg.fa -b Cinc_HIFIASM.hifiasm.p_ctg.bam -b Aligned.out.bam \
 -t blast.out -o Cinc_HIFIASM.hifiasm.p_ctg.blobtools-covplot

blobtools view -i Cinc_HIFIASM.hifiasm.p_ctg.blobtools-covplot.blobDB.json
blobtools covplot -i Cinc_HIFIASM.hifiasm.p_ctg.blobtools-covplot.blobDB.json --cov Cinc_HIFIASM.hifiasm.p_ctg.blobtools-covplot.Aligned.out.bam.cov

# the above steps need to be repeated for each iteration of cleaning
# Assembly version history:
# Cinc_v1.0.1: cov: 5, length 2kb, GC < 45, GC > 34, only arthopoda - no-hit - undef - unresolved
# Cinc_v1.1: length 2kb, GC < 50, GC > 30
# Cinc_v1.2: length 2kb, GC < 50, GC > 30, no Bacteria
# Cinc_v1.3: length 2kb, GC < 50, GC > 30, no Bacteria, Taxa with mean coverage in 10x and HiFi >= 1, Exclude Priapulida (mean 10x=131.9, std=0; mean HIFI =4.2, std=0)
# Cinc_v1.4: (diamond info available from here), length 2kb, GC < 50, GC > 30, Taxa with mean coverage in 10x and HiFi >= 1, Contigs with mean coverage in 10x >= 1, Contigs with mean coverage in HiFi >= 3.5, Exclude Bacteria, Exclude Priapulida, Exclude Brachiopoda, Exclude Ascomycota, Exclude unresolved, Exclude Annelida
# Cinc_v1.5: additionally set an upper threshold for coverage (3x mean): Contigs with mean coverage in 10x: 1 <= x <= 285.9; Contigs with mean coverage in HiFi: 3.5 <= x <= 34.8. Additionally exclude Streptophyta and circular contigs
# Cinc_v1.6-1.8: further trials
# Cinc_v1.9: going back to v1.5

# filter by length GC and coverage. Then filter by taxon. Then filter circular contigs
awk '{if ($2 > 1999 && $3 > 0.3 && $3 < 0.5 && $5 >= 3.5 && $6 >= 1 && $5 <= 34.8 && $6 <= 285.9) print $0}' Cinc_HIFIASM.hifiasm.p_ctg.blobtools.blobDB.bestsum.table.txt  | grep -v -f Cinc_v1.5.list | cut -f1 | grep -v "c" > Cinc_v1.5.names
seqtk subseq Cinc_HIFIASM.hifiasm.p_ctg.fa Cinc_v1.5.names > Cinc_v1.5.fa





## 2) PB745_02_Mani_PorcelainCrab_HiFi/ 
        
## #HiFi read location
## /scicore/home/ebertd/GROUP/Jonathon_genomeproj/PacBioData/PB745_02_Mani_PorcelainCrab_HiFi/r64069_20220420_183557/G1/outputs/m64069_220429_093835.hifi_reads.fasta.gz 
## /scicore/home/ebertd/GROUP/Jonathon_genomeproj/PacBioData/PB745_02_Mani_PorcelainCrab_HiFi/r64069_20220525_225441/F1/m64069_220530_114955.hifi_reads.fasta.gz  
## /scicore/home/ebertd/GROUP/Jonathon_genomeproj/PacBioData/PB745_02_Mani_PorcelainCrab_HiFi/r64069_20220525_225441/G1/m64069_220531_224633.hifi_reads.fasta.gz
        
## #Assembly location
## /scicore/home/ebertd/GROUP/HIFIASM_assemblies/Mani/HIFIASM/Mani_HIFIASM.hifiasm.p_ctg.gfa
# translate to fasta
awk '/^S/{print ">"$2;print $3}' /scicore/home/ebertd/GROUP/HIFIASM_assemblies/Mani/HIFIASM/Mani_HIFIASM.hifiasm.p_ctg.gfa > Mani_HIFIASM.hifiasm.p_ctg.fa

# BUSCO analysis
busco -i Mani_HIFIASM.hifiasm.p_ctg.fa -l /scicore/home/ebertd/angpas00/crab/hifiasm-assembly-Eric/Mani/HIFIASM/busco_downloads/lineages/arthropoda_odb10 -o Mani_busco -m genome --offline -f --cpu 8

# prepare files for blobtools analysis

# blast all contigs against ncbi and uniprot
blastn -db /home/pascal/bioinformatics/ncbi/nt \
       -query Mani_HIFIASM.hifiasm.p_ctg.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -out blast.out \
       -num_threads 20
       
diamond blastx \
      --query Mani_HIFIASM.hifiasm.p_ctg.fa \
      --db uniprot/reference_proteomes.dmnd \
      --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
      --sensitive \
      --max-target-seqs 1 \
      --evalue 1e-25 \
      --threads 42 \
      -F 15 -b4 -c1 \
      > diamond-scicore.all.out  

# mapping reads to assembly
minimap2 -ax map-hifi -t 54 Mani_HIFIASM.hifiasm.p_ctg.fa m64069_220429_093835.hifi_reads.fasta.gz m64069_220530_114955.hifi_reads.fasta.gz m64069_220531_224633.hifi_reads.fasta.gz  > Mani_HIFIASM.hifiasm.p_ctg.sam

samtools view -Sb -@ 10 Mani_HIFIASM.hifiasm.p_ctg.sam | samtools sort -@ 10 -T xyz -o Mani_HIFIASM.hifiasm.p_ctg.bam
rm Mani_HIFIASM.hifiasm.p_ctg.sam
samtools index -b Mani_HIFIASM.hifiasm.p_ctg.bam

# mapping RNA-seq reads to assembly
find /scicore/home/ebertd/GROUP/Jonathon_genomeproj/RNAseqData/Manimaculis_data/ -name *fastq.gz > RNA-seq_files.list
cat RNA-seq_files.list | tr '\n' ','

STAR --runThreadN 20 --runMode genomeGenerate --genomeDir . --genomeFastaFiles Mani_HIFIASM.hifiasm.p_ctg.fa
STAR --runThreadN 20 --genomeDir . --readFilesIn RNA-seq_files.list --readFilesCommand zcat

samtools view -Sb -@ 10 Aligned.out.sam | samtools sort -@ 10 -T xyz -o Aligned.out.bam
rm Aligned.out.sam
samtools index -b Aligned.out.bam


# run blobtools
blobtools create -i Mani_HIFIASM.hifiasm.p_ctg.fa -b Mani_HIFIASM.hifiasm.p_ctg.bam \
 -t blast.out -t diamond-scicore.all.out -o Mani_HIFIASM.hifiasm.p_ctg.blobtools

blobtools view -i Mani_HIFIASM.hifiasm.p_ctg.blobtools.blobDB.json
blobtools plot -i Mani_HIFIASM.hifiasm.p_ctg.blobtools.blobDB.json      
        
# blobtools in covplot mode (using RNA-seq as second source)      
blobtools create -i Mani_HIFIASM.hifiasm.p_ctg.fa -b Mani_HIFIASM.hifiasm.p_ctg.bam -b Aligned.out.bam \
 -t blast.out -o Mani_HIFIASM.hifiasm.p_ctg.blobtools-covplot

blobtools view -i Mani_HIFIASM.hifiasm.p_ctg.blobtools-covplot.blobDB.json
blobtools covplot -i Mani_HIFIASM.hifiasm.p_ctg.blobtools-covplot.blobDB.json --cov Mani_HIFIASM.hifiasm.p_ctg.blobtools-covplot.Aligned.out.bam.cov

# contaminant contig filtration was done in a stepwise manner as described above. Additional filters were to exclude no-hits in the taxon blasts which also had a GC below 33 and one relatively short contig that did not have all four nucleotides.

# filter by length GC and coverage. Then filter by taxon. Then filter circular contigs
awk '{if ($2 > 1999 && $3 > 0.3 && $3 < 0.5 && $5 >= 4.37 && $5 <= 30) print $0}'  Mani_HIFIASM.hifiasm.p_ctg.blobtools.blobDB.bestsum.table.txt | awk '{if ($6 !="no-hit" || $3 > 0.33) print $0}' | grep -v -f Mani_v1.6.list | cut -f1 | grep -v "c" > Mani_v1.6.names
seqtk subseq Mani_HIFIASM.hifiasm.p_ctg.fa Mani_v1.6.names > Mani_v1.6.fa
