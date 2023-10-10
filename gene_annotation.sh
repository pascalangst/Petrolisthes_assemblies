## Gene prediction and functional annotation of contaminant reduced crab assemblies

# creat database for rna-seq mapping
hisat2-build Mani_v1.7.fa.masked Mani_v1.7.fa

# map rna-seq reads to reference by executing the following for all read-sets
#sh runAlignConvertSortStats.sh $r1 $r2 . Mani_v1.7.fa.masked

# gene prediction using braker2 supplying the contaminant reduced genome and the mapped rna-reads
braker.pl --species=Petrolisthes_manimaculis --genome=Mani_v1.7.fa.masked --bam=EA1_10_CTTGTA_L006_R1_001.fastq_sorted.bam,EA1_10_CTTGTA_L006_R1_002.fastq_sorted.bam,EA1_10_CTTGTA_L006_R1_003.fastq_sorted.bam,EA1_17_TAGCTT_L006_R1_001.fastq_sorted.bam,EA1_17_TAGCTT_L006_R1_002.fastq_sorted.bam,EA1_17_TAGCTT_L006_R1_003.fastq_sorted.bam,EA1_19_GGCTAC_L006_R1_001.fastq_sorted.bam,EA1_19_GGCTAC_L006_R1_002.fastq_sorted.bam,EA1_19_GGCTAC_L006_R1_003.fastq_sorted.bam,EA1_22_GCCAAT_L006_R1_001.fastq_sorted.bam,EA1_22_GCCAAT_L006_R1_002.fastq_sorted.bam,EA1_22_GCCAAT_L006_R1_003.fastq_sorted.bam,EA1_46_AGTCAA_L006_R1_001.fastq_sorted.bam,EA1_46_AGTCAA_L006_R1_002.fastq_sorted.bam,EA1_46_AGTCAA_L006_R1_003.fastq_sorted.bam,EA1_48_GTCCGC_L006_R1_001.fastq_sorted.bam,EA1_48_GTCCGC_L006_R1_002.fastq_sorted.bam,EA1_48_GTCCGC_L006_R1_003.fastq_sorted.bam,EA1_4_TTAGGC_L006_R1_001.fastq_sorted.bam,EA1_4_TTAGGC_L006_R1_002.fastq_sorted.bam,EA1_4_TTAGGC_L006_R1_003.fastq_sorted.bam,EA1_4_TTAGGC_L006_R1_004.fastq_sorted.bam,EA1_50_CGTACG_L006_R1_001.fastq_sorted.bam,EA1_50_CGTACG_L006_R1_002.fastq_sorted.bam,EA1_50_CGTACG_L006_R1_003.fastq_sorted.bam,EA1_51_ATGTCA_L006_R1_001.fastq_sorted.bam,EA1_51_ATGTCA_L006_R1_002.fastq_sorted.bam,EA1_51_ATGTCA_L006_R1_003.fastq_sorted.bam,EA1_53_GTTTCG_L006_R1_001.fastq_sorted.bam,EA1_53_GTTTCG_L006_R1_002.fastq_sorted.bam,EA1_53_GTTTCG_L006_R1_003.fastq_sorted.bam,EA1_53_GTTTCG_L006_R1_004.fastq_sorted.bam,EA1_8_TGACCA_L006_R1_001.fastq_sorted.bam,EA1_8_TGACCA_L006_R1_002.fastq_sorted.bam,EA1_8_TGACCA_L006_R1_003.fastq_sorted.bam,EA2_11_TTAGGC_L007_R1_001.fastq_sorted.bam,EA2_11_TTAGGC_L007_R1_002.fastq_sorted.bam,EA2_11_TTAGGC_L007_R1_003.fastq_sorted.bam,EA2_11_TTAGGC_L007_R1_004.fastq_sorted.bam,EA2_11_TTAGGC_L007_R1_005.fastq_sorted.bam,EA2_13_ACTTGA_L007_R1_001.fastq_sorted.bam,EA2_13_ACTTGA_L007_R1_002.fastq_sorted.bam,EA2_13_ACTTGA_L007_R1_003.fastq_sorted.bam,EA2_13_ACTTGA_L007_R1_004.fastq_sorted.bam,EA2_15_GATCAG_L007_R1_001.fastq_sorted.bam,EA2_15_GATCAG_L007_R1_002.fastq_sorted.bam,EA2_15_GATCAG_L007_R1_003.fastq_sorted.bam,EA2_26_CGATGT_L007_R1_001.fastq_sorted.bam,EA2_26_CGATGT_L007_R1_002.fastq_sorted.bam,EA2_26_CGATGT_L007_R1_003.fastq_sorted.bam,EA2_2_CAGATC_L007_R1_001.fastq_sorted.bam,EA2_2_CAGATC_L007_R1_002.fastq_sorted.bam,EA2_2_CAGATC_L007_R1_003.fastq_sorted.bam,EA2_55_AGTCAA_L007_R1_001.fastq_sorted.bam,EA2_55_AGTCAA_L007_R1_002.fastq_sorted.bam,EA2_55_AGTCAA_L007_R1_003.fastq_sorted.bam,EA2_57_ATGTCA_L007_R1_001.fastq_sorted.bam,EA2_57_ATGTCA_L007_R1_002.fastq_sorted.bam,EA2_57_ATGTCA_L007_R1_003.fastq_sorted.bam,EA2_59_ATTCCT_L007_R1_001.fastq_sorted.bam,EA2_59_ATTCCT_L007_R1_002.fastq_sorted.bam,EA2_59_ATTCCT_L007_R1_003.fastq_sorted.bam,EA2_62_GAGTGG_L007_R1_001.fastq_sorted.bam,EA2_62_GAGTGG_L007_R1_002.fastq_sorted.bam,EA2_62_GAGTGG_L007_R1_003.fastq_sorted.bam,EA2_62_GAGTGG_L007_R1_004.fastq_sorted.bam,EA2_64_GTGAAA_L007_R1_001.fastq_sorted.bam,EA2_64_GTGAAA_L007_R1_002.fastq_sorted.bam,EA2_64_GTGAAA_L007_R1_003.fastq_sorted.bam,EA2_6_GGCTAC_L007_R1_001.fastq_sorted.bam,EA2_6_GGCTAC_L007_R1_002.fastq_sorted.bam,EA2_6_GGCTAC_L007_R1_003.fastq_sorted.bam --softmasking --cores 8 --workingdir=braker.softmask.hisat2

# convert annotation file from gtf to gff3 using agat
conda activate agat
grep -v "transcript " braker.softmask.hisat2_v1.7.gtf | grep -v "GeneMark.hmm" > braker_v1.7.AUG.tmp.gtf

agat_convert_sp_gxf2gxf.pl -g braker_v1.7.AUG.tmp.gtf -o braker_v1.7.AUG.gff3

# remove NCBI identifyied contaminants
awk '(NR==FNR) { toRemove[$1]; next }
     /^>/ { p=1; for(h in toRemove) if ( h ~ $0) p=0 }
    p' remove_from_v1.7.list  Mani_v1.7.fa > Mani_v1.8.fa
grep -v -f <(sed -e 's/>//' remove_from_v1.7.list)  braker_v1.7.AUG.gff3 >  braker_v1.8.AUG.gff3

# remove NCBI identifyied adaptors
grep -o '^>.*' Mani_v1.8.fa | sed 's/>//' > adaptor_regions.txt # do modifications as below

#ptg001234l:33-
#ptg002302l:571-
#ptg002345l:1-189685
#ptg002456l:711-
#ptg002472l:1-51573
#ptg003233l:1-37228
#ptg003378l:1-63925
#ptg006725l:1-105413
#ptg008210l:1-93235
#ptg009014l:1-66921
#ptg010163l:859-
#ptg011889l:1-43910
#ptg012161l:1-99714
#ptg012202l:246-
#ptg013009l:1-29922
#ptg013161l:533-
#ptg013211l:1-63507
#ptg013349l:1-20637
#ptg014587l:1066-
#ptg014895l:353-
#ptg001236l:1-288308
#ptg001236l:288353-325773

mkdir tmp
cd tmp
while IFS= read -r region; do samtools faidx ../Mani_v1.8.fa "$region" > "$region.fa"; done < ../adaptor_regions.txt
cat *.fa > Mani_v1.8.1.fasta
sed -i 's/ptg001236l:288353/ptg100000l:288353/' Mani_v1.8.1.fasta
sed -i -E 's/^(>ptg[0-9l]*):.*$/\1/' Mani_v1.8.1.fasta
mv Mani_v1.8.1.fasta ../Mani_v1.9.fasta
cd ../
rm -r tmp

awk -F'\t' '!/^#/ {OFS="\t"; if($1 == "ptg001234l" && $4 > 33) $4 -= 32;    if($1 == "ptg001234l" && $5 > 33)  $5 -= 32;
                             if($1 == "ptg002302l" && $4 > 571) $4 -= 570;  if($1 == "ptg002302l" && $5 > 571) $5 -= 570;
                             if($1 == "ptg002456l" && $4 > 711) $4 -= 710;  if($1 == "ptg002456l" && $5 > 711) $5 -= 710;
                             if($1 == "ptg010163l" && $4 > 859) $4 -= 858;  if($1 == "ptg010163l" && $5 > 859) $5 -= 858;
                             if($1 == "ptg012202l" && $4 > 246) $4 -= 245;  if($1 == "ptg012202l" && $5 > 246) $5 -= 245;
                             if($1 == "ptg013161l" && $4 > 533) $4 -= 532;  if($1 == "ptg013161l" && $5 > 533) $5 -= 532;
                             if($1 == "ptg014587l" && $4 > 1066) $4 -= 1065; if($1 == "ptg014587l" && $5 > 1066) $5 -= 1065;
                             if($1 == "ptg014895l" && $4 > 353) $4 -= 352;  if($1 == "ptg014895l" && $5 > 353) $5 -= 352;
                             if($1 == "ptg001236l" && $5 > 288352) $1 = "ptg100000l"; if($1 == "ptg100000l" && $4 > 288352) $4 -= 288352;  if($1 == "ptg100000l" && $5 > 288352) $5 -= 288352;
                             ; print}' braker_v1.8.AUG.gff3 > braker_v1.9.AUG.gff3
                             
# create multi line reference and extract protein sequences
fold -w 60 Mani_v1.9.fasta > Mani_v1.9.fold.fa
agat_sp_extract_sequences.pl --gff braker_v1.9.AUG.gff3 -f Mani_v1.9.fold.fa -p -o braker_v1.9.AUG.proteins.fasta

# functional annotation using InterProScan5 (create single line fasta and remove terminal stop codon first)
conda activate interproscan
sed -E 's/^(>.*) gene=(.*) seq.*/\1 \2/' braker_v1.9.AUG.proteins.fasta | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | sed -e 's/\*$//' > Mani_v1.9.AUG.proteins.fasta
python ~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/iprscan-local.py -c 6 -i Mani_v1.9.AUG.proteins.fasta -m local -o funannotate.AUG.Mani_v1.9/annotate_misc/iprscan.xml
conda deactivate

# functional annotation using PFAM, UniProt, EggNog, MEROPS, dbCAN, BUSCO, Phobius, SignalP, and InterProScan5 in funannotate
conda activate mamba-funannotate
export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate outgroups --show_buscos
funannotate setup -b arthropoda
funannotate annotate --gff braker_v1.9.AUG.gff3 --fasta Mani_v1.9.fasta -s "Mani" -o funannotate.AUG.Mani_v1.9 --cpus 8 --busco_db arthropoda --iprscan funannotate.AUG.Mani_v1.9/annotate_misc/iprscan.xml
#-------------------------------------------------------
#[Oct 03 03:52 PM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.12
#[Oct 03 03:52 PM]: Running 1.8.11
#[Oct 03 03:52 PM]: No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'
#[Oct 03 03:52 PM]: Found existing output directory funannotate.AUG.Mani_v1.9. Warning, will re-use any intermediate files found.
#[Oct 03 03:52 PM]: Parsing annotation and preparing annotation files.
#[Oct 03 03:53 PM]: Found 40,316 gene models from GFF3 annotation
#[Oct 03 03:58 PM]: Adding Functional Annotation to Mani, NCBI accession: None
#[Oct 03 03:58 PM]: Annotation consists of: 40,316 gene models
#[Oct 03 03:58 PM]: 43,067 protein records loaded
#[Oct 03 03:58 PM]: Running HMMer search of PFAM version 35.0
#[Oct 03 04:55 PM]: 31,169 annotations added
#[Oct 03 04:55 PM]: Running Diamond blastp search of UniProt DB version 2022_01
#[Oct 03 04:58 PM]: 1,242 valid gene/product annotations from 2,069 total
#[Oct 03 04:58 PM]: Running Eggnog-mapper
#[Oct 03 06:37 PM]: Parsing EggNog Annotations
#[Oct 03 06:37 PM]: EggNog version parsed as 2.1.9
#[Oct 03 06:37 PM]: 47,642 COG and EggNog annotations added
#[Oct 03 06:37 PM]: Combining UniProt/EggNog gene and product names using Gene2Product version 1.76
#[Oct 03 06:37 PM]: 9,895 gene name and product description annotations added
#[Oct 03 06:37 PM]: Running Diamond blastp search of MEROPS version 12.0
#[Oct 03 06:37 PM]: 1,050 annotations added
#[Oct 03 06:37 PM]: Annotating CAZYmes using HMMer search of dbCAN version 10.0
#[Oct 03 06:42 PM]: 389 annotations added
#[Oct 03 06:42 PM]: Annotating proteins with BUSCO arthropoda models
#[Oct 03 06:44 PM]: 1,344 annotations added
#[Oct 03 06:44 PM]: Predicting secreted and transmembrane proteins using Phobius
#[Oct 03 06:53 PM]: Predicting secreted proteins with SignalP
#[Oct 03 07:21 PM]: 2,822 secretome and 6,852 transmembane annotations added
#[Oct 03 07:21 PM]: Parsing InterProScan5 XML file
#[Oct 03 07:22 PM]: Found 0 duplicated annotations, adding 111,489 valid annotations
#[Oct 03 07:22 PM]: Converting to final Genbank format, good luck!
#[Oct 03 07:27 PM]: Creating AGP file and corresponding contigs file
#[Oct 03 07:28 PM]: Writing genome annotation table.
#[Oct 03 07:30 PM]: Funannotate annotate has completed successfully!
#
#        We need YOUR help to improve gene names/product descriptions:
#           0 gene/products names MUST be fixed, see funannotate.AUG.Mani_v1.9/annotate_results/Gene2Products.must-fix.txt
#           76 gene/product names need to be curated, see funannotate.AUG.Mani_v1.9/annotate_results/Gene2Products.need-curating.txt
#           931 gene/product names passed but are not in Database, see funannotate.AUG.Mani_v1.9/annotate_results/Gene2Products.new-names-passed.txt
#
#        Please consider contributing a PR at https://github.com/nextgenusfs/gene2product
#
#-------------------------------------------------------

# rename IDs and remove g6053 ("forSubmission")
cd funannotate.AUG.Mani_v1.9/annotate_results
sed -E 's/^(>)file_[0-9]_file_[0-9]_([^ ]*).*/\1Pmani_v1.9_\2/' Mani.proteins.fa | gzip > Mani_v1.9_221214.proteins.fa.gz
sed -E 's/^(>)file_[0-9]_file_[0-9]_([^ ]*).*/\1Pmani_v1.9_\2/' Mani.mrna-transcripts.fa | gzip > Mani_v1.9_221214.mrna-transcripts.fa.gz
sed -E 's/^(>)file_[0-9]_file_[0-9]_([^ ]*).*/\1Pmani_v1.9_\2/' Mani.cds-transcripts.fa | gzip > Mani_v1.9_221214.cds-transcripts.fa.gz
sed -E 's/file_[0-9]_file_[0-9]_([^;]*)/Pmani_v1.9_\1/g' Mani.gff3 | gzip > Mani_v1.9_221214.gff3.gz
table2asn -M n -J -c w -euk -t template.sbt -i ../../Mani_v1.9.fasta -f Mani_v1.9_221214.forSubmission.gff3 -o Mani_v1.9_221214.forSubmission.sqn -Z -locus-tag-prefix Pmani -j "[organism=Petrolisthes manimaculis] [isolate=PB745_02]"






# repeat the above steps for the second species

# creat database for rna-seq mapping
hisat2-build Cinc_v1.5.fa.masked Cinc_v1.5.fa

# map rna-seq reads to reference by executing the following for all read-sets
#sh runAlignConvertSortStats.sh $r1 $r2 . Cinc_v1.5.fa.masked

# gene prediction using braker2 supplying the contaminant reduced genome and the mapped rna-reads
braker.pl --species=Petrolisthes_cinctipes --genome=Cinc_v1.5.fa.masked --bam=EA1_12_ATCACG_L006_R1_001.fastq_sorted.bam,EA1_12_ATCACG_L006_R1_002.fastq_sorted.bam,EA1_12_ATCACG_L006_R1_003.fastq_sorted.bam,EA1_14_ACAGTG_L006_R1_001.fastq_sorted.bam,EA1_14_ACAGTG_L006_R1_002.fastq_sorted.bam,EA1_1_ACTTGA_L006_R1_001.fastq_sorted.bam,EA1_1_ACTTGA_L006_R1_002.fastq_sorted.bam,EA1_1_ACTTGA_L006_R1_003.fastq_sorted.bam,EA1_21_GATCAG_L006_R1_001.fastq_sorted.bam,EA1_21_GATCAG_L006_R1_002.fastq_sorted.bam,EA1_21_GATCAG_L006_R1_003.fastq_sorted.bam,EA1_29_CGATGT_L006_R1_001.fastq_sorted.bam,EA1_29_CGATGT_L006_R1_002.fastq_sorted.bam,EA1_29_CGATGT_L006_R1_003.fastq_sorted.bam,EA1_41_AGTTCC_L006_R1_001.fastq_sorted.bam,EA1_41_AGTTCC_L006_R1_002.fastq_sorted.bam,EA1_41_AGTTCC_L006_R1_003.fastq_sorted.bam,EA1_43_GTGAAA_L006_R1_001.fastq_sorted.bam,EA1_43_GTGAAA_L006_R1_002.fastq_sorted.bam,EA1_43_GTGAAA_L006_R1_003.fastq_sorted.bam,EA1_45_ATTCCT_L006_R1_001.fastq_sorted.bam,EA1_45_ATTCCT_L006_R1_002.fastq_sorted.bam,EA1_45_ATTCCT_L006_R1_003.fastq_sorted.bam,EA1_54_GAGTGG_L006_R1_001.fastq_sorted.bam,EA1_54_GAGTGG_L006_R1_002.fastq_sorted.bam,EA1_54_GAGTGG_L006_R1_003.fastq_sorted.bam,EA1_54_GAGTGG_L006_R1_004.fastq_sorted.bam,EA1_56_CCGTCC_L006_R1_001.fastq_sorted.bam,EA1_56_CCGTCC_L006_R1_002.fastq_sorted.bam,EA1_56_CCGTCC_L006_R1_003.fastq_sorted.bam,EA1_7_CAGATC_L006_R1_001.fastq_sorted.bam,EA1_7_CAGATC_L006_R1_002.fastq_sorted.bam,EA1_7_CAGATC_L006_R1_003.fastq_sorted.bam,EA2_16_TGACCA_L007_R1_001.fastq_sorted.bam,EA2_16_TGACCA_L007_R1_002.fastq_sorted.bam,EA2_16_TGACCA_L007_R1_003.fastq_sorted.bam,EA2_18_ATCACG_L007_R1_001.fastq_sorted.bam,EA2_18_ATCACG_L007_R1_002.fastq_sorted.bam,EA2_18_ATCACG_L007_R1_003.fastq_sorted.bam,EA2_20_ACAGTG_L007_R1_001.fastq_sorted.bam,EA2_20_ACAGTG_L007_R1_002.fastq_sorted.bam,EA2_20_ACAGTG_L007_R1_003.fastq_sorted.bam,EA2_3_TAGCTT_L007_R1_001.fastq_sorted.bam,EA2_3_TAGCTT_L007_R1_002.fastq_sorted.bam,EA2_3_TAGCTT_L007_R1_003.fastq_sorted.bam,EA2_3_TAGCTT_L007_R1_004.fastq_sorted.bam,EA2_52_GTCCGC_L007_R1_001.fastq_sorted.bam,EA2_52_GTCCGC_L007_R1_002.fastq_sorted.bam,EA2_52_GTCCGC_L007_R1_003.fastq_sorted.bam,EA2_58_CCGTCC_L007_R1_001.fastq_sorted.bam,EA2_58_CCGTCC_L007_R1_002.fastq_sorted.bam,EA2_58_CCGTCC_L007_R1_003.fastq_sorted.bam,EA2_5_CTTGTA_L007_R1_001.fastq_sorted.bam,EA2_5_CTTGTA_L007_R1_002.fastq_sorted.bam,EA2_5_CTTGTA_L007_R1_003.fastq_sorted.bam,EA2_60_AGTTCC_L007_R1_001.fastq_sorted.bam,EA2_60_AGTTCC_L007_R1_002.fastq_sorted.bam,EA2_60_AGTTCC_L007_R1_003.fastq_sorted.bam,EA2_60_AGTTCC_L007_R1_004.fastq_sorted.bam,EA2_67_CGTACG_L007_R1_001.fastq_sorted.bam,EA2_67_CGTACG_L007_R1_002.fastq_sorted.bam,EA2_67_CGTACG_L007_R1_003.fastq_sorted.bam,EA2_69_GTTTCG_L007_R1_001.fastq_sorted.bam,EA2_69_GTTTCG_L007_R1_002.fastq_sorted.bam,EA2_69_GTTTCG_L007_R1_003.fastq_sorted.bam,EA2_9_GCCAAT_L007_R1_001.fastq_sorted.bam,EA2_9_GCCAAT_L007_R1_002.fastq_sorted.bam,EA2_9_GCCAAT_L007_R1_003.fastq_sorted.bam --softmasking --cores 8 --workingdir=braker.softmask.hisat2

# merge overlapping genes and remove such with internal stop codons
conda activate agat
grep -v "transcript " braker.softmask.hisat2_v1.9.gtf | grep -v "GeneMark.hmm" > braker_v1.9.AUG.tmp.gtf

agat_convert_sp_gxf2gxf.pl -g braker_v1.9.AUG.tmp.gtf -o braker_v1.9.AUG.gff3

# remove NCBI identifyied adaptors
grep -o '^>.*' Cinc_v1.5.fa | sed 's/>//' > adaptor_regions.txt # do modifications as below

#ptg000421l:1-142362
#ptg001100l:394-
#ptg001292l:1-37608
#ptg001778l:399-
#ptg003101l:1-186530
#ptg005025l:1-67838
#ptg005762l:1-77525
#ptg006792l:436-
#ptg009405l:1-44053
#ptg009566l:1-30296
#ptg010060l:1-124103
#ptg013271l:1266-
#ptg017086l:86-
#ptg017777l:1-46808
#ptg018034l:1-44528
#ptg019000l:1-25793
#ptg019013l:1-26315
#ptg021049l:1-23519
#ptg007836l:1-13258
#ptg007836l:13286-

mkdir tmp
cd tmp
while IFS= read -r region; do samtools faidx ../Cinc_v1.5.fa "$region" > "$region.fa"; done < ../adaptor_regions.txt
cat *.fa > Cinc_v1.5.1.fasta
sed -i 's/ptg007836l:13286/ptg100000l:13286/' Cinc_v1.5.1.fasta
sed -i -E 's/^(>ptg[0-9l]*):.*$/\1/' Cinc_v1.5.1.fasta
mv Cinc_v1.5.1.fasta ../Cinc_v1.10.fasta
cd ../
rm -r tmp

awk -F'\t' '!/^#/ {OFS="\t"; if($1 == "ptg001100l" && $4 > 394) $4 -= 393;   if($1 == "ptg001100l" && $5 > 394) $5 -= 393;
                             if($1 == "ptg001778l" && $4 > 399) $4 -= 398;   if($1 == "ptg001778l" && $5 > 399) $5 -= 398;
                             if($1 == "ptg006792l" && $4 > 436) $4 -= 435;   if($1 == "ptg006792l" && $5 > 436) $5 -= 435;
                             if($1 == "ptg013271l" && $4 > 1266) $4 -= 1265; if($1 == "ptg013271l" && $5 > 1266) $5 -= 1265;
                             if($1 == "ptg017086l" && $4 > 86) $4 -= 85;     if($1 == "ptg017086l" && $5 > 86) $5 -= 85;
                             if($1 == "ptg007836l" && $5 > 13285) $1 = "ptg100000l"; if($1 == "ptg100000l" && $4 > 13285) $4 -= 13285;  if($1 == "ptg100000l" && $5 > 13285) $5 -= 13285;
                             ; print}' braker_v1.9.AUG.gff3 > braker_v1.10.AUG.gff3 # modify start of ptg100000l

# create multi line reference and extract protein sequences
fold -w 60 Cinc_v1.10.fasta > Cinc_v1.10.fold.fa
agat_sp_extract_sequences.pl --gff braker_v1.10.AUG.gff3 -f Cinc_v1.10.fold.fa -p -o braker_v1.10.AUG.proteins.fasta

# functional annotation using InterProScan5 (create single line fasta and remove terminal stop codon first)
conda activate interproscan
sed -E 's/^(>.*) gene=(.*) seq.*/\1 \2/' braker_v1.10.AUG.proteins.fasta | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' | sed -e 's/\*$//' > Cinc_v1.10.AUG.proteins.fasta
python ~/Ht_methylation_project/analysis/iso-seq/long_read_protocol/iprscan-local.py -c 6 -i Cinc_v1.10.AUG.proteins.fasta -m local -o funannotate.AUG.Cinc_v1.10/annotate_misc/iprscan.xml
conda deactivate

# functional annotation using PFAM, UniProt, EggNog, MEROPS, dbCAN, BUSCO, Phobius, SignalP, and InterProScan5 in funannotate
conda activate mamba-funannotate
export AUGUSTUS_CONFIG_PATH=$HOME/miniconda3/envs/mamba-funannotate/config/
export FUNANNOTATE_DB=$HOME/bioinformatics/funannotate
export GENEMARK_PATH=$HOME/miniconda3/envs/mamba-funannotate/bin/
export TRINITYHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/trinity-2.8.5
export EVM_HOME=$HOME/miniconda3/envs/mamba-funannotate/opt/evidencemodeler-1.1.1
export PASAHOME=$HOME/miniconda3/envs/mamba-funannotate/opt/pasa-2.4.1

funannotate annotate --gff braker_v1.10.AUG.gff3 --fasta Cinc_v1.10.fasta -s "Cinc" -o funannotate.AUG.Cinc_v1.10 --cpus 8 --busco_db arthropoda --iprscan funannotate.AUG.Cinc_v1.10/annotate_misc/iprscan.xml
#-------------------------------------------------------
#[Oct 04 08:31 AM]: OS: CentOS Linux 7, 24 cores, ~ 264 GB RAM. Python: 3.8.12
#[Oct 04 08:31 AM]: Running 1.8.11
#[Oct 04 08:31 AM]: No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'
#[Oct 04 08:31 AM]: Found existing output directory funannotate.AUG.Cinc_v1.10. Warning, will re-use any intermediate files found.
#[Oct 04 08:31 AM]: Parsing annotation and preparing annotation files.
#[Oct 04 08:32 AM]: Found 44,545 gene models from GFF3 annotation
#[Oct 04 08:39 AM]: Adding Functional Annotation to Cinc, NCBI accession: None
#[Oct 04 08:39 AM]: Annotation consists of: 44,545 gene models
#[Oct 04 08:39 AM]: 47,633 protein records loaded
#[Oct 04 08:39 AM]: Existing Pfam-A results found: funannotate.AUG.Cinc_v1.10/annotate_misc/annotations.pfam.txt
#[Oct 04 08:39 AM]: 33,584 annotations added
#[Oct 04 08:39 AM]: Running Diamond blastp search of UniProt DB version 2022_01
#[Oct 04 08:39 AM]: 1,388 valid gene/product annotations from 2,204 total
#[Oct 04 08:39 AM]: Existing Eggnog-mapper results found: funannotate.AUG.Cinc_v1.10/annotate_misc/eggnog.emapper.annotations
#[Oct 04 08:39 AM]: Parsing EggNog Annotations
#[Oct 04 08:39 AM]: EggNog version parsed as 2.1.9
#[Oct 04 08:39 AM]: 53,963 COG and EggNog annotations added
#[Oct 04 08:39 AM]: Combining UniProt/EggNog gene and product names using Gene2Product version 1.76
#[Oct 04 08:39 AM]: 10,683 gene name and product description annotations added
#[Oct 04 08:39 AM]: Existing MEROPS results found: funannotate.AUG.Cinc_v1.10/annotate_misc/annotations.merops.txt
#[Oct 04 08:39 AM]: 1,052 annotations added
#[Oct 04 08:39 AM]: Existing CAZYme results found: funannotate.AUG.Cinc_v1.10/annotate_misc/annotations.dbCAN.txt
#[Oct 04 08:39 AM]: 406 annotations added
#[Oct 04 08:39 AM]: Existing BUSCO2 results found: funannotate.AUG.Cinc_v1.10/annotate_misc/annotations.busco.txt
#[Oct 04 08:39 AM]: 1,378 annotations added
#[Oct 04 08:39 AM]: Existing Phobius results found: funannotate.AUG.Cinc_v1.10/annotate_misc/phobius.results.txt
#[Oct 04 08:39 AM]: Predicting secreted proteins with SignalP
#[Oct 04 09:01 AM]: 3,128 secretome and 7,498 transmembane annotations added
#[Oct 04 09:01 AM]: Parsing InterProScan5 XML file
#[Oct 04 09:01 AM]: Found 0 duplicated annotations, adding 122,817 valid annotations
#[Oct 04 09:01 AM]: Converting to final Genbank format, good luck!
#[Oct 04 09:07 AM]: Creating AGP file and corresponding contigs file
#[Oct 04 09:09 AM]: Writing genome annotation table.
#[Oct 04 09:09 AM]: putative transcript from ncbi:file_1_file_1_g2600.t1 has no ID
#(ncbi:file_1_file_1_g2600.t1 None ncbi:file_1_file_1_g2600.t1)
#[Oct 04 09:09 AM]: putative transcript from ncbi:file_1_file_1_g36486.t1 has no ID
#(ncbi:file_1_file_1_g36486.t1 None ncbi:file_1_file_1_g36486.t1)
#[Oct 04 09:11 AM]: Funannotate annotate has completed successfully!
#
#        We need YOUR help to improve gene names/product descriptions:
#           0 gene/products names MUST be fixed, see funannotate.AUG.Cinc_v1.10/annotate_results/Gene2Products.must-fix.txt
#           75 gene/product names need to be curated, see funannotate.AUG.Cinc_v1.10/annotate_results/Gene2Products.need-curating.txt
#           926 gene/product names passed but are not in Database, see funannotate.AUG.Cinc_v1.10/annotate_results/Gene2Products.new-names-passed.txt
#
#        Please consider contributing a PR at https://github.com/nextgenusfs/gene2product
#
#-------------------------------------------------------

# rename IDs and remove g2600 and g36486 ("forSubmission")
cd funannotate.AUG.Cinc_v1.10/annotate_results
sed -E 's/^(>)file_[0-9]_file_[0-9]_([^ ]*).*/\1Pcinc_v1.10_\2/' Cinc.proteins.fa | gzip > Cinc_v1.10_221214.proteins.fa.gz
sed -E 's/^(>)file_[0-9]_file_[0-9]_([^ ]*).*/\1Pcinc_v1.10_\2/' Cinc.mrna-transcripts.fa | gzip > Cinc_v1.10_221214.mrna-transcripts.fa.gz
sed -E 's/^(>)file_[0-9]_file_[0-9]_([^ ]*).*/\1Pcinc_v1.10_\2/' Cinc.cds-transcripts.fa | gzip > Cinc_v1.10_221214.cds-transcripts.fa.gz
sed -E 's/file_[0-9]_file_[0-9]_([^;]*)/Pcinc_v1.10_\1/g' Cinc.gff3 | gzip > Cinc_v1.10_221214.gff3.gz
table2asn -M n -J -c w -euk -t template.sbt -i ../../Cinc_v1.10.fasta -f Cinc_v1.10_221214.forSubmission.gff3 -o Cinc_v1.10_221214.forSubmission.sqn -Z -locus-tag-prefix Pcinc -j "[organism=Petrolisthes cinctipes] [isolate=PB745_01]"
