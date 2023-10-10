# get divergence estimates from single copy orthologs
cd /home/pascal/crab/hifiasm-assembly-Eric/dnds

# prepare input
cp ../Cinc/HIFIASM/funannotate_Cinc_v1.9/annotate_results/Cinc_v1.9_221214.cds-transcripts.fa.gz .
gunzip Cinc_v1.9_221214.cds-transcripts.fa.gz 
cp ../Mani/HIFIASM/funannotate_Mani_v1.7/annotate_results/Mani_v1.7_221214.cds-transcripts.fa.gz .
gunzip Mani_v1.7_221214.cds-transcripts.fa.gz

# cut orthologs
for i in `awk '{print $2}' <(awk 'NF==3' ../genespace/orthofinder/OrthoFinder/Results_Dec16/Orthogroups/Orthogroups.txt | grep Pcinc  | grep Pmani)`; do awk 'BEGIN{RS=">"} /'$i'/{print ">" $0}' Cinc_v1.9_221214.cds-transcripts.fa ;done > Cinc_v1.9_221214.orthologs.fasta
for i in `awk '{print $3}' <(awk 'NF==3' ../genespace/orthofinder/OrthoFinder/Results_Dec16/Orthogroups/Orthogroups.txt | grep Pcinc  | grep Pmani)`; do awk 'BEGIN{RS=">"} /'$i'/{print ">" $0}' Mani_v1.7_221214.cds-transcripts.fa ;done > Mani_v1.7_221214.orthologs.fasta

# translate gene IDs from Mani to Cinc
paste <(awk '{print $2}' <(awk 'NF==3' ../genespace/orthofinder/OrthoFinder/Results_Dec16/Orthogroups/Orthogroups.txt | grep Pcinc  | grep Pmani)) <(awk '{print $3}' <(awk 'NF==3' ../genespace/orthofinder/OrthoFinder/Results_Dec16/Orthogroups/Orthogroups.txt | grep Pcinc  | grep Pmani)) | while read b a; do sed -i "s/$a/$b/" Mani_v1.7_221214.orthologs.fasta; done

# make one file for each Cinc CDS
mkdir -p fasta_references/fasta_{Cinc,Mani}
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' Cinc_v1.9_221214.orthologs.fasta
for f in *cinc*; do cp "$f" "$(echo "fasta_references/fasta_Cinc/$f" | sed -e 's/>//')";done

# add the corresponding Mani CDS
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' Mani_v1.7_221214.orthologs.fasta
for f in *cinc*; do mv "$f" "$(echo "fasta_references/$f" | sed -e 's/>//;s/.fasta//')";done # as input list for Ka/Ks calculation

# make one file per CDS with just the Mani version
awk -F "|" '/^>/ {close(F) ; F = $1".fasta"} {print >> F}' Mani_v1.7_221214.orthologs.fasta
for f in *cinc*; do mv "$f" "$(echo "fasta_references/fasta_Mani/$f" | sed -e 's/>//')";done

# estimate divergence, mask poorly aligned sequences, and estimate divergence again
cd fasta_references
ls *cinc* | parallel -j 1 "Rscript ../originals/EstimateKaKs2_modified.Rscript {}"
ls *aligned* | parallel -j 6 "Rscript ../originals/mask_alignment.Rscript {}"
ls *aligned_masked* | parallel -j 1 "Rscript ../originals/mask_paml.Rscript {}"


# visualize gene families

library(VennDiagram)
library(cowplot)

p_venn <- draw.pairwise.venn(
    22433+13422, 20376+13422, 13422,
    c(expression(italic("P. cinctipes")), expression(italic("P. manimaculis"))),
    fill = c("white", "white"),
    alpha = c(0.7, 0.7),
    ind = FALSE, cex = 3, cat.cex = 3, cat.pos = c(335, 25), cat.dist = 0.06
)

ggdraw(p_venn) +
    theme(
        plot.background = element_rect(fill = NA),
        plot.margin = margin(12, 12, 12, 12)
    )
