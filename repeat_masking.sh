## Repeat masking for gene prediction of contaminant reduced crab assemblies

# remove contig which does not contain any G nucleotides (i.e., ATC only)
cat Mani_v1.6.fa | paste - - | grep -v -F ptg016427l | tr "\t" "\n" > Mani_v1.7.fa

# de novo repeat identification using RepeatModeler
BuildDatabase -name Mani Mani_v1.7.fa
RepeatModeler -database Mani -pa 24 -LTRStruct

# annotation and masking using RepeatMasker
RepeatMasker -pa 24 -lib Mani-families.fa Mani_v1.7.fa -xsmall

#==================================================
#file name: Mani_v1.7.fa
#sequences:          7486
#total length:  923773821 bp  (923773821 bp excl N/X-runs)
#GC level:         38.24 %
#bases masked:  542843944 bp ( 58.76 %)
#==================================================
#               number of      length   percentage
#               elements*    occupied  of sequence
#--------------------------------------------------
#Retroelements       251001    100201730 bp   10.85 %
#   SINEs:                0            0 bp    0.00 %
#   Penelope           5898      1578852 bp    0.17 %
#   LINEs:           192945     75638857 bp    8.19 %
#    CRE/SLACS           82        59106 bp    0.01 %
#     L2/CR1/Rex      80012     26418587 bp    2.86 %
#     R1/LOA/Jockey   25188     20240227 bp    2.19 %
#     R2/R4/NeSL       1297      1072855 bp    0.12 %
#     RTE/Bov-B       57721     17970553 bp    1.95 %
#     L1/CIN4          1082       270928 bp    0.03 %
#   LTR elements:     58056     24562873 bp    2.66 %
#     BEL/Pao             0            0 bp    0.00 %
#     Ty1/Copia       22421      3488085 bp    0.38 %
#     Gypsy/DIRS1     35279     21026105 bp    2.28 %
#       Retroviral      356        48683 bp    0.01 %
#
#DNA transposons      78739     25485922 bp    2.76 %
#   hobo-Activator    14775      3311690 bp    0.36 %
#   Tc1-IS630-Pogo    26736     10293149 bp    1.11 %
#   En-Spm                0            0 bp    0.00 %
#   MuDR-IS905            0            0 bp    0.00 %
#   PiggyBac           4075      3425279 bp    0.37 %
#   Tourist/Harbinger  2970      1758172 bp    0.19 %
#   Other (Mirage,     6819      2465531 bp    0.27 %
#    P-element, Transib)
#
#Rolling-circles       1630       265603 bp    0.03 %
#
#Unclassified:       1286958    329235969 bp   35.64 %
#
#Total interspersed repeats:   454923621 bp   49.25 %
#
#
#Small RNA:            9658      1992403 bp    0.22 %
#
#Satellites:           1114       148604 bp    0.02 %
#Simple repeats:     1020918     77734145 bp    8.41 %
#Low complexity:      75865      7779568 bp    0.84 %
#==================================================

# repeat the above steps for the second species

# de novo repeat identification using RepeatModeler
BuildDatabase -name Cinc Cinc_v1.5.fa
RepeatModeler -database Cinc -pa 24 -LTRStruct

# annotation and masking using RepeatMasker
RepeatMasker -pa 24 -lib Cinc-families.fa Cinc_v1.5.fa

#==================================================
#file name: Cinc_v1.5.fa
#sequences:          9390
#total length: 1491129859 bp  (1491129859 bp excl N/X-runs)
#GC level:         39.63 %
#bases masked: 1115386838 bp ( 74.80 %)
#==================================================
#               number of      length   percentage
#               elements*    occupied  of sequence
#--------------------------------------------------
#Retroelements       256439    202879600 bp   13.61 %
#   SINEs:                0            0 bp    0.00 %
#   Penelope          12476      5778472 bp    0.39 %
#   LINEs:           209515    149410244 bp   10.02 %
#    CRE/SLACS          428       866820 bp    0.06 %
#     L2/CR1/Rex      64027     38415218 bp    2.58 %
#     R1/LOA/Jockey   42758     61460597 bp    4.12 %
#     R2/R4/NeSL       1461       915305 bp    0.06 %
#     RTE/Bov-B       52879     25024372 bp    1.68 %
#     L1/CIN4            20         2113 bp    0.00 %
#   LTR elements:     46924     53469356 bp    3.59 %
#     BEL/Pao          1719      1775059 bp    0.12 %
#     Ty1/Copia         580       922309 bp    0.06 %
#     Gypsy/DIRS1     42230     49701106 bp    3.33 %
#       Retroviral     2395      1070882 bp    0.07 %
#
#DNA transposons      92610     29283451 bp    1.96 %
#   hobo-Activator    10113      2858343 bp    0.19 %
#   Tc1-IS630-Pogo    20527      5664433 bp    0.38 %
#   En-Spm                0            0 bp    0.00 %
#   MuDR-IS905            0            0 bp    0.00 %
#   PiggyBac           2261      2582346 bp    0.17 %
#   Tourist/Harbinger  5922      4384674 bp    0.29 %
#   Other (Mirage,    18406      4784824 bp    0.32 %
#    P-element, Transib)
#
#Rolling-circles       3964       986700 bp    0.07 %
#
#Unclassified:       1384990    806688934 bp   54.10 %
#
#Total interspersed repeats:  1038851985 bp   69.67 %
#
#
#Small RNA:            4811      1737726 bp    0.12 %
#
#Satellites:           2127       250510 bp    0.02 %
#Simple repeats:     799000     67321864 bp    4.51 %
#Low complexity:      57159      6238053 bp    0.42 %
#==================================================
