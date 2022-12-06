#############
### XAUT2 ###
#############

# biopsies
biopsies_annotations <- c(
"0" = "CD4T",                      # CD4(med), IL7R+, CD3E+, KLRB1+
"1" = "Plasmablast",               # IGKV+, IGLV+
"2" = "CD8T",                      # CD8+ (CITE and GEX), CCL4/5+, IFNG+
"3" = "Epithelial",                # PHGR1+, KRT8+, FABP1+
"4" = "Mesenchymal Stromal",       # COL3A1+, COL1A2+
"5" = "Macrophage",                # C1QA/B+
"6" = "B",                         # MS4A1+, CD74+, CD83+
"7" = "Cycling",                   # HIST1H4C+, TUBA1B+, STMN1+
"8" = "Endothelial",               # SPARCL1+, PECAM1+, VWF+
"9" = "Glial"                      # CLU+, CRYAB+, SPARC+
)

# Markers from Martin et al., Fawkner et al., Corridoni et al., and Boland et al.
biopsies_T <- c(
"0" = "CD4+ RM",                   # CD69+, FOS+, FOSB+, KLF6+, JUN+
"1" = "CD8+ Effector 1",           # GZMA+, IFNG+, ITGAE+, CD103 (CITE)
"2" = "CD4+ Tregs/Naive",          # CTLA4+, FOXP3+, SELL+, CCR7+
"3" = "CD8+ Effector 2",           # GZMK+, CCL4/5+
"4" = "CD4-/CD8- NK",              # CD56+ (CITE), TRDC+, IL2RB+, NKG7+, TYROBP+
"5" = "CD4+ Memory",               # CD45RO+ (CITE), CD45RA- (CITE)
"6" = "CD4-/CD8- Cycling"          # IL7R+, STMN1+, TUBA1B+
)

biopsies_T <- c(
"0" = "CD8+ Effector 1",                     # IFNG+, GZMB+, CCL5+, GZMA+
"1" = "CD4+ RM",                             # CD69+, FOS+, FOSB+, KLF6+, JUN+
"2" = "CD4+ Naive",                          # SELL+, CCR7+, LEF1+, TCF7+
"3" = "CD4+ RM",                             # CD69+, FOS+, FOSB+, KLF6+, JUN+
"4" = "CD4+ Tregs",                          # CTLA4+, BATF+, FOXP3+
"5" = "CD8+ Effector 2",                     # GZMK+, CCL4+, CCL5+, GZMA+
"6" = "CD4+ Tregs",                          # CTLA4+, BATF+, FOXP3+
"7" = "CD8+ Effector 1",                     # IFNG+, GZMB+, CCL5+, GZMA+
"8" = "CD8+ Effector 1",                     # IFNG+, GZMB+, CCL5+, GZMA+
"9" = "T.D. Gamma Delta",                    # TRDC+ TRDV1+, TRDV2+, TRGV9+
"10" = "Activated Gamma Delta",              # TRDC+, TRGV9+, GZMB+, CCL4+, CCL3+
"11" = "CD4+ Memory",                        # CD45RO+ (CITE), CD45RA- (CITE)
"12" = "Memory Gamma Delta"                  # TRDC+, IL7R+
)

biopsies_T_marks <- c("IL7R", "TRDV1", "TRDV2", "TRGV9", "TRDC", "GZMB", "CCL4", "CCL3", "IFNG", "CCL5", "GZMA", "GZMK", "ADT_CD45RO", "ADT_CD45RA", "SELL", "CCR7", "LEF1", "TCF7", "CD69", "FOS", "FOSB", "KLF6", "JUN", "CTLA4", "BATF", "FOXP3")

# All markers from Dr. Kattah
biopsies_MES <- c(
"0" = "S1",                        # CCL8+, TMSB4X+, CTSC+ (all S1 markers)
"1" = "S2",                        # POSTN+, VSTM2A+, COL4A6+ (all S2 markers)
"2" = "MF",                        # FLNA+, MYH11+, ACTA2+ (all MF markers)
"3" = "S3/S4",                     # NFIA+, GSN+, TNFSF13B+, CD74+ (most S3 and S4 markers)
"4" = "Pericytes"                  # CSPG4+, RGS5+
)

# All markers from Parikh et al., Fawkner et al., and Dr. Kattah
biopsies_EPI <- c(
"0" = "Colonocytes",               # CA1+, PHGR1+, SELENBP1+, S100A6+, CEACAM1/7-, AQP8-
"1" = "CT Colonocytes",            # AQP8+, SLC26A3+, CEACAM1/7+
"2" = "Undifferentiated 1",        # S100A6+, ADH1C+, TMSB10+
"3" = "BEST4/OTOP2 Cells",         # BEST4+, OTOP2+, CA7+,
"4" = "Goblets",                   # SPINK4+, MUC2+, TFF3+
"5" = "EECs",                      # CRYBA2+, SCGN+, PCSK1N+
"6" = "Undifferentiated 2"         # TMSB4X+, EEF1A1+, B2M+
)

# Markers from Dr. Kattah or Fawkner et al.
biopsies_ENDO <- c(
"0" = "Blood Endothelial",         # PECAM1+, VWF+, FLT1+, CD36+
"1" = "Venous",                    # ACKR1+, VWF+, PRCP+, JUNB-, ZFP36-
"2" = "Lymphatic"                  # LYVE1+, CCL21+, PROX1+, RELN+
)

# Markers from Boland et al. and Zilionis et al.
biopsies_B <- c(
"0" = "Naive 1",                   # BANK1+
"1" = "Naive 2",                   # CD19+ (GEX), IGHD+
"2" = "Plasmablast",               # IGKV+, IGHA+, IGHV+, JCHAIN+
"3" = "Memory"                     # CD19+ (CITE), IGHD-
)

# pbmcs
pbmcs_annotations <- c(
"0" = "CD4T",                      # CD4+ (CITE), TCR+, IL7R+, CD3E+
"1" = "Classical Monocytes",       # CD14+, FCGR3A-, S100A8/9+ 
"2" = "CD4T",                      # CD4+ (CITE), TCR+, IL7R+, CD3E+
"3" = "NK",                        # CD4-/CD8- (CITE), NKG7+, FCGR3A+ 
"4" = "B",                         # MS4A1+, BCR+, CD79A+
"5" = "CD8T",                      # CD8+ (CITE), TCR+, CD3E+
"6" = "Classical Monocytes",       # CD14+, FCGR3A-, S100A8/9+ 
"7" = "CD4T",                      # CD4+ (CITE), TCR+, IL7R+, CD3E+
"8" = "Classical Monocytes",       # CD14+, FCGR3A-, S100A8/9+ 
"9" = "Platelets",                 # PPBP+, PF4+
"10" = "Nonclassical Monocytes",   # FCGR3A+, CD14-
"11" = "APCs",                     # HLA-DPB1+, HLA-DQA1+, CD74+
"12" = "Neutrophils",              # FCGR3B+, CD16+ (CITE only), low nCount_RNA, CD19-, CD3-
"13" = "Classical Monocytes",      # CD14+, FCGR3A-, S100A8/9+
"14" = "Plasmablast",              # JCHAIN+, IGKC+, IGHA1+
"15" = "pDCs",                     # GZMB+, IRF7+, CD123+ (CITE), PTCRA+, LILRA4+, IL3RA+
"16" = "Cycling"                   # STMN1+, TUBB+, HIST1H4C+, TUBA1B+
)

# Markers from Boland et al., Martin et al., and Fawkner et al.
pbmcs_B <- c(
"0" = "Naive 1",                   # CD19+, IGHD+, FCER2+, CD79B(low)
"1" = "Memory",                    # CD19+, IGHD-, IGHG1+
"2" = "Naive 2",                   # CD19+, IGHD+, FCER2+, CD79B(high)
"3" = "Plasmablast"                # IGKV+, IGLV+
)

pbmcs_T <- c(
"0" = "CD4+ Memory",                      # CD45RO+ (CITE), CD45RA- (CITE)
"1" = "CD4+ Naive 1",                     # SELL+, CCR7+, LEF1+, TCF7+
"2" = "CD4+ Naive 1",                     # SELL+, CCR7+, LEF1+, TCF7+
"3" = "CD4+ Naive 1",                     # SELL+, CCR7+, LEF1+, TCF7+
"4" = "CD4+ Memory",                      # CD45RO+ (CITE), CD45RA- (CITE)
"5" = "CD8+ Naive",                       # SELL+, CCR7+, LEF1+, TCF7+
"6" = "CD4+ Naive 1",                     # SELL+, CCR7+, LEF1+, TCF7+
"7" = "CD4+ Memory",                      # CD45RO+ (CITE), CD45RA- (CITE)
"8" = "CD4+ RM",                          # CD69+, FOS+, FOSB+, JUN+, KLF6+
"9" = "CD8+ Effector",                    # NKG7+, GZMB+, CCL4/5+
"10" = "CD4+ RM",                         # CD69+, FOS+, FOSB+, JUN+, KLF6+
"11" = "CD4+ Naive 2",                    # IL7R+ (CITE only), TNFAIP3+, NFKBIA+, TSC22D3+ (last three genes: anti-inflammatory related functions)
"12" = "CD8+ Effector",                   # NKG7+, GZMB+, CCL4/5+
"13" = "CD8+ Effector",                   # NKG7+, GZMB+, CCL4/5+
"14" = "CD4-/CD8- Gamma Delta",           # TRDC+, TRDV2+, TRGV9+, NKG7+
"15" = "CD8+ Naive",                      # SELL+, CCR7+, LEF1+, TCF7+
"16" = "CD4+ Memory",                     # CD45RO+ (CITE), CD45RA- (CITE)
"17" = "CD8+ Naive",                      # SELL+, CCR7+, LEF1+, TCF7+
"18" = "CD4-/CD8- Cycling",               # MKI67+, TUBA1B+, STMN1+
"19" = "CD4+ Naive 1"                     # SELL+, CCR7+, LEF1+, TCF7+
)

pbmcs_T_marks <- c("CD8A","adt_CD4.2","STMN1", "TUBA1B", "MKI67", "TRDC","TRDV1","TRDV2","TRGV9", "adt_CD45RO", "adt_CD45RA", "SELL", "CCR7", "LEF1", "TCF7", "adt_CD127--IL-7Ra", "TNFAIP3", "NFKBIA", "TSC22D3", "CD69", "FOS", "FOSB", "JUN", "KLF6", "NKG7", "GZMB", "CCL4", "CCL5")

pbmcs_MONO <- c(
"0" = "Classical Monocytes",              # CD14+
"1" = "Classical Monocytes",              # CD14+
"2" = "Classical Monocytes",              # CD14+
"3" = "Classical Monocytes",              # CD14+
"4" = "Classical Monocytes",              # CD14+
"5" = "Classical Monocytes",              # CD14+
"6" = "Nonclassical Monocytes",           # CD16+, CD14-
"7" = "Classical Monocytes",              # CD14+
"8" = "Neutrophils"                       # FCGR3B+, CD14-, FCGR3A-
)

pbmcs_MONO_marks <- c("CD14","FCGR3B","FCGR3A")