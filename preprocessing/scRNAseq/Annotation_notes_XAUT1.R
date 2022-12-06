#############
### XAUT1 ###
#############

# Coarse annotations
biopsies_annotations <- c(
"0" = "CD4T",                 # CD4+ (CITE), IL7R+, TRAC+
"1" = "CD8T",                 # CD8+ (CITE), TRAC+, CCL4/5+, GZMA+
"2" = "Plasmablast",          # IGKC+, IGHA1+, IGLC2+
"3" = "Mito High",            # MT-COS2/3+ (no DEGs)
"4" = "Mesenchymal Stromal",  # COL1A2+, VIM+, COL3A1+, and MSC phenotype confirmed with CITEseq (Horwitz et al., Curr. Opin. Hematol, 2006)
"5" = "B",                    # MS4A1+, BANK1+, CD19+ (CITE), IgM+ (CITE)
"6" = "Epithelial Mito High", # PIGR+, KRT8+, MT-COS2/3+
"7" = "Plasmablast",          # IGKC+, IGHA1+, IGLC2+
"8" = "Myeloid",              # LYZ+, CST3+, CD74+
"9" = "Epithelial",           # PIGR+, KRT8+, EpCAM+ (also enterocyte markers FABP1+, CKB+, AQP8+)
"10" = "Endothelial",         # PECAM1+, VWF+, SPARCL1+, PLVAP+
"11" = "Cycling",             # TUBA1B+, HMGB2+, 
"12" = "NK",                  # KLRB1+, IL7R+
"13" = "Mast",                # GATA2+, TPSAB1+, TPSB2+
"14" = "Glial"                # S100A4+, FGL2+, GFRA3+
)

# Mesenchymal Stromal subcompartment (all markers from Kinchen et al. unless otherwise noted)
biopsies_MES <- c(
"0" = "S1",                      # CCL8+, ADAMDEC1+, APOE+
"1" = "S2",                      # POSTN+, SOX6+, FOXF1+, VSTM2A+, COL4A6+
"2" = "Activated Fibroblasts",   # MMP1+, MMP3+, LUM+, CHI3L1+, PDPN+ (only Martin et al.)
"3" = "S3",                      # EFEMP1+, CCDC80+, OGN+, GSN+
"4" = "MF",                      # ACTG2+, ACTA2+, MYH11+, MYOCD+
"5" = "Pericyte"                 # RGS5+, CSPG4+, NDUFA4L2+ (Kinchen et al., Martin et al.)
)

# T subcompartment
biopsies_T <- c(
"0" = "CD4+ Resident Memory",            # IL7R+, CD69+, NFKBIA+, KLF6+, JUN+, KLRB1+, FOS+
"1" = "CD4+ Naive",                      # SOX4+, SC5D+, CCR7+, LEF1+, SELL+
"2" = "CD8+ Resident Memory",            # CD8A+, CD69+, ITGA1+ NKG7+, KLF6+, JUN+, FOS+
"3" = "CD8+ Effector",                   # GZMK+, IFNG+, CCL4+, CD8A+, CCL5+, CST7+
"4" = "CD4+ Tregs",                      # CTLA4+, TNFRSF4+, IL32+, FOXP3+, IL2RA+
"5" = "CD8-/CD4- Gamma Deltas"           # CD3E+, TRDC+, TRGV9+, GDTCR+ (CITE)
)

biopsies_T_marks <- c("ADT_CD4.2","CD8A", "CD3E", "SOX4", "SC5D", "CCR7", "LEF1", "SELL", "IL7R", "CD69", "NFKBIA", "KLF6", "JUN", "KLRB1", "FOS", "CTLA4", "TNFRSF4", "IL32", "FOXP3", "TRDC", "TRGV9", "ADT_TCR-G--TCR-D", "CST7", "NKG7", "IFNG", "CCL4", "CCL5", "ITGA1")

# B subcompartment
biopsies_B <- c(
"0" = "Memory",               # CXCR3-6+ (CITE, abcam website), GPR183+
"1" = "Naive 1",              # IGHD+, CD19+, IgD+ (CITE), IgM+ (CITE) (Boland et al., abcam website)
"2" = "Naive 2",              # IGHD+, CD19+ (Boland et al.)
"3" = "Cycling Memory",       # IGHD-, CD19+ (GEX and CITE), GAPDH+ (Boland et al.)
"4" = "Cycling Plasmablasts"  # XBP1+, RAN+ (transcription factor), EIF5A+ (translation initiation factor), GAPDH+ (Boland et al.)
)

# Epithelial subcompartment
biopsies_EPI <- c(
"0" = "Undifferentiated 1",       # ADH1C+, UGT2B17+, EEF1A1[high] (Parikh et al.)
"1" = "Undifferentiated 2",       # CA1+, SELENBP1[high] (Parikh et al.)
"2" = "Enterocytes",              # FABP1+, AQP8+, IL32+ (Parikh et al., Mike's markers)
"3" = "Goblet",                   # SPINK4+, TFF3+, MUC2+ Parikh et al., Mike's markers)
"4" = "Cycling Undifferentiated", # STMN1+, HMGN2+, CENPM+, TACC3+, HMGB3+ (small amount of transit activating (TA) marker expression, Fawkner et al.)
"5" = "BEST4/OTOP2 Cells",        # BEST4+, OTOP2+, CA7+ (Parikh et al., Mike's markers)
"6" = "EECs"                      # PTMS+, PCSK1N+, SCGN+, CRYBA2+ (Parikh et al., Mike's markers)
)

# Endothelial
biopsies_ENDO <- c(
"0" = "Mesenchymal",         # COL3A1+, RGS5+, ACTA2+ (Mike, mesenchymal stromal markers)
"1" = "Blood Endothelial",   # PECAM1+, VWF+, FLT1+, CD36+ (Mike, Fawkner et al.)
"2" = "Venous",              # MADCAM1+, ACKR1+, PRCP+ (Fawkner et al.)
"3" = "Lymphatic",           # LYVE1+, PROX1+, CCL21+ (Mike, Fawkner et al.)
"4" = "Arterial"             # HEY1+, GJA5+, IGFBP3+ (Fawkner et al.)
)

# Myeloid
biopsies_MYE <- c(
"0" = "MoMacDC",          # APOE+, CTSD+, ACP5+ (Zilionis et al.)
"1" = "Macrophages",      # FGL2+, HES1+, SLC40A1+, MRC1+ (Fawkner et al., Martin et al.)
"2" = "Monocytes",        # S100A8/S100A9+, FCN1+, CLEC12A (low) (Fawkner et al.)
"3" = "DCs",              # MARCKSL1+, BIRC3+, TXN+, CCL19+, LAMP3+ (Zilionis et al., Martin et al.)
"4" = "pDCs"              # GZMB+, JCHAIN+, IRF7+, CD123+ (CITE), PTCRA+, LILRA4+, IL3RA+ (Fawkner et al., Martin et al., Zilionis et al.)
)

### Leukocytes ###

# Coarse annotations
pbmcs_annotations <- c(
"0" = "CD4T",                    # IL7R+, IL32+, CD4+ (not shown)
"1" = "Classical monocytes",     # CD14+, FCGR3A-, S100A8/9+
"2" = "CD8T",                    # CD8B+, GZMK+, CCL5+
"3" = "NK",                      # NKG7+, GNLY+
"4" = "B cells",                 # MS4A1+, IGHM+, HLA-DRA+
"5" = "Nonclassical monocytes",  # CD14-, FCGR3A+
"6" = "Platelets",               # PF4+, PPBP+
"7" = "Plasmablast",             # JCHAIN+, IGKC+
"8" = "Classical monocytes",     # FCGR3B- (not neutrophils), LYZ+, CST3+
"9" = "Cycling"                  # TUBA1B+, HIST1H4C+, HMGN2+
)

# T subcompartment
pbmcs_T <- c(
"0" = "CD4+ Naive",             # IL7R+, KLF6+, FOS+ (Martin et al.)
"1" = "CD4+ Naive",             # IL7R+, KLF6+, FOS+ (Martin et al.)
"2" = "CD8+ Effector",          # PRF1+, GZMA+, GMZB+, CCL4+ (Martin et al.)
"3" = "CD8+ Effector",          # PRF1+, GZMA+, GMZK+, CCL4+ (Martin et al.)
"4" = "CD8+ MAIT",              # NCR3+, KLRB1+, TRAV1-2+ (Corridoni et al.)
"5" = "CD4+ Tregs",             # FOXP3+, CTLA4+, IL2RA+ (Martin et al.)
"6" = "CD4-/CD8- Gamma Deltas"  # TRDC+, TRGV9+, TRDV2+, CD8-, CD4-(Martin et al., Garcillan et al.)
)

pbmcs_T_marks <- c("CD4", "CD8A", "CD3E", "FOXP3", "CTLA4", "IL2RA", "IL7R", "KLF6", "FOS", "NCR3", "KLRB1", "TRAV1-2", "PRF1", "GZMA", "GZMB", "CCL4", "TRDC", "TRDV1", "TRDV2", "TRGV9")

# B subcompartment
pbmcs_B <- c(
"0" = "Naive",           # CD79B+, PAX5+, IGHD+, FCER2+ (Martin et al.)
"1" = "Plasmablast",     # GAPDH+, JCHAIN+, CD74+ (Martin et al.)
"2" = "Memory"           # MS4A1+, CD19+, IGHD- (Boland et al)
)

pbmcs_MONO <- c(
"0" = "Classical Monocytes",             # CD14+, S100A8/9+
"1" = "Nonclassical Monocytes",          # CD16+, CD14-
"2" = "Basophils",                       # FCER1A+
"3" = "Classical Monocytes"              # CD14+, S100A8/9+
)

pbmcs_MONO_marks <- c("FCER1A","S100A8","S100A9","CD14","FCGR3A")
