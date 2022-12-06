require(reshape2)
require(ggplot2)
require(Seurat)
require(plyr)
require(dplyr)
require(tidyr)
require(harmony)
require(gtools)
require(scales)
require(dittoSeq)
require(viridis)
require(purrr)
require(stringr)

### Common functions I use in scRNAseq data analysis ###

# customFindMarkers:
# SUMMARY      =   FUNCTION    = customFindMakers is used to perform pseudo-bulk analysis in parts
# inputs:
# SOBJ         = Seurat object = seurat object to perform analysis with
# PREFIX       =    STRING     = beginning of file name
# COMPARE      =    STRING     = which comparisons to run. By default, compares each condition to "HC". The string passed to COMPARE will be printed to the resulting image files.
# output:
# DATA.FRAME  = combined list of top genes for subsets of interest
customFindMarkers <- function(sobj,prefix,compare="HC"){
    Idents(sobj) <- sobj$condition
    conditions <- c("UCV","UCNB")
    GEX_list <- list()
    for (curr_cond in conditions){
        if (compare == "HC"){
            baseline <- "HC"
        } else {
            baseline <- conditions[is.na(pmatch(conditions,curr_cond))]
        }
        curr_marks <- FindMarkers(sobj, 
                                  group.by="condition", 
                                  ident.1=curr_cond, 
                                  ident.2=baseline, 
                                  test.use="poisson", 
                                  assay="RNA", 
                                  latent.vars="LIBRARY", 
                                  logfc.threshold=0.2, 
                                  min.pct=0.2, 
                                  only.pos=T)
        curr_marks$cluster <- paste0(curr_cond," versus ",baseline)
        curr_marks$gene <- rownames(curr_marks)
        rownames(curr_marks) <- NULL
        GEX_list[[paste0(curr_cond," versus ",baseline)]] <- curr_marks
    }
    merged <- reduce(GEX_list,rbind)
    top_marks <- unique((merged %>% group_by(cluster) %>% top_n(n=10,wt=avg_logFC))$gene)

    # Subset sobj to include relevant cells and annotate comparisons
    ssobj <- subset(sobj,ident=c("HC"),invert=T)
    if (baseline == "HC"){
        compared_conditions <- c(
            "UCV" = paste0("UCV versus ",baseline),
            "UCNB" = paste0("UCNB versus ",baseline)
        )
    } else {
        compared_conditions <- c(
           "UCV" = "UCV versus UCNB",
           "UCNB" = "UCNB versus UCV"
        )
    }
    ssobj@meta.data$condition_renamed <- compared_conditions[ssobj@meta.data$condition]

    pdf(paste(prefix,"condition_vs",compare,"dotplot.pdf",sep="_"),
              width=8,height=10,useDingbats=F)
    print(DotPlot(ssobj,group.by="condition_renamed",cols="RdYlBu",features=top_marks) +
                  coord_flip() +
                  theme(axis.text.x=element_text(angle=45,hjust=1)))
    dev.off()
    png(paste(prefix,"condition_vs",compare,"dotplot.png",sep="_"),
              width=8,height=10,units="in",res=150)
    print(DotPlot(ssobj,group.by="condition_renamed",cols="RdYlBu",features=top_marks) +
                  coord_flip() +
                  theme(axis.text.x=element_text(angle=45,hjust=1)))
    dev.off()
    return(merged)
}

# TCR_BCR_Metadata
# SUMMARY  = Computes expanded TCR/BCR counts
# inputs:
# SOBJ     = Seurat object = seurat object to perform analysis with
# DATA     =    STRING     = string specifying which of TCR or BCR data to use
# output:
# Seurat object with expanded TCR/BCR counts stored in metadata 
TCR_BCR_Metadata <- function(sobj,data="TCR"){
  duplicates <- paste0(data,"_duplicates")
  cdr_column <- paste0(data,".cdr3s_aa")
  counts <- paste0(data,"_counts")
  log_counts <- paste0(data,"_log_counts")
  sobj[[duplicates]] <- ifelse(is.na(sobj@meta.data[[cdr_column]]),paste0("no ",data), ifelse(duplicated(sobj@meta.data[[cdr_column]]) | duplicated(sobj@meta.data[[cdr_column]], fromLast=T),paste0("expanded ",data),data))
  counts_tbl <- table(sobj@meta.data[[cdr_column]])
  sobj[[counts]] <- as.numeric(counts_tbl[sobj@meta.data[[cdr_column]]])
  sobj[[log_counts]] <- log2(sobj@meta.data[[counts]] + 1)
  return(sobj)
}

# XAUT2_metadata
# SUMMARY    = adds in patient-level metadata for XAUT2 coarse data object
# inputs:
# SOBJ       = Seurat object = seurat object to add metadata to
# PREFIX     =    STRING     = beginning of file name
# ADD_META   =   BOOLEAN     = whether to add in the metadata from the external file specified below
# RM_PATIENT =   BOOLEAN     = whether to remove patient 16
# TCR_BCR    =   BOOLEAN     = whether to add in TCR/BCR metadata (see TCR_BCR_Metadata function)
# output:
# seurat object with updated metadata (if requested) and png/pdf images of each patient level plot needed for XAUT2 analysis
XAUT2_metadata <- function(sobj,prefix,add_meta=F,rm_patient=F,TCR_BCR=F){
  # Add in metadata
  if (add_meta){
    print("Adding patient-level metadata...")
    cluster_annotations <- read.table('final_clusters.tsv', sep='\t', header=T, row.names=1, stringsAsFactors=F)
    metadata <- read.table('/krummellab/data1/DSCoLab/XAUT2/10x/XAUT2_biopsy_pbmc_metadata.tsv',sep='\t',header=T)
    cluster_annotations <- setNames(cluster_annotations$assignment, rownames(cluster_annotations))
    names(cluster_annotations) <- paste0('P.', gsub(':CLUST', '.C', gsub('SCG', 'S', names(cluster_annotations))))
    assert_that(all(sort(names(cluster_annotations)) == sort(unique(sobj$SAMPLE.by.SNPs))))
    print("All patient annotations accounted for.")
    sobj$CoLabs_ID <- paste0(cluster_annotations[sobj$SAMPLE.by.SNPs],  gsub('.*POOL', '', sobj$orig.ident))
    sobj$CoLabs_patient <- gsub('XAUT2-|-BLD.*', '', sobj$CoLabs_ID)
    sobj$imx_pbmc_id <- gsub('-SCG.*',"",sobj$CoLabs_ID)
    all_meta <- join(sobj@meta.data,metadata)
    rownames(all_meta) <- rownames(sobj@meta.data)
    sobj@meta.data <- all_meta
    rm(all_meta)
    print("Done.")
  }
  # Remove HS16 (can affect TCR/BCR metadata, so this is run separately beforehand)
  if (rm_patient){
    print("Removing HS16...")
    Idents(sobj) <- sobj$CoLabs_patient
    sobj <- subset(sobj,idents=c("HS16"),invert=T)
    sobj <- ProcessSobj(sobj,prefix=prefix)
  }
  # Add in TCR/BCR metadata
  if (TCR_BCR){
    print("Computing TCR/BCR counts...")
    sobj <- TCR_BCR_Metadata(sobj)
    sobj <- TCR_BCR_Metadata(sobj,data="BCR")
  }
  print("Generating plots...")
  png(paste0(prefix,"_patient_umap.png"),width=8,height=6,units="in",res=100)
  print(dittoDimPlot(sobj,var="CoLabs_patient"))
  dev.off()
  png(paste0(prefix,"_disease_umap.png"),width=8,height=6,units="in",res=100)
  print(dittoDimPlot(sobj,var="disease"))
  dev.off()
  png(paste0(prefix,"_cancer_umap.png"),width=8,height=6,units="in",res=100)
  print(dittoDimPlot(sobj,var="cancer"))
  dev.off()
  png(paste0(prefix,"_CPI_umap.png"),width=8,height=6,units="in",res=100)
  print(dittoDimPlot(sobj,var="cpi"))
  dev.off()
  png(paste0(prefix,"_suppression_umap.png"),width=8,height=6,units="in",res=100)
  print(dittoDimPlot(sobj,var="suppression"))
  dev.off()
  png(paste0(prefix,"_TCR_umap.png"),width=16,height=6,units="in",res=100)
  print(dittoDimPlot(sobj,var="TCR_duplicates",split.by="TCR_duplicates")+NoLegend())
  dev.off()
  png(paste0(prefix,"_BCR_umap.png"),width=16,height=6,units="in",res=100)
  print(dittoDimPlot(sobj,var="BCR_duplicates",split.by="BCR_duplicates")+NoLegend())
  dev.off()
  png(paste0(prefix,"_log2_TCR_counts.png"),width=8,height=6,units="in",res=100)
  print(dittoDimPlot(sobj,var="TCR_log_counts"))
  dev.off()
  png(paste0(prefix,"_log2_BCR_counts.png"),width=8,height=6,units="in",res=100)
  print(dittoDimPlot(sobj,var="BCR_log_counts"))
  dev.off()
  png(paste0(prefix,"_library_split_umap.png"),width=12,height=8,units="in",res=100)
  print(dittoDimPlot(sobj,split.by="LIBRARY",var="coarse_annotations",split.ncol=4))
  dev.off()
  pdf(paste0(prefix,"_patient_umap.pdf"),width=8,height=6,useDingbats=FALSE)
  print(dittoDimPlot(sobj,var="CoLabs_patient"))
  dev.off()
  pdf(paste0(prefix,"_disease_umap.pdf"),width=8,height=6,useDingbats=FALSE)
  print(dittoDimPlot(sobj,var="disease"))
  dev.off()
  pdf(paste0(prefix,"_cancer_umap.pdf"),width=8,height=6,useDingbats=FALSE)
  print(dittoDimPlot(sobj,var="cancer"))
  dev.off()
  pdf(paste0(prefix,"_CPI_umap.pdf"),width=8,height=6,useDingbats=FALSE)
  print(dittoDimPlot(sobj,var="cpi"))
  dev.off()
  pdf(paste0(prefix,"_suppression_umap.pdf"),width=8,height=6,useDingbats=FALSE)
  print(dittoDimPlot(sobj,var="suppression"))
  dev.off()
  pdf(paste0(prefix,"_TCR_umap.pdf"),width=16,height=6,useDingbats=FALSE)
  print(dittoDimPlot(sobj,var="TCR_duplicates",split.by="TCR_duplicates")+NoLegend())
  dev.off()
  pdf(paste0(prefix,"_BCR_umap.pdf"),width=16,height=6,useDingbats=FALSE)
  print(dittoDimPlot(sobj,var="BCR_duplicates",split.by="BCR_duplicates")+NoLegend())
  dev.off()
  pdf(paste0(prefix,"_log2_TCR_counts.pdf"),width=8,height=6,useDingbats=FALSE)
  print(dittoDimPlot(sobj,var="TCR_log_counts"))
  dev.off()
  pdf(paste0(prefix,"_log2_BCR_counts.pdf"),width=8,height=6,useDingbats=FALSE)
  print(dittoDimPlot(sobj,var="BCR_log_counts"))
  dev.off()
  pdf(paste0(prefix,"_library_split_umap.pdf"),width=12,height=8,useDingbats=FALSE)
  print(dittoDimPlot(sobj,split.by="LIBRARY",var="coarse_annotations",split.ncol=4))
  dev.off()
  print("Done.")
}
# XAUT2_FreqDF
# SUMMARY      = computes cell type frequencies on a per-patient basis against a specified metadata column for analysis in XAUT2. Generalized form of XAUT1_FreqDF in the case where only one metadata column of interest is chosen.
# inputs:
# SOBJ         = Seurat object = seurat object to perform analysis with
# PREFIX       =    STRING     = beginning of file name if a file output is desired. By default, will not produce an output file.
# ANNOT_NAME   =    STRING     = metadata column name containing the cell type annotations
# META_COL     =    STRING     = metadata column name of interest for frequencies analysis in XAUT2
# output:
# A dataframe with the computed frequencies for the metadata column of interest normalized for each patient with respect to each of their cell type annotations. If a string is passed to PREFIX, then this dataframe is written to a csv file as well.
XAUT2_FreqDF <- function(sobj,prefix=NULL,annot_name,meta_col){
  sobj@meta.data$CoLabs_sample <- paste(sobj@meta.data[,"CoLabs_patient"],
                                  sobj@meta.data[,meta_col],
                                  sep=";")
  sep_cols <- c("patient",meta_col)
  vis_data <- table(sobj@meta.data[,c("CoLabs_sample",annot_name)])
  names(dimnames(vis_data)) <- c("patient","annotation")
  piv_table <- dcast(data.frame(vis_data),
                     patient~annotation,
                     value="Freq")
  rownames(piv_table) <- piv_table$patient
  piv_table$patient <- NULL
  frq_tbl <- piv_table/rowSums(piv_table)
  frq_data <- melt(as.matrix(frq_tbl), varnames=names(dimnames(vis_data)))
  frq_data <- frq_data %>% rename(frequency = value)
  frq_data <- separate(frq_data,col="patient",sep=";",into=sep_cols)
  if (is.character(prefix)){
    write.table(frq_data,file=paste0(prefix,"_",meta_col,"_frequencies.csv"),sep=",",quote=F,row.names=F)
  }
  return(frq_data)
}

# XAUT1_Freq_DF
# SUMMARY      = similar to XAUT2_FreqDF. Computes cell type frequencies on a per-patient basis against one or two specified metadata columns for XAUT1.
# inputs:
# SOBJ         = Seurat object = seurat object to perform analysis with
# PREFIX       =    STRING     = beginning of file name if a file output is desired. By default, will not produce an output file.
# ANNOT_NAME   =    STRING     = metadata column name containing the cell type annotations
# BIOPSY       =   BOOLEAN     = whether or not to separate the data by sample (i.e. side of the colon) instead of just by patient.
# output:
# A dataframe with the computed frequencies for the metadata column of interest normalized for each patient with respect to each of their cell type annotations. If a string is passed to PREFIX, then this dataframe is written to a csv file as well.
XAUT1_FreqDF <- function(sobj,prefix=NULL,annot_name,biopsy=TRUE){
  if(biopsy==TRUE){
    sobj@meta.data$CoLabs_sample <- paste(sobj@meta.data$CoLabs_patient,
                                          sobj@meta.data$colon_biopsy,
                                          sobj@meta.data$condition,
                                          sep="_")
    sep_cols <- c("patient","LR","condition")
    prefix = paste0(prefix,"_LR")
  } else {
    sobj@meta.data$CoLabs_sample <- paste(sobj@meta.data$CoLabs_patient,
                                    sobj@meta.data$condition,
                                    sep="_")
    sep_cols <- c("patient","condition")
  }
  vis_data <- table(sobj@meta.data[,c("CoLabs_sample",annot_name)])
  names(dimnames(vis_data)) <- c("patient","annotation")
  piv_table <- dcast(data.frame(vis_data),
                     patient~annotation,
                     value="Freq")
  rownames(piv_table) <- piv_table$patient
  piv_table$patient <- NULL
  frq_tbl <- piv_table/rowSums(piv_table)
  frq_data <- melt(as.matrix(frq_tbl), varnames=names(dimnames(vis_data)))
  frq_data <- frq_data %>% rename(frequency = value)
  frq_data <- separate(frq_data,col="patient",sep="_",into=sep_cols)
  if (is.character(prefix)){
    write.table(frq_data,file=paste0(prefix,"_frequencies.csv"),sep=",",quote=F,row.names=F)
  }
  return(frq_data)
}

GenePerPatient <- function(sobj,prefix,gene_list=NULL,subtype=NULL,cutoff=0.99){
  # -- Goal: create bar plots showing expression of gene per patient for a subset of our
  # data or a set of cell types
  # ARGUMENTS:
  # -- sobj is the seurat object of interest
  # -- prefix is the string we want all of the files to begin with
  # -- gene_list contains the name(s) of the gene(s) of interest. Handles up to two genes of interest
  # -- subtype is the cell phenotype of interest. Default is all cell types
  # -- cutoff is the minimum value of the log normalized expression of a gene of interest
  # that is considered "positive" for a given gene
  # RETURN this function returns the cell percentages broken down by patient for the specified cell subset

  ### Section 1: Identify cells of interest

  # Ensure there are no more than two genes in gene_list
  if (length(gene_list) > 2){
    stop(paste0("There are more than two genes in gene_list. Please only use one or two genes."))
  }

  genes <- try(FetchData(sobj,vars=gene_list))
  if (inherits(genes,'try-error')){
    stop(paste0("At least one of the genes of interest was not found. Ensure the gene names are correct and try again."))
  }
  # Identify compartment
  pattern <- "T|B|EPI|MYE|ENDO|MES|MONO"
  split <- unlist(strsplit(prefix,split="_"))
  compartment <- grep(pattern=pattern,x=split,value=TRUE)
  if (rlang::is_empty(compartment)){
    compartment <- "coarse"
    Idents(sobj) <- sobj$coarse_annotations
  } else {
    Idents(sobj) <- sobj$fine_annotations
  }
  
  # Subset down to cell type (if requested)
  if (is.character(subtype)){
    sobj <- subset(sobj,ident=subtype)
    subtype_fixed <- paste(unlist(strsplit(subtype,split="/")),collapse=" ")
    prefix <- paste(prefix,gsub(" ","_",subtype_fixed),sep="_")
    subtype <- paste0(subtype," ")
  } else {
    subtype <- ""
  }
  
  # Start creating gene expression breakdowns
  # If there are two genes...
  if (length(gene_list)==2){
    # Create double positive cell subset, if it exists
    expr <- FetchData(object = sobj, vars = gene_list)
    DP_cells <- try(sobj[,which((expr[[1]] > cutoff) & (expr[[2]] > cutoff))],silent=TRUE)
    # If there are no double positive cells, break
    if(inherits(DP_cells,'try-error')){
      stop(paste0("There are no ", gene_list[1],", ", gene_list[2]," double positive cells in the current subset: ", paste(compartment, subtype, sep = " ")))
    }
    # If there are double positive cells, create cell labels accordingly
    gene1_cells <- sobj[,which(expr[[1]] > cutoff)]
    gene2_cells <- sobj[,which(expr[[2]] > cutoff)]
    sobj$cells_of_interest <- ifelse(rownames(sobj@meta.data) %in% rownames(DP_cells@meta.data),
                                       paste(gene_list[1],gene_list[2],sep=" and "),
                                     ifelse(rownames(sobj@meta.data) %in% rownames(gene1_cells@meta.data),
                                              gene_list[1],
                                            ifelse(rownames(sobj@meta.data) %in% rownames(gene2_cells@meta.data),
                                                     gene_list[2],
                                                     paste0(gene_list[1],"-/",gene_list[2],"-"))))
  # If there is only one gene...
  } else {
    # Find the cells that express it
    expr <- FetchData(object = sobj, vars = gene_list)
    cells_of_interest <- try(sobj[,which(expr[[1]] > cutoff)],silent=TRUE)
    # If no cells are found, break
    if(inherits(cells_of_interest,'try-error')){
      stop(paste0("There are no ", gene_list[1]," positive cells in the ",compartment," subcompartment."))
    }
    # If there is at least one cell, create cell labels
    sobj$cells_of_interest <- ifelse(rownames(sobj@meta.data) %in% rownames(cells_of_interest@meta.data),
                                     paste0(gene_list[1]," positive"),
                                     paste0(gene_list[1]," negative"))
  }

  ### Section 2 : Compute percentages and plot
  
  # patient_short = CoLabs_patient - "XAUT1-"
  sobj$patient_short <- sapply(strsplit(sobj$CoLabs_patient,split="-"),`[`,2) 
  # Get cell counts for each patient and convert to dataframe format
  vis_data <- table(sobj@meta.data[,c("patient_short","cells_of_interest")])
  piv_table <- dcast(data.frame(vis_data),patient_short~cells_of_interest,value="Freq")
  rownames(piv_table) <- piv_table$patient_short
  piv_table$patient_short <- NULL
  # Convert counts to patient-normalized percentages and remove double negative cells
  prc_tbl <- (piv_table/rowSums(piv_table))*100
  prc_tbl <- prc_tbl %>% select(-contains("-/"))
  # Convert data to long format, adjust data types, then write to csv
  prc_data <- melt(as.matrix(prc_tbl), varnames=names(dimnames(vis_data)))
  prc_data <- prc_data %>% rename(percentage = value)
  prc_data$cells_of_interest <- factor(prc_data$cells_of_interest)
  prc_data$patient_short <- factor(prc_data$patient_short, levels=mixedsort(levels(prc_data$patient_short),decreasing=FALSE)) 
  prc_data <- prc_data[mixedorder(prc_data$patient_short),]
  write.table(prc_data,file=paste0(prefix,"_",paste(gene_list,collapse="_"),"_breakdown.csv"),
              sep=",",quote=FALSE,row.names=FALSE)
  # Begin plotting if needed
  # pdf(paste0(prefix,"_",paste(gene_list,collapse="_"),"_breakdown.pdf"),width=8,height=4,useDingbats=F)
  # print(ggplot(prc_data,aes(x=patient_short,y=percentage,fill=patient_short)) +
    # geom_bar(stat="identity") +
    # facet_wrap(~cells_of_interest) +
    # theme_dark() +
    # theme(axis.text.x=element_text(angle=45,hjust=1),
    #       plot.title = element_text(hjust = 0.5)) +
    # NoLegend() +
    # ggtitle(paste0("Percent ",subtype,"Cells Expressing ",paste(gene_list,collapse=" and ")," in ",compartment," Subcompartment")) +
    # scale_fill_viridis(discrete=T))
  # dev.off()  
  return(prc_data)
}

# CustomGenes
# SUMMARY      = converts results of DGE analysis to format for use in XAUT1 gene comparisons.
# inputs: 
# ANNOTATIONS  = NAMED VECTOR = vector containing cluster : annotation key value pairs
# GEX_FILE     =    STRING    = name of tsv file containing DGE analysis results
# output:
# dataframe in format to use for XAUT1 gene comparisons
CustomGenes <- function(annotations,GEX_file){
  GEX_markers <- read.table(GEX_file,header=T,sep="\t",stringsAsFactors=F)
  GEX_markers$cluster <- as.factor(GEX_markers$cluster)
  GEX_markers$annotations <- annotations[GEX_markers$cluster]
  cluster_gene <- subset(GEX_markers,select=c("annotations","gene"))
  reform <- as.data.frame(cluster_gene %>% 
      group_by(annotations) %>% 
      mutate(row = row_number()) %>% 
      tidyr::pivot_wider(names_from = annotations, values_from = gene) %>% 
      select(-row))
  return(reform)
}

# AnnotatedPlots
# SUMMARY     = Shorthand for creating annotated dotplot and umap
# inputs:
# SOBJ        = Seurat object = seurat object on which to perform analysis
# RES         =     FLOAT     = resolution at which to make plots (based on DGE analysis and annotations)
# PREFIX      =    STRING     = beginning of file name
# ANNOTATIONS = NAMED VECTOR  = vector containing cluster : annotation key value pairs
# markers     =    STRING     = name of tsv file containing DGE analysis results
# useDitto    =    BOOLEAN    = whether or not to use DittSeq to make color-blindness friendly plots
# output:
# png and pdf files of annotated dotplot and umap for the given seurat object. Saves the seurat object to .RData file if it is not present in the current working directory and returns the seurat object itself.
AnnotatedPlots <- function(sobj,res=0.2,prefix,annotations=NULL,markers=NULL,useDitto=TRUE){
    if (is.character(markers) & length(markers) == 1){
        GEX <- read.table(markers,sep="\t",header=TRUE,stringsAsFactors=FALSE)
        top_markers <- unique((GEX %>% group_by(cluster) %>% top_n(n=5,wt=avg_logFC))$gene)
    } else if (!is.data.frame(markers) & !is.null(markers)){
        top_markers <- unique(markers)
    }
    louvain_col <- paste0("louvain_res",res)
    pattern <- "merged_SCG\\d+_\\d+_(?:T|B|EPI|MYE|ENDO|MES|MONO)"
    if(grepl(pattern,prefix)){
        annot_col_name <- "fine_annotations"
    } else {
        annot_col_name <- "coarse_annotations"
    }
    if (! annot_col_name %in% colnames(sobj@meta.data)){
      sobj@meta.data[[annot_col_name]] <- annotations[sobj@meta.data[[louvain_col]]]
    }
    if (useDitto){
      pdf(paste0(prefix,"_ditto_res_",res,"_annotated_umap.pdf"),
          width=8,
          height=6,
          useDingbats=F)
      print(dittoDimPlot(object=sobj,var=annot_col_name,do.label=T,legend.show=F))
      dev.off()

      png(paste0(prefix,"_ditto_res_",res,"_annotated_umap.png"),
          width=8,
          height=6,
          units="in",
          res=150)
      print(dittoDimPlot(object=sobj,var=annot_col_name,do.label=T,legend.show=F))
      dev.off()
    } else {
      pdf(paste0(prefix,"_res_",res,"_annotated_umap.pdf"),
          width=8,
          height=6,
          useDingbats=F)
      print(DimPlot(sobj, group.by=annot_col_name))
      dev.off()

      png(paste0(prefix,"_res_",res,"_annotated_umap.png"),
          width=8,
          height=6,
          units="in",
          res=150)
      print(DimPlot(sobj, group.by=annot_col_name))
      dev.off()
    }
    pdf(paste0(prefix,"_res_",res,"_annotated_dotplot.pdf"),
        width=8,
        height=14,
        useDingbats=F)
    print(DotPlot(sobj,
                  group.by=annot_col_name,
                  cols="RdYlBu",
                  features=top_markers) +
          coord_flip() +
          theme(axis.text.x=element_text(angle=45,hjust=1)))
    dev.off()

    png(paste0(prefix,"_res_",res,"_annotated_dotplot.png"),
        width=8,
        height=14,
        units="in",
        res=150)
    print(DotPlot(sobj,
                  group.by=annot_col_name,
                  cols="RdYlBu",
                  features=top_markers) +
          coord_flip() +
          theme(axis.text.x=element_text(angle=45,hjust=1)))
    dev.off()

    if (! file.exists(paste0(prefix,"_annotated.RData"))){
      save(sobj,file=paste0(prefix,"_annotated.RData"))
    }
    return(sobj)
}

# loadRData
# SUMMARY  = Simple loader function for sobjs
# input:
# FILENAME = name of .RData file to pull Seurat object from
# output:
# The seurat object of interest
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

# RNAStats
# SUMMARY = Produce both nCount and nFeature RNA violin plots for a given sobj
# inputs:
# SOBJ    = Seurat object = seurat object to perform analysis with
# COL     =     FLOAT     = resolution to use in visualization
# PREFIX  =    STRING     = beginning of file name. If not specified, no plots are made
# output:
# returns a list of the statistics calculated for nCount_RNA and nFeature_RNA grouped by the specified resolution according to COL. If PREFIX is specified, violin plots are made for each.
RNAStats <- function(sobj,col=0.2,prefix="FALSE"){
    # Pass in numeric value for clusters, otherwise pass in full
    # metadata column name
    if (is.numeric(col) & prefix != "FALSE"){
      res_column <- paste0("louvain_res",col)
      prefix <- paste0(prefix,"_res_",col)
    } else {
      res_column <- col
      prefix <- paste0(prefix,"_",col)
    }
    # Create summary stats
    count_stats <- sobj@meta.data %>% 
        group_by(!!! rlang::syms(res_column)) %>% 
        summarise(sum=sum(nCount_RNA),
                  mean=mean(nCount_RNA),
                  sd=sd(nCount_RNA),
                  median=median(nCount_RNA),
                  min=min(nCount_RNA),
                  max=max(nCount_RNA))
    feature_stats <- sobj@meta.data %>% 
        group_by(!!! rlang::syms(res_column)) %>% 
        summarise(sum=sum(nFeature_RNA),
                  mean=mean(nFeature_RNA),
                  sd=sd(nFeature_RNA),
                  median=median(nFeature_RNA),
                  min=min(nFeature_RNA),
                  max=max(nFeature_RNA))
    # Create Violin plots for counts and features if requested
    if(prefix!="FALSE"){
        pdf(file=paste0(prefix,"_RNA_counts.pdf"),
            width=8,
            height=6,
            useDingbats=FALSE)
        print(VlnPlot(sobj,
                      group.by=res_column,
                      features="nCount_RNA",
                      pt.size=0))
        dev.off()
        pdf(file=paste0(prefix,"_RNA_features.pdf"),
            width=8,
            height=6,
            useDingbats=FALSE)
        print(VlnPlot(sobj,
                      group.by=res_column,
                      features="nFeature_RNA",
                      pt.size=0))
        dev.off()
    }
    return(list(count_stats,feature_stats))
}

saturate <- function(vec, sat=0, binary=FALSE){
  ###
  # DESCRIPTION
  # A Function to convert a vector of scores into a saturated vectore of scores. A saturated vector is one where all values below the
  # provided "saturation" (percentile of data) are set to 0. If the binary flag is specified, all values greater than or equal to the
  # saturation will be set to 1.
  #
  # INPUTS
  # vec: A numeric vector of scores
  # sat: A value to saturate the vector at (float (0.0-1.0) or percent (1.0-100.0))
  # binary: A flag to indicate if we should make the output vector a binary one.
  #
  # OUTPUT
  # A vector of saturated scores
  ###
  sat = if (sat > 1.0) sat/100 else sat
  z <- quantile(vec, sat)
  for (i in 1:length(vec)){
    if (vec[i] < z) {
      vec[i] = 0
    } else if(binary) {
      vec[i] = 1
    }
  }
  vec
}

# ProcessSobj
# SUMMARY   = merged R seurat object processing function
# inputs:
# SOBJ      = Seurat object = seurat object to perform analysis with
# PREFIX    =    STRING     = beginning of file name
# PLOT_GEX  =    BOOLEAN    = whether or not to perform DGE analysis on RNA data
# PLOT_CITE =    BOOLEAN    = whether or not to perform DGE analysis on protein data
# DEG_RES   =     FLOAT     = resolution at which to perform DGE analysis 
# output:
# png and pdf images of the umaps (resolutions from 0.2 to 2 in increments of 0.2), a tsv file summarizing the DGE analysis, and the seurat object
ProcessSobj <- function(sobj,prefix,plot_GEX=TRUE,plot_cite=FALSE,DEG_res=0.2){
    if (! file.exists(paste0(prefix, '.RData'))){
        sobj <- FindVariableFeatures(sobj, 
                                     selection.method = "vst", 
                                     nfeatures = 2000)
        sobj <- ScaleData(object = sobj, 
                          vars.to.regress = c("percent.mt",
                                              "nCount_RNA", 
                                              "percent.ribo", 
                                              "S.Score", 
                                              "G2M.Score"))
        sobj <- RunPCA(object = sobj)
        pdf(paste0(prefix, '_harmony_convergence.pdf'))
        sobj <- RunHarmony(sobj, 
                           "LIBRARY", 
                           plot_convergence = TRUE, 
                           max.iter.harmony=30, 
                           max.iter.cluster=30)
        dev.off()
        
        # Run clustering analysis
        sobj <- FindNeighbors(sobj,
                              reduction="harmony",
                              verbose=TRUE)

        # run umap with reduction="harmony"
        sobj <- RunUMAP(sobj,
                        dims = 1:30,
                        n.neighbors = 30,
                        min.dist = 0.3,
                        spread = 1,
                        verbose = FALSE,
                        reduction="harmony")

        # Plot umap at various resolutions
        for (res in seq(from = 0.2, to = 2, by = 0.2)){
          sobj <- FindClusters(sobj, verbose = TRUE,
                               algorithm = 1,
                               resolution = res)
          sobj@meta.data[[paste0('louvain_res', res)]] <- Idents(sobj)
          
          png(filename=paste0(prefix, '_res_', res, '_umap.png'), width = 8, height = 6, units = "in", res = 100)
          print(DimPlot(sobj, group.by=paste0('louvain_res', res), label=T) + NoLegend())
          dev.off()
        }
        save(sobj,file=paste0(prefix,".RData"))
    }
    # Compute and plot DEGs per cluster (GEX) if requested
    if (plot_GEX == TRUE){
        Idents(sobj) <- sobj@meta.data[[paste0("louvain_res",DEG_res)]]
        GEX_markers <- FindAllMarkers(sobj,
                                      test.use='poisson',
                                      latent.vars=c('LIBRARY','LIBRARY.TYPE'),
                                      assay='RNA',
                                      logfc.threshold=0.2,
                                      min.pct=0.2,
                                      only.pos=TRUE)
        write.table(GEX_markers,
                    file=paste0(prefix,"_GEX_markers.tsv"),
                    quote=FALSE,
                    sep="\t",
                    row.names=FALSE)
        top_markers <- unique((GEX_markers %>% 
                               group_by(cluster) %>% 
                               top_n(n=5,wt=avg_logFC))$gene)

        pdf(paste0(prefix, "_res_", DEG_res, "_dotplot.pdf"),
            width=8,
            height=12,
            useDingbats=FALSE)
        print(DotPlot(sobj,
                      features=top_markers,
                      cols="RdYlBu") +
              coord_flip())
        dev.off()
    }
    # If requested, do DEG analysis for CITEseq
    if (plot_cite == TRUE){
      col <- paste0("louvain_res",DEG_res)
      sobj <- plot_CITE(sobj,prefix,meta_col=col)
    }
    sobj
}

# plot_CITE
# SUMMARY   = computes and plots results from DGE analysis on protein data
# inputs:
# SOBJ      = Seurat object = seurat object to perform analysis with
# PREFIX    =    STRING     = beginning of file name
# META_COL  =    STRING     = name of metadata column to perform DGE analysis against
# output:
# png and pdf images of the dotplot summarizing DGE analysis, a tsv file containing the DGE analysis results, and the seurat object
plot_CITE <- function(sobj,prefix,meta_col){
  Idents(sobj) <- sobj@meta.data[[meta_col]]
  CITE_markers <- FindAllMarkers(sobj,
                                 assay='ADT',
                                 logfc.threshold=0.2,
                                 min.pct=0.2,
                                 only.pos=TRUE)
  write.table(CITE_markers,
              file=paste0(prefix,"_CITE_markers.tsv"),
              quote=FALSE,
              sep="\t",
              row.names=FALSE)
  top_markers <- unique((CITE_markers %>% 
                         group_by(cluster) %>% 
                         top_n(n=5,wt=avg_logFC))$gene)
  pdf(paste0(prefix,"_CITE_dotplot.pdf"),
      width=8,
      height=12,
      useDingbats=FALSE)
  if (max(nchar(levels(Idents(sobj)))) > 2){
    print(DotPlot(sobj,
                  assay="ADT",
                  features=top_markers,
                  cols="RdYlBu") +
          coord_flip() +
          theme(axis.text.x=element_text(angle=45,hjust=1)))
  } else {
    print(DotPlot(sobj,
                  assay="ADT",
                  features=top_markers,
                  cols="RdYlBu") +
          coord_flip())
  }
  dev.off()
  png(paste0(prefix,"_CITE_dotplot.png"),
      width=8,
      height=12,
      units="in",
      res=150)
  if (max(nchar(levels(Idents(sobj)))) > 2){
    print(DotPlot(sobj,
                  assay="ADT",
                  features=top_markers,
                  cols="RdYlBu") +
          coord_flip() +
          theme(axis.text.x=element_text(angle=45,hjust=1)))
  } else {
    print(DotPlot(sobj,
                  assay="ADT",
                  features=top_markers,
                  cols="RdYlBu") +
          coord_flip())
  }
  dev.off()
  return(sobj)
}

