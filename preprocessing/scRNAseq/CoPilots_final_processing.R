#########################################################################
# 
# Author: Chris Andersen
# GOAL: Creating all the final files for CoPilots (6/18/2021)
# DESCRIPTION: The following script can be used to perform a variety of
# analyses performed for IBD CoPilots XAUT1 and XAUT2. The input data
# is assumed to be a fully annotated set of Seurat objects (.RData format)
# that only need to be analyzed. The details of many of the analyses
# are thoroughly described in the assocated "RNA_seq_functions.R" script
# that accompanies this one. See that script for the list of required
# packages and be sure to install each of them before running this script.
# This script is generally set up to perform analysis for the
# subcompartments (e.g. T cells, B cells, etc.) for each data set, but
# not the coarse data. This script performs the following analyses:
# - Compute relative cell frequencies per annotation
# - Perform DGE analysis on specified metadata columns (see code below)
# - Perform DGE analysis on CITE-Seq data only
# - Compute relative abundances of cells expressing ITGA4 and ITGB7 (or
#   MADCAM1)
# - Add in metadata from external file
# - Generate dotplot and UMAP visualizations (only applicable for certain
#   subsets as written)
# 
#########################################################################

# myFindAllMarkers
# SUMMARY  = Shorthand function for performing DGE analysis for a given metadata column
# inputs:
# SOBJ     = Seurat object = seurat object on which to perform DGE analysis
# PREFIX   =     STRING    = name to begin each file with
# SUFFIX   =     STRING    = name to end each file with (before the file extension)
# METAVAR  =     STRING    = name of metadata column which DGE analysis is to to be performed against
# output:
# png and pdf files of dotplot, a tsv file containing detailed DGE results, and the seurat object
myFindAllMarkers <- function(sobj,prefix,suffix,metavar){
    Idents(sobj) <- sobj[[metavar]]
    curr_GEX <- FindAllMarkers(sobj,test.use='poisson',latent.vars=c('LIBRARY'),assay='RNA',logfc.threshold=0.2,min.pct=0.2,only.pos=TRUE)
    write.table(curr_GEX,file=paste0(prefix,"_",suffix,".tsv"),sep="\t",quote=F,row.names=F)
    top_GEX <- unique((curr_GEX %>% group_by(cluster) %>% top_n(n=10,wt=avg_logFC))$gene)
    png(paste0(prefix,"_",suffix,".png"),width=8,height=8,units="in",res=150)
    print(DotPlot(sobj,features=top_GEX,cols="RdYlBu")+coord_flip()+theme(axis.text.x=element_text(angle=45,hjust=1)))
    dev.off()
    pdf(paste0(prefix,"_",suffix,".pdf"),width=8,height=8,useDingbats=FALSE)
    print(DotPlot(sobj,features=top_GEX,cols="RdYlBu")+coord_flip()+theme(axis.text.x=element_text(angle=45,hjust=1)))
    dev.off()
    return(sobj)
}

# Pull in my custom tool set for RNAseq
source("/krummellab/data1/chrisquatjr/candersen_scripts/RNA_seq_functions.R")

# Specify which dataset to use and which analysis to perform
project <- "XAUT2"
do_dot_umap <- FALSE
update_metadata <- FALSE
do_freq <- FALSE
do_A4B7 <- FALSE
do_pseudobulk <- FALSE
do_CITE <- TRUE
if (project == "XAUT1"){
    datasets <- c("merged_SCG1_10","merged_SCG11_14")
    prefix <- "FINAL"
} else {
    datasets <- c("merged_SCG1_8","merged_SCG9_12")
    prefix <- "Analyzed"
}
# Pull in all annotations relevant to this project
source(paste0("/krummellab/data1/chrisquatjr/candersen_scripts/Annotation_notes_",project,".R"))
# For biopsies then pbmcs...
for (dataset in datasets){
    if (dataset %in% c("merged_SCG1_10","merged_SCG1_8")){
        data_name <- "biopsies"
    } else {
        data_name <- "pbmcs"
    }
    datadir <- paste0("/krummellab/data1/DSCoLab/",project,"/10x/",dataset,"/subclusters/")
    setwd(datadir)
    print(paste0("Starting analysis for ",data_name," in project ",project,"."))
    for (subset in Sys.glob("*")){
        
        # Navigate to each subcompartment, identify its code (e.g. T, B, EPI, etc.)
        setwd(paste0(datadir,"/",subset))
        if (subset != "Monocytes"){
            subset_name <- toupper(unlist(strsplit(subset,split="_"))[1])[1]
        } else {
            subset_name <- "MONO"
        }
        # print which cell type is being processed, define the appropriate names
        print(paste0("Processing ",subset_name," cells..."))
        file_prefix <- paste(prefix,dataset,subset_name,sep="_")
        object_name <- paste(file_prefix,"annotated.RData",sep="_")

        # If the expected RData object is found, load it, otherwise move to next subset
        if (file.exists(object_name)){
            sobj <- loadRData(object_name)
        } else {
            print(paste0("No object found for ",subset_name," cells in ",dataset,", skipping."))
            next
        }

        # Make sure fine annotations column is a factor
        if (! is.factor(sobj@meta.data$fine_annotations)){
            annots <- paste0(data_name,"_",subset_name)
            sobj@meta.data$fine_annotations <- factor(sobj@meta.data$fine_annotations,levels=unique(unname(get(annots))))
        }

        if (do_CITE){
            if (dataset == "merged_SCG11_14"){
                next
            } else {
                tryCatch({
                    sobj <- plot_CITE(sobj,prefix=file_prefix,meta_col="fine_annotations")
                },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            }
        }

        # If prompted, perform pseudo-bulk analysis
        if (do_pseudobulk){
            if (project == "XAUT1"){
                # Pseudo-bulk analyses for XAUT1
                if (! "condition" %in% colnames(sobj@meta.data)){
                    print(paste0("Column \"condition\" not found in ",subset_name," metadata. Skipping."))
                    next
                }
                pseudo_HC <- customFindMarkers(sobj,prefix=file_prefix)
                pseudo_COND <- customFindMarkers(sobj,prefix=file_prefix,compare="condition")
                pseudo_marks <- rbind(pseudo_HC,pseudo_COND)
                write.table(pseudo_marks,file=paste0(file_prefix,"_pseudo_bulk_markers.tsv"),sep="\t",quote=F,row.names=F)
            } else {
                # Pseudo-bulk analyses for XAUT2
                # Each versus complement of each metavariable
                setwd(paste0(datadir,"/",subset,"/patient_level/pseudo_bulk/"))
                for (curr_meta in c("disease","cpi","suppression")){
                    sobj <- myFindAllMarkers(sobj,prefix=file_prefix,suffix=paste0(curr_meta,"_all"),metavar=curr_meta)
                    if (curr_meta == "cpi"){
                        Idents(sobj) <- sobj[[curr_meta]]
                        cpi_sobj <- subset(sobj,ident=c("none"),invert=TRUE)
                        cpi_sobj <- myFindAllMarkers(cpi_sobj,prefix=file_prefix,suffix="cpi_pd-1_vs_combo",metavar=curr_meta)
                        rm(cpi_sobj)
                        sobj@meta.data$cpi_treated <- ifelse(sobj@meta.data$cpi == "none","none","treated")
                        sobj <- myFindAllMarkers(sobj,prefix=file_prefix,suffix="cpi_treated_vs_none",metavar="cpi_treated")
                    }
                    if (curr_meta == "suppression"){
                        sobj@meta.data$steroid_ifx_none <- ifelse(sobj@meta.data$suppression == "none" | sobj@meta.data$suppression == "steroid",as.character(sobj@meta.data$suppression),"ifx_and_or_ada")
                        sobj <- myFindAllMarkers(sobj,prefix=file_prefix,suffix="suppression_ifx_grouped",metavar="steroid_ifx_none")
                        sobj@meta.data$steroid_vs_all <- ifelse(sobj@meta.data$suppression == "steroid","steroid","no steroid")
                        sobj <- myFindAllMarkers(sobj,prefix=file_prefix,suffix="steroid_vs_all",metavar="steroid_vs_all")
                    }
                }
            }
        }

        # If prompted, update the metadata and save it to a separate file
        if (update_metadata){
            setwd(paste0(datadir,"/",subset))
            if (project == "XAUT2"){
                updated_meta <- read.table("/krummellab/data1/DSCoLab/XAUT2/10x/XAUT2_biopsy_pbmc_metadata_update.tsv",sep="\t",header=T,row.names=NULL)
                sobj@meta.data$cpi <- NULL; sobj@meta.data$disease <- NULL
                if (data_name == "biopsies"){
                    all_meta <- join(sobj@meta.data,updated_meta,by="imx_tissue_id")
                    rownames(all_meta) <- rownames(sobj@meta.data)
                    all_meta$imx_pbmc_id <- NULL; sobj@meta.data$imx_pbmc_id <- NULL
                } else {
                    all_meta <- join(sobj@meta.data,updated_meta,by="imx_pbmc_id")
                    rownames(all_meta) <- rownames(sobj@meta.data)
                    all_meta$imx_tissue_id <- NULL; sobj@meta.data$imx_tissue_id <- NULL
                }
                sobj@meta.data <- all_meta
                write.table(sobj@meta.data,file=paste0(file_prefix,"_metadata.tsv"),sep="\t",quote=FALSE)
                save(sobj,file=paste0("NEW_",file_prefix,"_annotated.RData"))
            } else {
                # For XAUT1, just save metadata to tsv file
                write.table(sobj@meta.data,file=paste0(file_prefix,"_metadata.tsv"),sep="\t",quote=FALSE)
            }
        }

        # If prompted, compute cell frequencies using appropriate project-specific functions
        if (do_freq){
            setwd(paste0(datadir,"/",subset,"/patient_level/cell_fractions/"))
            if (project == "XAUT1"){
                if (data_name == "biopsies"){
                    LR_freq <- XAUT1_FreqDF(sobj,prefix=file_prefix,annot_name="fine_annotations")
                    freq <- XAUT1_FreqDF(sobj,prefix=file_prefix,annot_name="fine_annotations",biopsy=FALSE)
                } else {
                    freq <- XAUT1_FreqDF(sobj,prefix=file_prefix,annot_name="fine_annotations",biopsy=FALSE)
                }
            } else {
                for (var_of_int in c("cpi","disease","suppression")){
                    freq <- XAUT2_FreqDF(sobj,prefix=file_prefix,annot_name="fine_annotations",meta_col=var_of_int)
                }
            }
        }

        # If prompted, compute A4B7 frequencies per subtype (assumes but does not check that project == XAUT1)
        if(do_A4B7 == TRUE){
            setwd(paste0(datadir,"/",subset,"/patient_level/A4B7/"))
            # Save frequencies for A4B7 across compartment and within each fine annotation
            tryCatch({
                freq <- GenePerPatient(sobj,prefix=file_prefix,gene_list=c("ITGA4","ITGB7"))
                for (subAnnot in unique(levels(sobj@meta.data$fine_annotations))){
                    tryCatch({
                        freq <- GenePerPatient(sobj,prefix=file_prefix,subtype=subAnnot,gene_list=c("ITGA4","ITGB7"))
                    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                }
            },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                # For endothelial cells, also compute frequencies for MADCAM1
            if (subset_name == "ENDO"){
                tryCatch({
                    freq <- GenePerPatient(sobj,prefix=file_prefix,gene_list=c("MADCAM1"))
                    for (subAnnot in unique(levels(sobj@meta.data$fine_annotations))){
                        tryCatch({
                            freq <- GenePerPatient(sobj,prefix=file_prefix,subtype=subAnnot,gene_list=c("MADCAM1"))
                        },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                    }
                },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            }
            setwd("../")
        }

        # If prompted, create annotated dotplot and UMAP
        if (do_dot_umap){
            # T cells handled separately
            if (subset_name %in% c("T","MONO")){
                print("Dotplot and umap handled separately for this subset. Skipping.")
                next
            }
            setwd(paste0(datadir,"/",subset))
            annots <- paste0(data_name,"_",subset_name)
            marks <- paste0(file_prefix,"_GEX_markers.tsv")
            Idents(sobj) <- sobj$fine_annotations
            GEX_marks <- FindAllMarkers(sobj,test.use='poisson',latent.vars=c('LIBRARY'),assay='RNA',logfc.threshold=0.2,min.pct=0.2,only.pos=TRUE)
            write.table(GEX_marks,file=marks,sep='\t',quote=F)
            sobj <- AnnotatedPlots(sobj,prefix=file_prefix,annotations=get(annots),markers=marks)
            print("Done.")
        }
    }
}
