# Differential expression using MAST 
# Input: Annotated Seurat object 
#  


suppressPackageStartupMessages({
  library(dplyr)
  library(zinbwave)
  library(SingleCellExperiment)
  library(Matrix)
  library(Seurat)
  library(MAST)
  library(doParallel)
  library(BiocParallel)
  
})

# # define the directory for saving the result of this DGE test (should be an input for this function)
# DGE_DIR = '/mnt/ibm_lg/yangjoon.kim/UC_UCSF_Multiome/MAST/result/'

# threshold for the minimum fold-change
FCTHRESHOLD = log2(1.5)

# output 
# output_dir= '/mnt/ibm_lg/yangjoon.kim/UC_UCSF_Multiome/MAST/'
  
# raw input files 
NCORES = 50
#registerDoParallel(NCORES)
#register(DoparParam())

options(mc.cores = NCORES)



####### Ulcerative Colitis project specific parts #########
# input variables required to run MAST:
# 1) input, output directories
# 2) dataset names, condition, baseline_conds, contrast_conds
# datast names are like pairs of conditions(disease_status), such as "subset_condition_XAUT1_Blood_cond1_cond2"
# condition = name of the category that we'd like to compare with (i.e. conditions, disease_status, etc.)
# baseline_conds = a list of categories that are baselines for DGE (i.e. cond1)
# contrast_conds = a list of categories that we calculate the DGE from the baselines (i.e. cond2)
# 3) options = (1) whole_cells=TRUE: do not subset for each cell type, just get the DGE for all cells between "conditions".


# define the dataset names, input and output directories,  at the very top
# 1) define the input and output directories
output_dir = '/mnt/ibm_lg/yangjoon.kim/UC_UCSF_Multiome/MAST/subset_disease_XAUT2_Blood/' # write DGE .csvs here 
input_dir_prefix= '/mnt/ibm_lg/yangjoon.kim/UC_UCSF_Multiome/MAST/subset_disease_XAUT2_Blood/'
# 2) dataset names (individual h5ad files)
datasets = c('XAUT2_Blood_hc_cpi_colitis', 'XAUT2_Blood_hc_uc', 'XAUT2_Blood_cpi_colitis_uc')
# 3) conditions, baseline_conds, contrast_conds
condition = 'condition' # or "disease_status"
list_baseline_conds = c('hc','hc','uc')
list_contrast_conds = c('cpi_colitis','uc','cpi_colitis')
# 4) annotation_granularity (name of the cell-type annotation group. i.e. coarse or fine)
annotation = "fine"


# Run for all datasets (multiple h5ad objects which have a pair of "conditions/disease_status")
dgeAllDatasets <- function(input_dir_prefix = 'input_dir_prefix', output_dir = 'output_dir', datasets='datasets', # directories
                   condition = 'condition',list_baseline_conds = c('1'), list_contrast_conds = c('1'), # conditions
                   regr_ngenes = T, out_id=annotation, ann_col=annotation){
  
  for(t in 1:length(datasets)){
    # define the dataset
    dataset = datasets[t]
    # define the conditions (disease_states)
    condition = condition
    baseline_cond = list_baseline_conds[t]
    contrast_cond = list_contrast_conds[t]
    #contrast_cond = paste0(category,cond_name)
    
    dgeDataset(input_dir_prefix = input_dir_prefix,
             output_dir = output_dir,
             dataset =dataset, 
             ann_col = ann_col, 
             data_id = out_id, 
             regress_ngenes = regr_ngenes,
             condition = condition,
             baseline_cond = baseline_cond,
             contrast_cond = contrast_cond)
    print(paste0('DONE tissue ', datasets[t]))
    gc() 
  }
}

# Run all cell types for a given dataset 
dgeDataset <- function(input_dir_prefix = 'input_dir_prefix',
                       output_dir = 'output_dir',
                       dataset ='dataset', 
                       ann_col = 'fine', 
                       data_id = 'fine', 
                       regress_ngenes = T,
                       condition = 'condition',
                       baseline_cond = c('1'),
                       contrast_cond = c('2'),
                       AllCells=FALSE){
  # 1. read from organ directory 
  input_dir= paste0(input_dir_prefix, dataset, "/r_files/")
  sce <- load_sce_object(input_dir, organ = dataset )
  # 2. list all cell types for this organ
  colData(sce)[ann_col][,1] %>% unique() -> cell_type_list



  
  all_dge = data.frame() 
  
  # Regress number of detected genes 
  if(regress_ngenes){
    cov_var = " + cgeneson"
  }else{ 
    cov_var = ""
  }
  
    for(c in 1:length(cell_type_list)){
      
      sce_sub <- splitAndProcess(sce, ann_labels = ann_col, cell_type = cell_type_list[c], min_counts = 1, min_cells = 5)
      # check the dimension of the sce_sub - just to make sure whether we selected for HVGs or the whole transcriptome
      print(dim(sce_sub))
      print(cell_type_list[c])
      #browser() 
      n_factors <- colData(sce_sub)["condition"][,1] %>% table %>% length
      # run differential expression ONLY if there are 2 factors for this cell type
      if(n_factors>1){
        # NOTE: by default we correct for number of detected genes 
        
        fc_results <- dge_MAST(sce_sub, covariate = cov_var, condition, baseline_cond, contrast_cond) 
        # 4. Write files for each cell type
        # collapse white spaces 
        #celltype_file = paste(str_split(cell_type_list[c], " ")[[1]], collapse ="_")
        # remove special characters 
        #celltype_file =  str_replace(celltype_file,'/', '_')
        fc_results$cell_type  = cell_type_list[c]
        
        all_dge = rbind(all_dge, fc_results)
        print(dim(fc_results))
      }
    }
    write.csv(all_dge, file = paste0(output_dir, dataset,"_DGE_allcelltypes_" , data_id, ".csv") , quote = F, row.names = T)
  # }
  
}

# Run MAST for specific pairs of conditions for one cell-type (group)
dge_MAST <- function(core = c(), covariate = "", condition = 'condition',
                     baseline_cond=c('1'),
                     contrast_cond=c('2')){
  # Export count matrix
  counts = assay(core)
  
  # MAST requires log data
  tpm <- counts*1e6/colSums(counts)
  tpm <- log2(tpm+1)
  
  core_log = core 
  assay(core_log) <- tpm
  
  # MAST sca object structure
  sca <- SceToSingleCellAssay(core_log, class = "SingleCellAssay", check_sanity = FALSE)
  
  # define the baseline and the experimental group
  condition = condition
  baseline_cond = baseline_cond
  contrast_cond = contrast_cond
  
  #browser()
  cond<-factor(colData(sca)$condition) # 'condition' should be an input at the very top level
  cond<-relevel(cond, baseline_cond) # "baseline_cond" should also be an input at the very top level
  colData(sca)$condition<-cond
  
  contrast_cond = paste0(condition,contrast_cond)
  
  #zlmCond <- zlm(formula = as.formula(paste0("~condition", covariate)), sca=sca)
  zlmCond<-zlm(~condition + cgeneson, sca)
  # Change the doLRT name for each condition pair (this needs to be edited in future)
  summaryCond <- summary(zlmCond, doLRT=contrast_cond)
  
  summaryDt <- summaryCond$datatable
  # dt1 = summaryDt[contrast=="disease_statusCov19" & component=="H", .(primerid, `Pr(>Chisq)`)]
  # dt2 = summaryDt[contrast=="disease_statusCov19" & component=="logFC", .(primerid, coef, z)]
  # de_res = merge(dt1, dt2, by="primerid")
  # colnames(de_res) <- c("gene", "age.H_p", "age.logFC", 'age.logFC_z')
  # de_res$age.H_fdr <- p.adjust(de_res$age.H_p, "BH")
  
  # Here, we again hard-coded the constrast group (the Experimental group compared to the Control group)
  fcHurdle <- merge(summaryDt[contrast==contrast_cond & component=='H',.(primerid, `Pr(>Chisq)`)],
                    summaryDt[contrast==contrast_cond & component=='logFC', .(primerid, coef, ci.hi, ci.lo, z)], by='primerid')
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  
  fcHurdleSig <-fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD]
  #fcHurdleSig <-fcHurdle[fdr<=1.0 & abs(coef)>=0]
  # fdr <=1 and FCTHRESHOLD>=0
  
  return(fcHurdleSig)
}



# Reads exported files from scanpy and creates a SCE object 
# counts are assumed to be raw
# Since QC already took place in scanpy, here we only focus on creating the object 
load_sce_object <- function( file_path = "",
                             organ = ""){
  #### 1. Prepare data 
  # Assumes general data structure exported from scanpy 
  # read from .mtx (better)
  counts = readMM(paste0(file_path,'sparse_matrix.mtx'))
  # obs data.frame
  obs = read.csv( paste0(file_path, 'obs.csv'), header = T)
  # var data.frame
  vars = read.csv( paste0(file_path, 'var.csv'), row.names = 1)
  # obsm UMAP coords
  obsm = read.csv( paste0(file_path,'obsm.csv'), header = T)
  
  #check that vars is not empty
#  if(dim(vars)[2]==0){
#    vars$gene_symbol = row.names(vars)
#    #vars$gene_ids = row.names(vars)
#    #vars$feature_types = 'Gene Expression'
#  }
  
  # define the gene_symbol using the first column
  vars$gene_symbol = row.names(vars)
  
  # for CTA datasets, we replace the name of the virus transcripts
  # fix the virus name
  #vars$gene_symbol[grep('SARS',vars$gene_ids)] <- vars$gene_ids[grep('SARS',vars$gene_ids)]
  #fix the exmpy first element in vars (some bug from scanpy after concatenation)
  if(vars$gene_symbol[1]=="") 
    vars$gene_symbol[1] = "101"
  
  colnames(counts) <- vars$gene_symbol
  row.names(counts) <- obs[,1]
  
  row.names(obs) <- obs[,1]
  row.names(obsm) <- obs[,1]
  
  # Create SCE object 
  # for some reason is not default in readMM
  c <- as(counts, 'dgCMatrix') %>% as.matrix
  sce <- SingleCellExperiment(list(counts=t(c)),
                              colData=obs,
                              rowData=vars,
                              metadata=list(tissue= organ)
  )
  
  return(sce)
  
}

# From the organ object, subset for a specific cell type 
# Performs normalization and filtering of lowly expressed genes 
splitAndProcess <- function(sce =c() , ann_labels = 'cell_type_annotation', 
                            cell_type = 'hepatocyte',
                            min_counts = 2, # filter lowly expressed genes
                            min_cells = 5, 
                            n_top_genes = 0 # default, consider all genes in the object for DGE 
){
  
  # 2. Set up metadata for differntial gene expression
  
  cell_idx <- colData(sce)[ann_labels] %in% c(cell_type) %>% unlist()
  core <- sce[,cell_idx]
  
  # relaxed gene filter: 
  core = core[rowSums(assay(core)) > 0, ]
  # stronger filter: only genes expressed 
  # ideally we want to use highly variable genes. 
  # but we don't want to compute DGE from genes expressed in a few cells. 
  filter <- rowSums(assay(core)>min_counts)>min_cells
  core <- core[filter,]
  
  # further filter by high variance 
  assay(core) %>% log1p %>% rowVars -> vars
  names(vars) <- rownames(core)
  
  # Number of genes expressed (to use as covariate)
  cdr2 <-colSums(assay(core)>0)
  colData(core)$cgeneson <- scale(cdr2)
  
  # top variable genes -- quick and dirty tbh
  vars <- sort(vars, decreasing = TRUE)
  if(n_top_genes){
    core <- core[names(vars)[1:n_top_genes], ]
  }
  
  return(core)
}

# merge all DE results into a single data.frame 
# this function also cleans the cell type annotations 
merge_organs <- function(tissues = c('liver','lung','kidney','prostate','heart'), 
                         file_id = 'DGE_allcelltypes_ngenes_regression.csv', 
                         output_file_id = '', covariate = 'ngenes',
                         file_name = 'file_name', input_dir = 'input_dir',
                         write_out = TRUE){ 
  
  # Read DE output for MAST and merge into a single data.frame
  organs_dge = list()
  for(t in 1:length(tissues)){
    organs_dge[[t]] = read.csv(file =paste0(input_dir, tissues[t], '_', file_id), row.names = 1, header = T)
    organs_dge[[t]]$tissue = tissues[t]
  }
  master_df <- do.call(rbind, organs_dge )
  
  master_df %>% mutate(full_cell_type = paste0(tissue, ": " , cell_type)) -> master_df
  master_df %>% dplyr::rename(gene = primerid, pval = fdr, log2fc = coef) -> master_df 
  master_df$method = 'MAST'
  master_df$covariate = covariate 
  master_df %>% dplyr::select(gene, pval, log2fc, cell_type, tissue, method, covariate) -> master_df
  
  if(write_out){
    write.csv(master_df, file = paste0(output_dir, file_name, 'all_celltypes_DGE_',output_file_id,'_MAST.csv'), quote = F, row.names = F)
  }else{
    return(master_df)
  }
}

find_shared_substring <- function(words) {
  #extract substrings from length 1 to length of shortest word
  subs <- sapply(seq_len(min(nchar(words))), 
                 function(x, words) substring(words, 1, x), 
                 words=words)
  #max length for which substrings are equal
  neqal <- max(cumsum(apply(subs, 2, function(x) length(unique(x)) == 1L)))
  #return substring
  substring(words[1], 1, neqal)
}