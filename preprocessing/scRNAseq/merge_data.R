library(cowplot)
library(dplyr)
library(ggplot2)
library(harmony)
library(RColorBrewer)
library(reshape2)
library(Seurat)
library(yaml)

#options(future.globals.maxSize = 200 * 1024^3)  # 200Gb

setwd("/krummellab/data1/DSCoLab/XAUT1/10x/merged_SCG1_10")
prefix='merged_SCG1_10'
if (! file.exists(paste0(prefix, '_merged_harmonized.RData'))){
  if (! file.exists(paste0(prefix, '_merged_temp.RData'))){
    setwd("/krummellab/data1/immunox/XAUT1/data/single_cell_GEX/processed")

    sobjs <- list()
    for (sobj in Sys.glob('*/automated_processing/*scTransformed.RData')){
      sample_name <- dirname(dirname(sobj))
      cat(paste0('processing sample: ', sample_name, "\n"))
      load(sobj)
      sobjs[[sample_name]] <- get(sample_name)
      rm(list=sample_name)
      temp <- read.table(file=gsub('.RData', '_processed_DF.tsv', sobj),
                         sep='\t',
                         header=T,
                         row.names=1,
                         stringsAsFactors=F)
      temp <- temp[, 'DF.DROPLET.TYPE', drop=F]
      sobjs[[sample_name]] <- AddMetaData(sobjs[[sample_name]], temp)
      sobjs[[sample_name]] <- subset(sobjs[[sample_name]], subset= DF.DROPLET.TYPE == 'Singlet')
      sobjs[[sample_name]]@active.assay = 'RNA'
      sobjs[[sample_name]][['SCT']] <- NULL
      sobjs[[sample_name]] <- NormalizeData(sobjs[[sample_name]],
                                            assay='RNA',
                                            normalization.method = "LogNormalize",
                                            scale.factor = 10000)
      if ('ADT' %in% names(sobjs[[sample_name]]@assays)){
        cat(paste0('found ADT assay for sample: ', sample_name, "\n"))
        sobjs[[sample_name]] <- NormalizeData(sobjs[[sample_name]],
                                              assay='ADT',
                                              normalization.method = "CLR")
      }
      sobjs[[sample_name]]@meta.data$LIBRARY <- sample_name
    }
    setwd("/krummellab/data1/DSCoLab/XAUT1/10x/merged_SCG1_10")
    sobjs.orig <- sobjs
    gene_intersection <- lapply(sobjs, row.names) %>% Reduce(intersect, .)

    first_sobj <- sobjs[[names(sobjs)[1]]]
    sobjs[[names(sobjs)[1]]] <- NULL
    merged_data <- merge(x = first_sobj, y = unname(sobjs))

    merged_data <- FindVariableFeatures(merged_data, selection.method = "vst", nfeatures = 3000)
    merged_data <- ScaleData(merged_data, vars.to.regress=c('percent.mt', 'percent.ribo', 'S.Score', 'G2M.Score'),
                             verbose = FALSE)
    merged_data <- RunPCA(merged_data, npcs = 30, verbose = FALSE)

    save(merged_data, file=paste0(prefix, '_merged_temp.RData'))
  } else {
    load(paste0(prefix, '_merged_temp.RData'))
  }
  
  pdf(paste0(prefix, '_harmony_convergence.pdf'))
  merged_data <- RunHarmony(merged_data, assay.use='RNA', group.by.vars="LIBRARY", plot_convergence = TRUE,max.iter.harmony=30, max.iter.cluster=30)
  dev.off()
  save(merged_data, file=paste0(prefix, '_merged_harmonized.RData'))
} else {
  load(paste0(prefix, '_merged_harmonized.RData'))  
}

merged_data <- RunUMAP(merged_data,
                       reduction='harmony',
                       dims = 1:30,  # Num PCs to use
                       n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                       min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                       spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                       a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                       b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                       verbose = FALSE)

# Calculate the neighborhood graph
merged_data <- FindNeighbors(merged_data,
                             reduction='harmony',
                             dims = 1:30,  # Num PCs to use
                             k.param = 20,  # k for the knn algorithm
                             verbose = FALSE
                             )

merged_data$SAMPLE.by.SNPs <- paste(gsub('XAUT1-POOL-SCG', 'P.S', merged_data$LIBRARY), merged_data$SAMPLE.by.SNPs, sep='.C')

png(filename=paste0(prefix, '_library_umap.png'), width = 5, height = 5, units = "in", res = 300)
print(DimPlot(merged_data, group.by='LIBRARY') + NoLegend())
dev.off()

png(filename=paste0(prefix, '_split_library_umap.png'), width = 12, height = 6, units = "in", res = 300)
print(DimPlot(merged_data, split.by='LIBRARY', ncol=4) + theme(legend.position="none", axis.title=element_blank(), axis.text=element_blank()))
dev.off()

png(filename=paste0(prefix, '_sample_umap.png'), width = 5, height = 5, units = "in", res = 300)
print(DimPlot(merged_data, group.by='SAMPLE.by.SNPs', label=T) + NoLegend())
dev.off()

png(filename=paste0(prefix, '_split_sample_umap.png'), width = 18, height = 42, units = "in", res = 300)
print(DimPlot(merged_data, split.by='SAMPLE.by.SNPs', ncol=6) + theme(legend.position="none", axis.title=element_blank(), axis.text=element_blank()))
dev.off()


for (res in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2)){
  if (paste0('louvain_res', res) %in% colnames(merged_data@meta.data)){
    next
  }
  merged_data <- FindClusters(merged_data, verbose = TRUE,
                              algorithm = 1,
                              resolution = res)
  merged_data@meta.data[[paste0('louvain_res', res)]] <- merged_data@meta.data$seurat_clusters
  
  png(filename=paste0(prefix, '_louvain_res_', res, '.png'), width = 5, height = 5, units = "in", res = 300)
  print(DimPlot(merged_data, group.by=paste0('louvain_res', res), label=T) + NoLegend())
  dev.off()
}

save(merged_data, file=paste0(prefix, '_merged_processed.RData'))

metadata <- merged_data@meta.data
for (redn in c('pca', 'umap')){
  cell_embeddings <- as.data.frame(merged_data@reductions[[redn]]@cell.embeddings[, c(1: min(5, dim(merged_data@reductions[[redn]]@cell.embeddings)[2]))])
  metadata <- merge(metadata, 
                    cell_embeddings,
                    by=0)
  rownames(metadata) <- metadata$Row.names
  metadata$Row.names <- NULL
}

write.table(metadata,
            file=paste(prefix, 'metadata.tsv', sep='_'),
            sep='\t',
            row.names=T,
            col.names=T,
            quote=F)

final_res <- 0.2
plots <- list()
for (l in levels(metadata[[paste0('louvain_res', final_res)]])){
  colors <- rep('lightgrey', length(levels(metadata[[paste0('louvain_res', final_res)]])))
  colors[as.numeric(l)+1] <- 'red'
  alphas <- rep(0.25, length(levels(metadata[[paste0('louvain_res', final_res)]])))
  alphas[as.numeric(l)+1] <- 1
  plots[[l]] <- ggplot(metadata, aes_string(x="UMAP_1", y="UMAP_2", col=paste0('louvain_res', final_res))) + 
                  geom_point(size=0.05) + 
                  scale_color_manual(values=colors) +
                  scale_alpha_manual(values=alphas) +
                  annotate("text", label = l, x = -Inf, y = Inf, size = 10, colour = "red", hjust = 0, vjust = 1)+
                  theme_bw() +
                  theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
}

png(paste0(prefix, '_res_', final_res, '_per_cluster_dimplots.png'), width=15, height=(3 * ceiling(length(plots)/5)), units='in', res=150)
plot_grid(plotlist=plots, ncol=5)
dev.off()

plots <- list()
metadata$LIBRARY <- as.factor(metadata$LIBRARY)
libraries <- setNames(1:length(levels(metadata$LIBRARY)), levels(metadata$LIBRARY))
for (l in names(libraries)){
  colors <- rep('lightgrey', length(libraries))
  colors[libraries[l]] <- 'red'
  alphas <- rep(0.25, length(libraries))
  alphas[libraries[l]] <- 1
  plots[[l]] <- ggplot(metadata, aes(x=UMAP_1, y=UMAP_2, col=LIBRARY, alpha=LIBRARY)) + 
                  geom_point(size=0.05) + 
                  scale_color_manual(values=colors) +
                  scale_alpha_manual(values=alphas) +
                  annotate("text", label = l, x = -Inf, y = Inf, size = 5, colour = "red", hjust = 0, vjust = 1)+
                  theme_bw() +
                  theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
}

png(paste0(prefix, '_per_LIBRARY_dimplots.png'), width=15, height=(3 * ceiling(length(plots)/5)), units='in', res=150)
plot_grid(plotlist=plots, ncol=5)
dev.off()

Idents(merged_data) <- merged_data@meta.data[[paste0('louvain_res', final_res)]]
markers <- FindAllMarkers(merged_data, 
                          test.use='poisson', 
                          latent.vars='LIBRARY', 
                          assay='RNA', 
                          logfc.threshold=0.4, 
                          min.pct=0.2, 
                          only.pos=TRUE)
write.table(markers, file=paste0(prefix, '_res_', final_res, '_markers.tsv'), sep='\t', col.names=T, row.names=F, quote=F)

top_markers = (markers %>% group_by(cluster) %>% top_n(5, wt=avg_logFC))$gene
genes_to_consider = c()
for (gene in top_markers){
  if (!gene %in% genes_to_consider){
    genes_to_consider <- c(genes_to_consider, gene)
  }
}

png(paste(prefix, 'res', final_res, '_TopMarkerDotPlot.png', sep='_'), 
    height=50 * length(genes_to_consider), 
    width=200 + (100 * length(unique(as.vector(merged_data@meta.data[[paste0('louvain_res', final_res)]])))), 
    units='px', 
    res=150)
print(DotPlot(merged_data, 
              features=genes_to_consider, 
              group.by=paste0('louvain_res', final_res),
              assay='RNA',
              cols='RdYlBu') + coord_flip())
dev.off()



