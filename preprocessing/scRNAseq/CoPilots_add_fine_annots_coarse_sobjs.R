source("/krummellab/data1/chrisquatjr/candersen_scripts/RNA_seq_functions.R")

project <- "XAUT2"
if (project == "XAUT1"){
  prefix = "FINAL_"
} else {
  prefix = "Analyzed_"
}
dataset <- "merged_SCG9_12"
datadir <- paste0("/krummellab/data1/DSCoLab/",project,"/10x/",dataset)
setwd(datadir)
sobj <- loadRData(paste0(prefix,dataset,"_annotated.RData"))

# Add in placeholder column for all fine annotations
sobj@meta.data$fine_annotations <- "NOT CONSIDERED"

for (folder in dir(paste0(datadir,"/subclusters/"))){
  setwd(paste0(datadir,"/subclusters/",folder))
  if (folder == "Monocytes"){
    subset_name <- "MONO"
  } else {
    subset_name <- toupper(unlist(strsplit(folder,split="_"))[1])
  }
  subsobj <- loadRData(paste0(prefix,dataset,"_",subset_name,"_annotated.RData"))
  # Add in fine annotations for all subset's barcodes
  sobj@meta.data[rownames(sobj@meta.data) %in% rownames(subsobj@meta.data),"fine_annotations"] <- paste0(subset_name,".",as.character(subsobj@meta.data$fine_annotations))
}

setwd(datadir)
png(paste0(prefix,dataset,"_all_annotations_umap.png"),width=20,height=12,res=150,units="in")
print(dittoDimPlot(sobj,var="fine_annotations"))
dev.off()

pdf(paste0(prefix,dataset,"_all_annotations_umap.pdf"),width=20,height=12,useDingbats=F)
print(dittoDimPlot(sobj,var="fine_annotations"))
dev.off()

save(sobj,file=paste0(prefix,dataset,"_annotated.RData"))
write.table(sobj@meta.data,file=paste0(prefix,dataset,"_annotated_metadata.tsv"),sep="\t",quote=F)
