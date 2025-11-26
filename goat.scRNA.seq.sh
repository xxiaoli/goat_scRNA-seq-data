
##Load required R packages
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(sctransform)
library(cowplot)
library(DoubletFinder)
library(harmony)
library(ape)
library(clustree)

## Batch loading of raw data files
setwd("/single_cell/Seurat")
dir =c("B1_Horn_buds_1_EmptyDrops_CR_matrix",
	        "B1_Horn_buds_2_EmptyDrops_CR_matrix",
	        "B1_Horn_buds_3_EmptyDrops_CR_matrix",
	       "B1_Skin_1_EmptyDrops_CR_matrix",
	       "B1_Skin_2_EmptyDrops_CR_matrix",
	       "B1_Skin_3_EmptyDrops_CR_matrix",
	       "B7_Horn_buds_1_EmptyDrops_CR_matrix",
	       "B7_Horn_buds_2_EmptyDrops_CR_matrix",
	       "B7_Horn_buds_3_EmptyDrops_CR_matrix",
	        "B7_Skin_1_EmptyDrops_CR_matrix",
	       "B7_Skin_3_EmptyDrops_CR_matrix",
           "B21_2_Horn_buds_EmptyDrops_CR_matrix",
           "B21_3_Horn_buds_EmptyDrops_CR_matrix",
	       "B21_1_Skin_EmptyDrops_CR_matrix",
	       "B21_2_Skin_EmptyDrops_CR_matrix",
	       "B21_3_Skin_EmptyDrops_CR_matrix"
	       ) 

samples_name =c('H1.1', 'H1.2', 'H1.3', 'S1.1', 'S1.2','S1.3','H7.1', 'H7.2', 'H7.3', 'S7.1', 'S7.3', 'H21.2', 'H21.3', 'S21.1', 'S21.2','S21.3')
#Generation of Seurat objects
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
scRNAlist[[i]] <- CreateSeuratObject(counts, project= samples_name[i],
                                       min.cells=3)
  samInfo <- rep(NA,nrow(scRNAlist[[i]]@meta.data))
  scRNAlist[[i]] <- AddMetaData(object = scRNAlist[[i]], metadata = rep(samples_name[i],ncol(counts)), col.name = "samInfo")
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])

  #Calculate the proportion of mitochondrial genes
  if(T){    
    mt.genes <- c("ND6","ND5","ND4L","ND4","ND3","ND2","ND1","CYTB","COX3","COX2","COX1","ATP8","ATP6")
    mt.genes <- CaseMatch(mt.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], features=mt.genes)
  }
  #Calculate the proportion of ribosomal genes
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
}

setwd("/Results")
##Merge multiple scRNA-seq objects into a single Seurat object using the merge function
names(scRNAlist) <- samples_name
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNA
head(scRNA,n=10)
table(scRNA$orig.ident)
y <- table(scRNA$orig.ident)
saveRDS(scRNA,file = 'B21_B7_B1_orig.Raw_data.rds')

scRNA <- subset(scRNA, subset = nCount_RNA > 500 & nFeature_RNA > 500 & percent.mt < 10)
saveRDS(scRNA,file = 'B21_B7_B1_clean_data_filter_gene.rds')
#作图
samples_order =c('H1.1', 'H1.2', 'H1.3', 'S1.1', 'S1.2','S1.3','H7.1', 'H7.2', 'H7.3', 'S7.1', 'S7.3', 'H21.2', 'H21.3', 'S21.1', 'S21.2','S21.3')
scRNA@meta.data$orig.ident <- factor(scRNA@meta.data$orig.ident,levels = samples_order)
p1 <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA","percent.rb", "percent.mt", "percent.HB"), group.by="orig.ident", ncol = 5) + NoLegend()
ggsave("B21_B7_B1_clean_data_filter_gene_qc.pdf", plot = p1, width = 24, height = 6)

##
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
ScaleData(features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
pc.num=1:30
scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters(resolution=0.6)
saveRDS(scRNA,file = 'B21_B7_B1_clean_data_filter_gene_Nor.rds')
#作图
p1<- DimPlot(scRNA, reduction = "tsne", group.by = "orig.ident",raster=FALSE, pt.size = 0.25)
p2<- DimPlot(scRNA, reduction = "tsne", raster=FALSE, pt.size = 0.25)
p1+p2
ggsave("B21_B7_B1_clean_data_filter_gene_have_doublts_r0.6_tsne.pdf",,height=6, width=14)
p1<- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident",raster=FALSE, pt.size = 0.25)
p2<- DimPlot(scRNA, reduction = "umap",raster=FALSE, pt.size = 0.25)
p1+p2
ggsave("B21_B7_B1_clean_data_filter_gene_have_doublts_r0.6_umap.pdf",,height=6, width=14)

####Removal of double cells
sweep.res.list_scRNA <- paramSweep_v3(scRNA, PCs = 1:30, sct = FALSE)
sweep.stats_scRNA <- summarizeSweep(sweep.res.list_scRNA, GT = FALSE)
bcmvn_scRNA <- find.pK(sweep.stats_scRNA)
opt_pK <- as.numeric(as.vector(bcmvn_scRNA$pK[which.max(bcmvn_scRNA$BCmetric)]))
print(opt_pK)
annotations <- scRNA@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.06*nrow(scRNA@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
scRNA <- doubletFinder_v3(scRNA,PCs = 1:30,pN = 0.25,pK = opt_pK,nExp = nExp_poi,reuse.pANN = FALSE,sct = FALSE)
print(nExp_poi)
saveRDS(scRNA, file = "B21_B7_B1_clean_single_doublts_ing.rds")
name <- colnames(scRNA@meta.data)[11]
Idents(scRNA) <- name
scRNA <- subset(scRNA, idents = "Singlet")
saveRDS(scRNA, file = "B21_B7_B1_clean_afterdoublts.rds")

##############CCA
scRNA <- readRDS("B21_B7_B1_clean_afterdoublts.rds")

scRNAlist <- SplitObject(scRNA, split.by = "orig.ident")
head(scRNAlist,n=11)
merged.anchors <- FindIntegrationAnchors(object.list = scRNAlist, dims = 1:25)
merged <- IntegrateData(anchorset = merged.anchors, dims = 1:25)
DefaultAssay(merged) <- "integrated"
merged <- ScaleData(merged, features = rownames(merged))
merged <- RunPCA(merged, npcs = 25, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:25)
saveRDS(merged, file = "B21_B1_delet_H21.1_r0.6_CCA_integrat.ing.rds")
merged <- FindClusters(merged, resolution = 0.6)
merged <- RunUMAP(merged, reduction = "pca",dims = 1:25)
merged <- RunTSNE(merged, reduction = "pca",dims = 1:25)
saveRDS(merged, file = "B21_B7_B1_CCA_integrated.rds")

####Find marker
scRNA<-readRDS(file="B21_B7_B1_CCA_integrated.rds")
merged.markers <- FindAllMarkers(scRNA, only.pos = TRUE)
merged.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
write.csv(merged.markers,file = "B21_B7_B1_CCA_integrated_allmarker.csv")

###GO enrichment
library(Seurat)
library(BiocManager)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(msigdbr)
library("org.Hs.eg.db")
scRNA <- readRDS("B21_B7_B1_CCA_celltype_integrated.rds")
head(scRNA@meta.data)
DefaultAssay(scRNA) <- 'RNA'
Idents(scRNA) <- 'integrated_merge_cluster'
only.pos = TRUE
logfc.threshold = 0.25
all.clusters.markers <- FindAllMarkers(scRNA, only.pos = T, logfc.threshold = logfc.threshold)
write.csv(all.clusters.markers,file = "B21_B7_B1_CCA_integrated_allmarker.csv")
#head(all.clusters.markers)
top600 <- all.clusters.markers %>% group_by(cluster) %>% top_n(n = 600, wt = avg_log2FC) 
top600_gene<-unstack(top600,gene~cluster)
top600_entrez <- lapply(X = top600_gene, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})

write.csv(top600,file = "B21_B7_B1_CCA_integrated_top600_allmarker.csv")

### ALL biological process
all_bp <- compareCluster(top600_entrez, 
                         fun='enrichGO',
                         ont= 'all',
                         OrgDb='org.Hs.eg.db')
write.csv(as.data.frame(all_bp),"B21_B7_B1_CCA_integrated_top600_all_GO.csv")

cell_type_order <- c("Epidermal_lineage_cells","Dermal_lineage_cells","Pericytes","Smooth_muscle_cells","Endothelial_cells","SG_Luminal_cells","SG_Basal_cells","Melanocytes","Immune_cells","Schwann_cells")

all_bp_sorted <- all_bp %>%
  mutate(Cluster = factor(Cluster, levels = cell_type_order)) %>%
  arrange(Cluster)
p <- ggplot(all_bp_sorted, aes(y = Description, x = Cluster), includeAll = FALSE, showCategory = 5) +
  geom_point(aes(color = -log10(p.adjust), size = Count)) + 
  scale_color_gradientn(colours = c('green', 'orange', 'red')) +
  labs(color = expression(-log10(p.adjust)), size = "Count", 
       x = "Gene Number", y = " ", title = "GO Pathway Enrichment") +
  theme_bw() +
  theme(text = element_text(size = 28, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p <- p + scale_size(range = c(3, 15)) 
print(p)
ggsave("B21_B7_B1_CCA_integrated_top600_all_GO.pdf", plot = p, width = 16, height = 18)


#KEGG enrichment analysis
all_kegg <- compareCluster(top600_entrez,
                           fun='enrichKEGG',
                           pvalueCutoff=0.05, 
                           organism="hsa"
)

write.csv(as.data.frame(all_kegg),"B21_B7_B1_CCA_integrated_top600_all_KEGG.csv")

###cell order
all_kegg_sorted <- all_kegg %>%
  mutate(Cluster = factor(Cluster, levels = cell_type_order)) %>%
  arrange(Cluster)
p <- ggplot(all_kegg_sorted,aes(y=Description,x=Cluster),includeAll = FALSE,showCategory = 5)+
  geom_point(aes(color=-log10(p.adjust),size=Count))+
  scale_color_gradientn(colours = c('green', 'orange', 'red'))+
  labs(color=expression(-log10(p.adjust),size="Count"), 
       x="Gene Number",y=" ",title="KEGG Pathway Enrichment")+
  theme_bw()+
  theme(text = element_text(size=20, color = "black"))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  theme()
p <- p + scale_size(range = c(3, 10)) # 你可以根据需要调整大小范围
print(p)
ggsave("B21_B7_B1_CCA_integrated_top600_all_KEGG2.pdf",width = 12,height = 11)

#########Dermal subclusters
library(parallel)
library(monocle)
library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(ggsignif)
library(patchwork)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(msigdbr)
library("org.Hs.eg.db")
scRNA<-readRDS(file="B21_B7_B1_CCA_celltype_integrated.rds")
scRNA<- subset(scRNA,integrated_merge_cluster %in% c("Dermal_lineage_cells"))
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA = scRNA %>% RunHarmony("orig.ident",plot_convergence = FALSE)
scRNA <- scRNA %>%
  RunTSNE(reduction = "harmony", dims = 1:10) %>%
  RunUMAP(reduction = "harmony", dims = 1:10) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = 0.2) %>%
  identity()
saveRDS(scRNA, file = "r0.2_nfe2000_dim10_Horn_skin_Dermal_umap.rds")


##
allcolour=c("#FFA500","#9370DB","#98FB98","#FA8072","#87CEEB","#40E0D0")
scRNA$group <- factor(scRNA$group, levels = c("H1","H7","H21","S1","S7","S21")) 

p1<-DimPlot(scRNA, reduction = "umap",group.by = "group", cols=allcolour,label = TRUE,repel = T) 
p2<-DimPlot(scRNA, reduction = "umap", group.by = "RNA_snn_res.0.2",label = TRUE,repel = T) 
p1+p2
ggsave("r0.2_nfe2000_dim10_Horn_skin_Dermal_umap.pdf",height=3.5, width=8)
p <- p1 + p2
pdf("r0.2_nfe2000_dim10_Horn_skin_Dermal_celltype_umap.pdf", height = 3.5, width = 8)
print(p)
dev.off()
gene_order2 <- c("LUM","DCN","MFAP5","APOD","CRABP1","APCDD1","PTN","UCHL1","IBSP","BGLAP","COL2A1","CHAD","ZEB2","DLC1")

p <- DotPlot(scRNA, features = gene_order2, cols = c("#DCDCDC", "#DC143C"),
             col.min = 0.5, col.max = 5, dot.min = 0, dot.scale = 8) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
           axis.text.y = element_text(size = 12))

pdf("r0.2_nfe2000_dim10_Horn_skin_Dermal_DotPlot.pdf", height = 4, width = 10)
print(p)
dev.off()

###
Idents(scRNA) <- "RNA_snn_res.0.2"
levels(Idents(scRNA))
scRNA@meta.data$RNA_snn_res.0.2 <- scRNA@meta.data$RNA_snn_res.0.2
celltype <-          c("0"="OCPs", #Mesenchymal_stem_cells
                       "1"="Papi_Fib", #Papillary_Fib
                       "2"= "DP", #Dermal_papilla_cells
                       "3"= "Reti_Fib",#Reticular_Fib
                       "4"= "Osteoblasts", 
                       "5"= "ZEB2+_Fib", 
                       "6"= "Chondrocytes" 
                       )
names(celltype) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, celltype)
scRNA$celltype <- Idents(scRNA)
head(scRNA@meta.data)
saveRDS(scRNA, file = "Dermal_celltype.rds")

###monocle2
scRNA <- readRDS("Dermal_celltype.rds")
scRNA <- subset(scRNA,tissue_type %in% c("Horn_buds"))
scRNA <- subset(scRNA,celltype %in% c("Chondrocytes","Osteoblasts","OCPs","ZEB2+_Fib"))
N=length(colnames(scRNA))/2
N=round(N)
scRNA<-scRNA[,sample(x=colnames(scRNA),size = N,replace=F)]
data <- as(as.matrix(scRNA[["RNA"]]$counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA@meta.data)
fData <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 0.3 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
suppressWarnings(mycds <- orderCells(mycds))
saveRDS(mycds, file = "Dermal_celltype_monocle2.rds")
mycds@phenoData@data$group <- factor(mycds@phenoData@data$group, levels = c("H1", "H7", "H21"))
mycds@phenoData@data$celltype <- factor(mycds@phenoData@data$celltype, levels = c("ZEB2+_Fib","OCPs","Osteoblasts","Chondrocytes"))
p1 <- plot_cell_trajectory(mycds,color_by="Pseudotime",cell_size=0.2,show_backbone=TRUE)
p2 <- plot_cell_trajectory(mycds,color_by="celltype",cell_size=0.2,show_backbone=TRUE)
p3 <- plot_cell_trajectory(mycds, color_by = "State",cell_size=0.2,show_backbone=TRUE)
p<- p1+p2+p3
pC <- p1 + p2 + p3 + plot_layout(ncol = 3)
ggsave(pC,file="Dermal_celltype_monocle2.pdf",width = 12, height = 4)

###GO
BEAM_res=BEAM(mycds,branch_point = 1,cores = 8, , progenitor_method = "duplicate")

BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "Dermal_celltype_monocle2_branched_heatmap.rds")
pdf("Dermal_celltype_monocle2_branched_heatmap.pdf",width = 4,height = 4.5)
tmp1<-plot_genes_branched_heatmap(mycds[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 3, 
                                 cores = 8,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T 
)
tmp1$ph_res
dev.off()
gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
write.csv(allcluster_go,  file = "Dermal_celltype_monocle2_branched_heatmap_GO.csv",row.names = FALSE)

##########TF
####R
sce<-readRDS(file="Dermal_celltype_monocle2.rds")
sce<- subset(sce,tissue_type %in% c("Horn_buds"))
sce<- subset(sce,celltype %in% c("Chondrocytes","Osteoblasts","OCPs","ZEB2+_Fib")) 
sce_count <- GetAssayData(sce[["RNA"]], layer='counts')
write.csv(sce_count, file="sce_count_Horn_Skin.csv",quote=F)
sce_count_transposed <- t(sce_count)
 write.csv(sce_count_transposed, file="sce_count_transposedHorn_Skin.csv", quote=FALSE)

conda activate pyscenic2 
cat >change.py
import os,sys
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("sce_count_transposedHorn_Skin.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("sce_count.loom",x.X.transpose(),row_attrs,col_attrs);
ctrl+D
python change.py

######scenic
/Software/miniconda/miniconda3/envs/pyscenic2/bin/pyscenic grn --num_workers 20 --output /data/scenic/Horn_skin_AllCell_grn.tsv --method grnboost2  /data/scenic/sce_count.loom /data/.scenic/allTFs_hg38.txt

/Software/miniconda/miniconda3/envs/pyscenic2/bin/pyscenic ctx  /data/scenic/Horn_skin_AllCell_grn.tsv /data/scenic/Up_500bp_Down_100bp.regions_vs_motifs.rankings.feather --annotations_fname /data/scenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl  --expression_mtx_fname /data/scenic/sce_count.loom  --mode "dask_multiprocessing" --output /data/scenic/Horn_skin_AllCell_reg.csv --num_workers 20 --mask_dropouts

/Software/miniconda/miniconda3/envs/pyscenic2/bin/pyscenic aucell  /data/scenic/sce_count.loom /data/scenic/Horn_skin_AllCell_reg.csv --output /data/scenic/Horn_skin_aucell_goat.loom --num_workers 20




###cellchat
library(CellChat)
library(patchwork)
library(Seurat)
library(cpp11)
library(igraph)
library(ComplexHeatmap)
library(glue)
library(purrr)
library(ggplot2)  
scRNA <- readRDS("Dermal_celltype.rds")
scRNA1 <- readRDS("Epidermal_celltype.rds")
scRNA<-merge(scRNA,scRNA1)
saveRDS(scRNA, file = "Horn_skin_Dermal_Epidermal_celltype.rds")
scRNA <-subset(scRNA, celltype %in% c("Basal.Kc_I","Basal.Kc_II","Spinous.Kc_I","Spinous.Kc_II","Reti_Fib","Papi_Fib","DP","ZEB2+_Fib","OCPs","Osteoblasts","Chondrocytes"))
H1<-subset(scRNA,group %in% c("H1"))
H1<-subset(H1, celltype %in% c("Basal.Kc_I","Basal.Kc_II","Spinous.Kc_I","Spinous.Kc_II","Reti_Fib","Papi_Fib","DP","ZEB2+_Fib"))
H7<-subset(scRNA,group %in% c("H7"))
H21<-subset(scRNA,group %in% c("H21"))
S1<-subset(scRNA,group %in% c("S1"))
S7<-subset(S7, celltype %in% c("Basal.Kc_I","Basal.Kc_II","Spinous.Kc_I","Spinous.Kc_II","Reti_Fib","Papi_Fib","DP","ZEB2+_Fib"))
S7<-subset(scRNA,group %in% c("S7"))
S7<-subset(S7, celltype %in% c("Basal.Kc_I","Basal.Kc_II","Spinous.Kc_I","Spinous.Kc_II","Reti_Fib","Papi_Fib","DP","ZEB2+_Fib"))
S21<-subset(scRNA,group %in% c("S21"))
S21<-subset(S21, celltype %in% c("Basal.Kc_I","Basal.Kc_II","Spinous.Kc_I","Spinous.Kc_II","Reti_Fib","Papi_Fib","DP","ZEB2+_Fib"))
data1 <- GetAssayData(H1[['RNA']], layer = "data.1")
data2 <- GetAssayData(H1[['RNA']], layer = "data.2")
data.input <- cbind(data1, data2)
meta =  H1@meta.data
unique(meta$group)
data.input[1:6,1:3]
head(meta);table(meta$group) #含normal(NL)和diseases(LS)
identical(rownames(meta),colnames(data.input)) 
unique(meta$celltype) 
meta$celltype <- as.factor(meta$celltype)
meta$celltype <- droplevels(meta$celltype)
print(levels(meta$celltype))
meta$celltype <- droplevels(meta$celltype, exclude = setdiff(levels(meta$celltype), unique(meta$celltype)))
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = 'celltype') 

cellchat <- setIdent(cellchat, ident.use = 'celltype') 
levels(cellchat@idents)
table(cellchat@idents) 
CellChatDB <- CellChatDB.human #有人类(CellChatDB.human)和小鼠(CellChatDB.mouse)
p<- showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.us
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat,
                                       only.pos = TRUE, 
                                       thresh.pc = 0, 
                                       thresh.fc = 0,
                                       thresh.p = 0.05) 
head(cellchat@var.features$features) 
head(cellchat@var.features$features.info) 
cellchat <- identifyOverExpressedInteractions(cellchat)
head(cellchat@LR$LRsig) 
cellchat <- computeCommunProb(cellchat, raw.use = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10) 
df.net <- subsetCommunication(cellchat, slot.name = 'net')
head(df.net) 
cellchat <- computeCommunProbPathway(cellchat) 
df.netp <- subsetCommunication(cellchat, slot.name = 'net') 
head(df.netp)
write.csv(df.netp, file = "H1_Horn_skin_Dermal_Epidermal_celltype_net_all_cellchat.csv")
df.netp <- subsetCommunication(cellchat, slot.name = 'netP') 
head(df.netp)
write.csv(df.netp, file = "H1_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.csv")
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "H1_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.rds")

library(CellChat)
library(ComplexHeatmap)
library(patchwork)
library(purrr)
cellchat.H1 <- readRDS("H1_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.rds")
cellchat.H7 <- readRDS("H7_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.rds")
cellchat.H21 <- readRDS("H21_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.rds")
cellchat.S1 <- readRDS("S1_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.rds")
cellchat.S7 <- readRDS("S7_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.rds")
cellchat.S21 <- readRDS("S21_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.rds")
levels_list <- list(
  levels(cellchat.H1@idents),
  levels(cellchat.H7@idents),
  levels(cellchat.H21@idents),
  levels(cellchat.S1@idents),
  levels(cellchat.S7@idents),
  levels(cellchat.S21@idents)
)

all_identical <- all(sapply(levels_list[-1], identical, levels_list[[1]]))
all_identical
object.list <- list(S1 = cellchat.S1,
                    S7 = cellchat.S7,
                    S21 = cellchat.S21,
                    H1 = cellchat.H1, 
                    H7 = cellchat.H7,
                    H21 = cellchat.H21) 
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
cellchat
dplyr::glimpse(cellchat)
saveRDS(cellchat, file = "group_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.rds")
cellchat <- readRDS("group_Horn_skin_Dermal_Epidermal_celltype_all_cellchat.rds")
p1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
p2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p<-p1 + p2
ggsave(p, file="group_stage_Number_of_inferred_interactions_epi_dermal.pdf",height=2.5, width=4

###outgoing
object.list <- list(
  S1 = cellchat.S1,
  S7 = cellchat.S7,
  S21 = cellchat.S21,
  H1 = cellchat.H1, 
  H7 = cellchat.H7,
  H21 = cellchat.H21
)
for (i in seq_along(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}
pathway.union <- Reduce(union, lapply(object.list, function(x) x@netP$pathways))
pdf("group.signalingRole_heatmap_ootgoing.pdf", width = 5, height = 5)
par(mfrow = c(2, 3))
for (i in seq_along(object.list)) {
  ht <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 8.5)
  draw(ht, ht_gap = unit(0.5, "cm"))
}
dev.off()
###incoming

object.list <- list(
  S1 = cellchat.S1,
  S7 = cellchat.S7,
  S21 = cellchat.S21,
  H1 = cellchat.H1, 
  H7 = cellchat.H7,
  H21 = cellchat.H21
)
for (i in seq_along(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}
pathway.union <- Reduce(union, lapply(object.list, function(x) x@netP$pathways))
pdf("group.signalingRole_heatmap_incoming.pdf", width = 5, height = 5)
par(mfrow = c(2, 3))
for (i in seq_along(object.list)) {
  ht <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 8.5,color.heatmap = "GnBu")
  draw(ht, ht_gap = unit(0.5, "cm"))
}
dev.off()





###cross-species
###RPCA
ifnb <- readRDS("goat_deer.rds")
ifnb.list <- SplitObject(ifnb, split.by = "group")
gene_common <- Reduce(intersect, lapply(ifnb.list, function(x) rownames(x[["RNA"]])))
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- subset(x, features = gene_common)
  return(x)
})

# Step 1: Normalize and identify HVGs
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

# Step 2: Select integration features
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2000)

# Step 3: Scale & PCA
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  return(x)
})

# Step 4: RPCA
immune.anchors <- FindIntegrationAnchors(
  object.list = ifnb.list,
  anchor.features = features,
  reduction = "rpca"
)

immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 25, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:25)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:25)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)
saveRDS(immune.combined, file = "goat_deer_RPCA_integrated.rds")

p1<- DimPlot(immune.combined, reduction = "umap", group.by = "group")
p2<- DimPlot(immune.combined, reduction = "umap", group.by = "integrated_snn_res.0.8",label.size = 3.2,label = TRUE)
p<-p1+p2
ggsave("r0.8_nfe2000_dim25_RPCA_goat_deer_anlter_umap.pdf", height = 5, width = 12)

p<-DimPlot(immune.combined, reduction = "umap",split.by = "group",label=TRUE)
ggsave("r0.8_nfe2000_dim35_RPCA_goat_deer_anlter_group_umap.pdf", height = 5, width = 12)



#####find marker

immune.combined <- JoinLayers(immune.combined)
Idents(immune.combined) <- "integrated_snn_res.0.8"
DefaultAssay(immune.combined) <- "RNA"

immune.combined.markers <- FindAllMarkers(
  immune.combined,
  only.pos = TRUE,
  logfc.threshold = 0.25,  
  min.pct = 0.1,           
  min.diff.pct = 0.05       
)
write.csv(immune.combined.markers, file = "r0.5_nfe2000_dim35_RPCA_goat_deer_anlter_allmarker1.csv")




