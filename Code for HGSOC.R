#GSE158937 analysis
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
pbmc.data1 <- Read10X(data.dir = "Patient_1 working directory")
pbmc1 <- CreateSeuratObject(counts = pbmc.data1, project = "Patient_1", 
                            min.cells = 3,min.features = 100,names.delim = "_sc")
pbmc.data2 <- Read10X(data.dir = "Patient_2 working directory")
pbmc2 <- CreateSeuratObject(counts = pbmc.data2, project = "Patient_2", 
                            min.cells = 3,min.features = 100,names.delim = "_sc")
pbmc.data3 <- Read10X(data.dir = "Patient_3 working directory")
pbmc3 <- CreateSeuratObject(counts = pbmc.data3, project = "Patient_3", 
                            min.cells = 3,min.features = 100,names.delim = "_sc")
pbmc<-merge(pbmc1,c(pbmc2,pbmc3))
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 50) 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1500)
top10 <- head(VariableFeatures(pbmc), 10) 
head(pbmc$RNA@var.features,10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimHeatmap(pbmc, dims = 1:18, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:18)
JackStrawPlot(pbmc, dims = 1:18)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20) 
pbmc <- FindClusters(pbmc, resolution = 0.8)
head(Idents(pbmc), 10)
pbmc <- RunTSNE(pbmc, dims = 1:20)
DimPlot(pbmc,  pt.size = 1,reduction = "tsne")
DimPlot(pbmc, pt.size = 1, reduction = "tsne",label = TRUE)
LabelClusters(DimPlot(pbmc, reduction = "tsne"),id = 'ident')
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

library(SingleR)
hpca.se <- HumanPrimaryCellAtlasData()   #参考数据集
pred.hesc <- SingleR(test =pbmc@assays$RNA@counts, ref = hpca.se, labels = hpca.se$label.main)
new.cluster.ids <- c("Epithelial cell",
                     "Epithelial cell",
                     "Epithelial cell",
                     "Epithelial cell",
                     "Epithelial cell",
                     "Epithelial cell",
                     "Macrophage",
                     "Epithelial cell",
                     "Epithelial cell",
                     "Epithelial cell",
                     "T cell",
                     "Fibroblast",
                     "Epithelial cell",
                     "Endothelial cell",
                     "Monocyte",
                     "NK cell",
                     "Fibroblast",
                     "Epithelial cell",
                     "Macrophage",
                     "B cell")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()


#Scissor analysis
setwd("working directory")
bulk_dataset<-read.csv("expression.csv",header = T,row.names = 1)
bulk_dataset<-as.matrix(bulk_dataset)
bulk_phenotype<-read.csv("clinical.CSV",header = T,row.names = 1)
table(bulk_phenotype)
phenotype<-bulk_phenotype[,1]
tag <- c('Sensitive', 'Resistant')
infos <- Scissor(bulk_dataset, pbmc, phenotype, tag = tag, alpha = 0.05,
                  family = "binomial")
sc_dataset<-pbmc
length(infos$Scissor_pos)
Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos$Scissor_neg] <- 1
Scissor_select[infos$Scissor_pos] <- 2
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(sc_dataset, reduction = 'tsne', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))


#InferCNV analysis
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
library(Seurat)
library(infercnv)
library(miscTools)
expFile=pbmc@assays[["RNA"]]@counts
expFile<-data.frame(expFile)
geneFile<-read.csv('geneInfor.csv',header = T,row.names = 1)
groupinfo=read.csv('group.CSV',header = T)
groupFiles='groupFiles.txt'
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c("T_cells","Monocyte","NK_cell","Macrophage","B_cell")) 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir="infercnv_output",
                             cluster_by_groups=T, 
                             denoise=TRUE,
                             HMM=F,hclust_method="ward.D2", plot_steps=F)
# CNV score analysis
expr <- read.table("infercnv.observations.txt", header=T) %>% as.matrix()
expr.scale <- scale(t(expr))
tmp1 <- sweep(expr.scale, 2, apply(expr.scale, 2, min),'-')
tmp2 <- apply(expr.scale, 2, max) - apply(expr.scale,2,min)
expr_1 <- t(2*sweep(tmp1, 2, tmp2, "/")-1)
cnv_score <- as.data.frame(colSums(expr_1 * expr_1))
colnames(cnv_score)="cnv_score"
cnv_score <- rownames_to_column(cnv_score, var='cell')

#GO and KEGG analysis
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(enrichplot)
library(data.table)
library(org.Mm.eg.db)
setwd("working directory")
sig_deg<-read.csv("DEGs.csv",header=T,row.name=1)
sig_deg<-sig_deg[sig_deg$p_val_adj<0.05,]
ego_ALL <- enrichGO(gene          = row.names(sig_deg),
                    OrgDb         = ' org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
p1<-dotplot(ego_ALL,split="ONTOLOGY")+facet_grid(ONTOLOGY~.,scales = "free")
genelist <- bitr(row.names(sig_deg), 
                 fromType="SYMBOL",
                 toType="ENTREZID", 
                 OrgDb='org.Hs.eg.db')
genelist <- pull(genelist, ENTREZID)
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
p2<-dotplot(ekegg,showCategory = 10)

#

#GSE154600 analysis
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
pbmc.data1 <- Read10X(data.dir = "Chemoresistance_1 working directory")
pbmc1 <- CreateSeuratObject(counts = pbmc.data1, project = "Chemoresistance_1", 
                            min.cells = 3,min.features = 100,names.delim = "_sc")
pbmc.data2 <- Read10X(data.dir = "Chemoresistance_2 working directory")
pbmc2 <- CreateSeuratObject(counts = pbmc.data2, project = "Chemoresistance_2", 
                            min.cells = 3,min.features = 100,names.delim = "_sc")
pbmc.data3 <- Read10X(data.dir = "Chemosensitivity_3 working directory")
pbmc3 <- CreateSeuratObject(counts = pbmc.data3, project = "Chemosensitivity_3", 
                            min.cells = 3,min.features = 100,names.delim = "_sc")
pbmc.data4 <- Read10X(data.dir = "Chemosensitivity_4 working directory")
pbmc4<- CreateSeuratObject(counts = pbmc.data4, project = "Chemosensitivity_4", 
                           min.cells = 3,min.features = 100,names.delim = "_sc")
pbmc<-merge(pbmc1,c(pbmc2,pbmc3,pbmc4))
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25) 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1500)
top10 <- head(VariableFeatures(pbmc), 10) 
head(pbmc$RNA@var.features,10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimHeatmap(pbmc, dims = 1:18, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:18)
JackStrawPlot(pbmc, dims = 1:18)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:25) 
pbmc <- FindClusters(pbmc, resolution = 0.8)
head(Idents(pbmc), 10)
pbmc <- RunTSNE(pbmc, dims = 1:25)
DimPlot(pbmc,  pt.size = 1,reduction = "tsne")
DimPlot(pbmc, pt.size = 1, reduction = "tsne",label = TRUE)
LabelClusters(DimPlot(pbmc, reduction = "tsne"),id = 'ident')
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

library(SingleR)
hpca.se <- HumanPrimaryCellAtlasData()   #参考数据集
pred.hesc <- SingleR(test =pbmc@assays$RNA@counts, ref = hpca.se, labels = hpca.se$label.main)
new.cluster.ids <- c("Macrophage","T cell","Epithelial cell", "Epithelial cell", "T cell","Fibroblast",
                     "T cell","Macrophage","Epithelial cell","Epithelial cell","Epithelial cell","T cell",
                     "T cell","B cell","Endothelial cell","B cell","Fibroblast","Macrophage",
                     "Macrophage","B cell","Macrophage","Fibroblast","Epithelial cell", "Epithelial cell", 
                     "B cell","T cell","Macrophage", "B cell", "B cell","Fibroblast")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()


