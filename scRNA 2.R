library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
# import the Rdata
load("SRA779509_SRS3805247.sparse.RData")

# split the name in two part, one with the correct one and one with the _ENSG part that i want delete
splitnames = strsplit(rownames(sm), "_ENSG") 

# sapply simplify the split in to a vector, in this case I create a vector with the first part of the 
# all the splits using '[',1
splitnames = sapply(splitnames, `[`, 1)

# I need to replace the _ with - because Seurat give me a warning that say:
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Error in validObject(.Object) : 
#  invalid class “LogMap” object: Duplicate rownames not allowed
splitnames = gsub("_", "-", splitnames)

# Here, I want that all the gene symbols in the splitnames vector are unique because having duplicate 
# row names in your matrix can cause issues when creating the Seurat object
unique_gene_symbols = make.unique(splitnames)

# I add the correct row names and I create the Seurat object
rownames(sm) = unique_gene_symbols
BN_data = CreateSeuratObject(counts = sm)

# Check if the row names are replaced correctly
rownames(BN_data)
head(BN_data@meta.data, 5)
gene_proliferative = c("TUBG1","CDK1","CENPE","RRM2","CCNA2","CDC20","ASPM","NUSAP1","MKI67","TPX2","CCNB2","CCNB1","CLSPN","TOP2A","ZWINT","TYMS","FEN1","PCNA","RRM1","MCM7","CCND3")
BN_data <- subset(BN_data, features = setdiff(rownames(BN_data), gene_proliferative))

BN_data[["percent.mt"]] <- PercentageFeatureSet(BN_data, pattern = "^MT-")
BN_data[["percent.rbp"]] <- PercentageFeatureSet(BN_data, pattern = "^RP[LS]")
VlnPlot(BN_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4)
VlnPlot(BN_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4, pt.size=0, cols = "grey")

plot1 <- FeatureScatter(BN_data, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = "darkgrey") +
  theme(legend.position = "none")
plot2 <- FeatureScatter(BN_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "grey") +
  theme(legend.position = "none")
plot3 <- FeatureScatter(BN_data, feature1 = "nCount_RNA", feature2 = "percent.rbp", cols = "darkgrey") + 
  theme(legend.position = "none")
plot1 + plot2 + plot3

BN_data <- subset(BN_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
BN_data <- NormalizeData(BN_data, normalization.method = "LogNormalize", scale.factor = 10000)

apply(BN_data[["RNA"]]$data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
head(gene.expression, n=50)
VlnPlot(BN_data, features = c("MALAT1","GAPDH"))

CellCycleScoring(BN_data, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> BN_data

BN_data <- FindVariableFeatures(BN_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BN_data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(BN_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(BN_data)
BN_data <- ScaleData(BN_data, features = all.genes)
BN_data@assays$RNA

BN_data <- RunPCA(BN_data, features = VariableFeatures(object = BN_data))

print(BN_data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(BN_data, dims = 1:2, reduction = "pca")
DimPlot(BN_data, reduction = "pca")

elbow_plot <- ElbowPlot(BN_data, ndims = 50)

# Creazione del Data Frame per le linee verticali
vlines_data <- data.frame(
  xintercept = c(10, 16),
  Label = c("PC = 10", "PC = 16")
)

# Aggiunta delle linee verticali con la legenda
elbow_plot <- elbow_plot +
  geom_vline(data = vlines_data, aes(xintercept = xintercept, color = Label), 
             linetype = "dashed", size = 1) +
  scale_color_manual(values = c("PC = 10" = "red", "PC = 16" = "blue")) +
  guides(color = guide_legend(""))

pc.touse <- (BN_data$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.75))
pc.touse

BN_data <- FindNeighbors(BN_data, dims = 1:16)
BN_data <- FindClusters(BN_data, resolution = 0.5)
DimPlot(BN_data, reduction = "pca")
DimPlot(BN_data,reduction="pca", dims=c(4,9))

BN_data <- RunTSNE(BN_data, dims=1:16)
DimPlot(BN_data, reduction = "tsne")

BN_data <- RunUMAP(BN_data, dims = 1:16)
umap_plot =DimPlot(BN_data, reduction = "umap")
cell_counts <- table(Idents(BN_data))
cell_counts_df <- as.data.frame(cell_counts)
names(cell_counts_df) <- c("Cluster", "CellCount")

# Estrai le coordinate UMAP e le identità dei cluster
umap_coords <- Embeddings(BN_data, "umap")
cluster_ids <- Idents(BN_data)

# Combina le coordinate con le identità in un dataframe
umap_df <- data.frame(UMAP_1 = umap_coords[,1], UMAP_2 = umap_coords[,2], Cluster = cluster_ids)

# Calcola i centri dei cluster
cluster_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ Cluster, data = umap_df, FUN = mean)

# Unisci i centri con il numero di cellule
cluster_centers <- merge(cluster_centers, cell_counts_df, by = "Cluster")

# Aggiungi le annotazioni al plot UMAP
umap_plot <- umap_plot +
  geom_text(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = CellCount),
            color = "black", size = 5, hjust = 0.5, vjust = -0.5)


VlnPlot(BN_data,features="nCount_RNA")
VlnPlot(BN_data,features="nFeature_RNA")
VlnPlot(BN_data,features="percent.mt")
VlnPlot(BN_data,features="percent.rbp")

library(ggplot2)
BN_data@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")

cluster2.markers <- FindMarkers(BN_data, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
head(cluster2.markers, n = 5)

cluster2_01.markers <- FindMarkers(BN_data, ident.1 = 2, ident.2 = c(0, 1), min.pct = 0.25)
head(cluster2_01.markers, n = 5)

BN_data.markers <- FindAllMarkers(BN_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BN_data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n = 65)

BN_data.markers %>%
  group_by(cluster) %>%
  slice_min(n = 2, order_by = p_val_adj )

BN_data.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(BN_data, features = top10$gene) + NoLegend()

#################################################################################################################
cluster0.markers <- FindMarkers(BN_data, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster0.markers <- cluster0.markers[order(-cluster0.markers$avg_log2FC),]
head(cluster0.markers, n = 20)
# marker gene primario IL7R, altro gene marker è CD3D. In questo cluster abbiamo Naive CD4+ T cells

#################################################################################################################
cluster1.markers <- FindMarkers(BN_data, ident.1 = 1, min.pct = 0.25, test.use = "wilcox")
cluster1.markers <- cluster1.markers[order(-cluster1.markers$avg_log2FC),]
head(cluster1.markers, n = 20)
# marker gene primario MS4A1 (CD20). In questo cluster abbiamo CD20+ B cells

#################################################################################################################
cluster2.markers <- FindMarkers(BN_data, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
cluster2.markers <- cluster2.markers[order(-cluster2.markers$avg_log2FC),]
head(cluster2.markers, n = 20)
# marker gene primari HBA1, HBA2, HBB. In questo cluster abbiamo Erythroid cells and erythroid like cells

#################################################################################################################
cluster3.markers <- FindMarkers(BN_data, ident.1 = 3, min.pct = 0.25, test.use = "wilcox")
cluster3.markers <- cluster3.markers[order(-cluster3.markers$avg_log2FC),]
head(cluster3.markers, n = 20)
# marker gene primario IL32, S100A10, in questo cluster abbiamo cellule T della memoria centrale

#################################################################################################################
cluster4.markers <- FindMarkers(BN_data, ident.1 = 4, min.pct = 0.25, test.use = "wilcox")
cluster4.markers <- cluster4.markers[order(-cluster4.markers$avg_log2FC),]
head(cluster4.markers, n = 20)
# marker gene primario HBD, In questo cluster abbiamo Erythroid cells and erythroid like cells

#################################################################################################################
cluster5.markers <- FindMarkers(BN_data, ident.1 = 5, min.pct = 0.25, test.use = "wilcox")
cluster5.markers <- cluster5.markers[order(-cluster5.markers$avg_log2FC),]
head(cluster5.markers, n = 20)
# marker gene NGK7, in questo cluster abbiamo cellule NK

##################################################################################################################
cluster6.markers <- FindMarkers(BN_data, ident.1 = 6, min.pct = 0.25, test.use = "wilcox")
cluster6.markers <- cluster6.markers[order(-cluster6.markers$avg_log2FC),]
head(cluster6.markers, n = 20)
# gene marker primario LYZ, in questo cluster abbiamo CD14+ monocytes

##################################################################################################################
cluster7.markers <- FindMarkers(BN_data, ident.1 = 7, min.pct = 0.25, test.use = "wilcox")
cluster7.markers <- cluster7.markers[order(-cluster7.markers$avg_log2FC),]
head(cluster7.markers, n = 20)
# gene marker PRC1, in questo cluster ci sono cellule staminali ematopoietiche

##################################################################################################################
cluster8.markers <- FindMarkers(BN_data, ident.1 = 8, min.pct = 0.25, test.use = "wilcox")
cluster8.markers <- cluster8.markers[order(-cluster8.markers$avg_log2FC),]
head(cluster8.markers, n = 20)
# gene marker CD24, in questo cluster abbiamo B cells CD24+

##################################################################################################################
cluster9.markers <- FindMarkers(BN_data, ident.1 = 9, min.pct = 0.25, test.use = "wilcox")
cluster9.markers <- cluster9.markers[order(-cluster9.markers$avg_log2FC),]
head(cluster9.markers, n = 20)
# gene marker VPREB1, in questo cluster abbiamo precursori delle B cells

##################################################################################################################
cluster10.markers <- FindMarkers(BN_data, ident.1 = 10, min.pct = 0.25, test.use = "wilcox")
cluster10.markers <- cluster10.markers[order(-cluster10.markers$avg_log2FC),]
head(cluster10.markers, n = 20)
# marker gene CD34, in questo cluster abbimo cellule staminali ematopoietiche CD34+

##################################################################################################################
cluster11.markers <- FindMarkers(BN_data, ident.1 = 11, min.pct = 0.25, test.use = "wilcox")
cluster11.markers <- cluster11.markers[order(-cluster11.markers$avg_log2FC),]
head(cluster11.markers, n = 20)
# marker gene  TRDC, in questo cluster ci sono le cellule T gamma delta

##################################################################################################################
cluster12.markers <- FindMarkers(BN_data, ident.1 = 12, min.pct = 0.25, test.use = "wilcox")
cluster12.markers <- cluster12.markers[order(-cluster12.markers$avg_log2FC),]
head(cluster12.markers, n = 20)
# marker gene  IGHG1, in questo cluster abbiamo plasmacellule

##################################################################################################################
cluster13.markers <- FindMarkers(BN_data, ident.1 = 13, min.pct = 0.25, test.use = "wilcox")
cluster13.markers <- cluster13.markers[order(-cluster13.markers$avg_log2FC),]
head(cluster13.markers, n = 20)
# marker gene  LILRA4, in questo cluster abbiamo cellule dentritiche plasmacitoidi

##################################################################################################################

BN_data <- RenameIdents(BN_data, '2' = 'Combined24', '4' = 'Combined24')

DimPlot(BN_data, reduction = "tsne")
DimPlot(BN_data, reduction = "umap")
FeaturePlot(BN_data, features = c("CCR7", "MS4A1", "HBA1", "IL7R","HBD", "NKG7", "LYZ", "PRC1", "CD24", "VPREB1", "CD34", "TRDC", "IGHG1", "LILRA4"))
DotPlot(BN_data, features = c("CCR7", "MS4A1", "HBA1", "IL7R", "NKG7", "LYZ", "PRC1", "CD24", "VPREB1", "CD34", "TRDC", "IGHG1", "LILRA4"))

new.cluster.ids <- c("Erythroid-like and Erythroid precursor cells","Naive CD4+ T cells","CD20+ B cells","Central memory T cells","Natural Killer cells", "CD14+ Monocytes","HSC","CD24+ B cells","Precursor B cells","CD34+ HSC","Gamma delta T cells","Plasma cells", "Plasmacytoid dentritic cells")
names(new.cluster.ids) <- levels(BN_data)
BN_data <- RenameIdents(BN_data, new.cluster.ids)
DimPlot(BN_data, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(BN_data, reduction = "tsne", label = TRUE, pt.size = 0.5)
