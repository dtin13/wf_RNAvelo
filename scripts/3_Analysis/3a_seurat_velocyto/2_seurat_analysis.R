library(ggplot2)
library(dplyr)
library(sctransform)
library(data.table)
library(SoupX)
library(DropletUtils)
library(clustree)
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(RColorBrewer)

##. PreQC
# Mitochondrial genes - check levels of expression for mt genes. 
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
# Ribosomal genes - check levels of expression for rb genes.
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL][[:digit:]]")
#If we obtain NaN values...
seurat$percent.mt[is.nan(seurat$percent.mt)] <- 0
seurat$percent.ribo[is.nan(seurat$percent.ribo)] <- 0

# QC: violin plots.
p1 <- VlnPlot(seurat, features = c("nFeature_spliced"), pt.size = 0.25) + ggtitle("Nº features") + theme(legend.position="bottom") 
p2 <- VlnPlot(seurat, features = c("nCount_spliced"), pt.size = 0.25)  + ggtitle("Nº counts") + theme(legend.position="bottom")
p3 <- VlnPlot(seurat, features = c("percent.mt"), pt.size = 0.25) + ggtitle("Mitochondrial %") + theme(legend.position="bottom")
p4 <- VlnPlot(seurat, features = c("percent.ribo"), pt.size = 0.25) + ggtitle("Ribosomal %") + theme(legend.position="bottom")

##. PostQC
#Pre-filtering stats calculus.
stats_name <- c("N_cells","Count_median","Feature_median","Mit_median","Ribo_median")
stats_pre <- c(length(colnames(seurat)), median(seurat@meta.data[["nCount_spliced"]]), median(seurat@meta.data[["nFeature_spliced"]]), median(seurat@meta.data[["percent.mt"]]), median(seurat@meta.data[["percent.ribo"]]))
# We should apply the filterings once the QC plots (GenePlot and Violin plots) have been checked.
# Feature filter (select your own numbers)
cells_seurat <- FetchData(object = seurat, vars = "nFeature_spliced")     
seurat <- seurat[, which(x = cells_seurat > 250 & cells_seurat < 4000)]                          
# Mitochondrial filter (select your own numbers)
mit_seurat <- FetchData(object = seurat, vars = "percent.mt")
seurat <- seurat[, which(x = mit_seurat < 3.5)]
# QC: violin plots - After filter.
pp1 <- VlnPlot(seurat, features = c("nFeature_spliced"), pt.size = 0.25) + ggtitle("Nº features") + theme(legend.position="bottom") 
pp3 <- VlnPlot(seurat, features = c("percent.mt"), pt.size = 0.25) + ggtitle("Mitochondrial %") + theme(legend.position="bottom")
# Post-filter stats calculus.
stats_post <- c(length(colnames(seurat)), median(seurat@meta.data[["nCount_spliced"]]), median(seurat@meta.data[["nFeature_spliced"]]), median(seurat@meta.data[["percent.mt"]]), median(seurat@meta.data[["percent.ribo"]]))

stats_gral<-rbind(stats_name, stats_pre, stats_post)
write.table(stats_gral, "qc_stats", col.names = FALSE)
# Save RDS: we can use this object to generate all the rest of the data.
saveRDS(seurat, file = "seurat_postqc.rds")

##. Normalization, Scaling and HVG identification.
seurat <- SCTransform(seurat, assay='spliced', new.assay.name='RNA',verbose=FALSE)

##. DOWNSTREAM ANALYSIS.
# Perform linear dimensional reduction , by default only the variable features are taken as input.
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 50) # This result could all be saved in a table.

# Determine the dimensionality of the dataset
pdf('elbowplot.pdf')
ElbowPlot(seurat, ndims = 50) + theme(legend.position="bottom") 
dev.off()

# Run non-linear dimensional reduction (UMAP/tSNE) [We decide to choose (e.g.) 15 PCs based on elbow plot]
seurat <- RunUMAP(object = seurat, dims = 1:15)
seurat <- RunTSNE(object = seurat, dims = 1:15)

# Construction of KNN graph based on the euclidean distance in PCA space (k=20 by default).Refine edge weights with Jaccard similarity.
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:15)

# Cluster cells with Louvain algorithm (by default). The resolution parameter sets the ‘granularity’ of the downstream clustering. We perform clustering with different resolution values to test. cluster robustness
seurat <- FindClusters(object = seurat, resolution = c(0.2, 0.4, 0.6, 0.7, 0.8, 1, 1.2))

clustree_plot <- clustree(seurat)


DefaultAssay(object = seurat) <- "RNA"
Idents(seurat) <- seurat$RNA_snn_res.0.6 

pdf('tsne_seurat.pdf')
DimPlot(object = seurat, reduction = "tsne", label = TRUE)
dev.off()

pdf('umap_seurat.pdf')
DimPlot(object = seurat, reduction = "umap", label = TRUE)
dev.off()


saveRDS(seurat, file = "seurat_clusters.rds")
# Find markers for every cluster compared to all remaining cells, report only the positive ones

seurat.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(seurat.markers, file = 'seurat.markers.txt', col.names = TRUE, row.names = TRUE, sep = '\t')

top10seurat <- seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

heatmap <- DoHeatmap(object = seurat, features = top10seurat$gene, slot = 'scale.data') + NoLegend() + theme (axis.text.y = element_text(size=2.5)) + scale_fill_gradientn(colors = c("blue","black","red"))

# Subset (select clusters of interest)
seurat<-subset(seurat, idents =c("1", "2", "3","4","5"))
levels(seurat)
saveRDS(seurat, "seurat_subset.rds")
seurat <- SCTransform(seurat, assay='spliced', new.assay.name='RNA',verbose=FALSE)
saveRDS(seurat, "subset_normalized.rds")
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = 50)
pdf('elbowplot.pdf')
ElbowPlot(seurat, ndims = 50) + theme(legend.position="bottom") 
dev.off()
seurat <- RunUMAP(object = seurat, dims = 1:20)
seurat <- RunTSNE(object = seurat, dims = 1:20)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:20)
seurat <- FindClusters(object = seurat, resolution = c(0.2, 0.4, 0.6, 0.7, 0.8, 1, 1.2))
clustree_plot <- clustree(seurat)
DefaultAssay(object = seurat) <- "RNA"
Idents(seurat) <- seurat$RNA_snn_res.0.2
pdf('tsne_seurat.pdf')
DimPlot(object = seurat, reduction = "tsne", label = TRUE)
dev.off()

pdf('umap_seurat.pdf')
DimPlot(object = seurat, reduction = "umap", label = TRUE)
dev.off()

saveRDS(seurat, file = "subset_clusters.rds")
subset.markers <- FindAllMarkers(object = seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
subset.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(subset.markers, file = 'subset.markers.txt', col.names = TRUE, row.names = TRUE, sep = '\t')
top10subset <- subset.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
heatmap <- DoHeatmap(object = seurat, features = top10subset$gene, slot = 'scale.data') + NoLegend() + theme (axis.text.y = element_text(size=2.5)) + scale_fill_gradientn(colors = c("blue","black","red"))

