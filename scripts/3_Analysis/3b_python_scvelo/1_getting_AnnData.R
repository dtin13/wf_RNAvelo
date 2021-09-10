## FROM SEURAT(R) TO ANNDATA(PYTHON)

m.RNA <-GetAssayData(object = seurat, assay = ’RNA’, slot = "counts")
m.spliced <- GetAssayData(object = seurat, assay = ’spliced’, slot = "counts")
m.unspliced <-GetAssayData(object = seurat, assay = ’unspliced’, slot = "counts")
m.spliced <- m.spliced[rownames(seurat),]
m.unspliced <- m.unspliced[rownames(seurat),]

writeMM(t(m.spliced),’splicedMM.mtx’)
writeMM(t(m.unspliced),’unsplicedMM.mtx’)
writeMM(t(m.RNA),’RNAMM.mtx’)

write.csv(rownames(seurat), file = "genenames.csv", row.names = FALSE, quote = F)

scvelo_seurat.obs <- as.data.frame(cbind(seurat@meta.data$Sample,

seurat@meta.data$RNA_snn_res.0.2))

colnames(scvelo_seurat.obs) <- c("orig.ident", "clusters")
scvelo_seurat.obsm <- as.data.frame(Embeddings(seurat, reduction = "umap"))

write.csv(scvelo_seurat.obs, file = "seurat.obs.csv", row.names = FALSE)
write.csv(scvelo_seurat.obsm, file = "seurat.obsm.csv", row.names = FALSE)
