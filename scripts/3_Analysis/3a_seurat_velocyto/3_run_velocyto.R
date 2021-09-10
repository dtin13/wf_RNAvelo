#Velocyto

seurat <- RunVelocity(object = seurat, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = seurat)))
names(x = ident.colors) <- levels(x = seurat)
cell.colors <- ident.colors[Idents(object = seurat)]
names(x = cell.colors) <- colnames(x = seurat)
velo <- show.velocity.on.embedding.cor(emb = Embeddings(object = seurat, reduction = "umap"), vel = Tool(object = seurat, 
                                                                                             slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
saveRDS(seurat, "runvelo.rds")
