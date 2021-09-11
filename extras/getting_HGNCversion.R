## Modifying genenames (HUMAN): From ENSEMBL to HGNC version

library(biomaRt)
library(dplyr)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes <- read.table("~/genes.txt", quote="\"", comment.char="")
genes$V2 <- sub("[.]$", "", genes$V1)
colnames(genes) <- c("original","ensembl_gene_id_version")
w <- getBM(c("ensembl_gene_id_version", "hgnc_symbol", "clone_based_ensembl_gene"), "ensembl_gene_id_version", genes$ensembl_gene_id_version, mart, useCache = FALSE)
w$symbol <- coalesce(w$hgnc_symbol, w$clone_based_ensembl_gene)
genes$symbol <- w[match(genes$ensembl_gene_id_version, w$ensembl_gene_id_version), 4]
genes$symbol <- toupper(genes$symbol)
write.table(genes,file = "NEWgenes.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

mito.genes <- c("ENSG00000210049.1.","ENSG00000211459.2.","ENSG00000210077.1.","ENSG00000210082.2.","ENSG00000209082.1.","ENSG00000198888.2.","ENSG00000210100.1.","ENSG00000210107.1.","ENSG00000210112.1.","ENSG00000198763.3.","ENSG00000210117.1.","ENSG00000210127.1.","ENSG00000210135.1.","ENSG00000210140.1.","ENSG00000210144.1.","ENSG00000198804.2.","ENSG00000210151.2.","ENSG00000210154.1.","ENSG00000198712.1.","ENSG00000210156.1.","ENSG00000228253.1.","ENSG00000198899.2.","ENSG00000198938.2.","ENSG00000210164.1.","ENSG00000198840.2.","ENSG00000210174.1.","ENSG00000212907.2.","ENSG00000198886.2.","ENSG00000210176.1.","ENSG00000210184.1.","ENSG00000210191.1.","ENSG00000198786.2.","ENSG00000198695.2.","ENSG00000210194.1.","ENSG00000198727.2.","ENSG00000210195.2.","ENSG00000210196.2.")
mito.genes2 <- intersect(mito.genes, rownames(seurat))
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, features = mito.genes2)

#MODIFY FINAL SEURAT.MARKERS TABLE
NEWgenes <- read.delim("~/Escritorio/scBladder/kb_human/kb_analysis/NEWgenes.txt", header=FALSE)
seurat.markers$HGNC_symbol <- NEWgenes[match(seurat.markers$gene, NEWgenes$V1), 3]
