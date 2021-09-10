## PREVIOUS STEPS

# 1. Install packages
if (!require(devtools)) {
install.packages("devtools")
}
devtools::install_github("BUStools/BUSpaRse")
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github("velocyto-team/velocyto.R")
if (!require(BiocManager)) {
install.packages("BiocManager")
}
BiocManager::install(c("DropletUtils", "BSgenome.Mmusculus.UCSC.mm10",
"AnnotationHub", "SingleR"))

library(BUSpaRse)
library(Seurat)
library(SeuratWrappers)
library(BSgenome.Mmusculus.UCSC.mm10)
library(AnnotationHub)
library(zeallot)
library(DropletUtils)
library(tidyverse)
library(GGally)
library(velocyto.R)
library(SingleR)
library(scales)
library(Matrix)

# 2. Reference genome preparation to transcript detection
ah <- AnnotationHub()
query(ah, pattern = c("Ensembl", "97", "Mus musculus", "EnsDb"))
edb <- ah[["AH73905"]] #Anotación Ensembl 97
get_velocity_files(edb, L = 91, Genome = BSgenome.Mmusculus.UCSC.mm10,
out_path = "./output",
isoform_action = "separate")
kallisto index -i ./output/mm_cDNA_introns_97.idx ./output/cDNA_introns.fa
