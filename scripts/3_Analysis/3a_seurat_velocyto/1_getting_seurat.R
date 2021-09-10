## CREATING SEURAT OBJECTS TO POSTERIOR BIONF ANALYSIS

# 1. Spliced and unspliced matrix generation 
t <- readr::read_delim(’matrix.mtx’, delim = ’ ’, skip = 3, col_names = FALSE)
spliced <- Matrix::sparseMatrix(
i = t$X1,
j = t$X2,
x = t$X3,
dims = c(55487, 6794880)
)
Matrix::writeMM(spliced, ’spliced.mtx’)
unspliced <- Matrix::sparseMatrix(
i = t$X1,
j = t$X2,
x = t$X4,
dims = c(55487, 6794880)
)
Matrix::writeMM(unspliced, ’unspliced.mtx’)

# 2. Creates Seurat objects
seurat <- Read10X(data.dir = "~/s_97/Solo.out/Gene/filtered", gene.column = 2)
seurat <- CreateSeuratObject(counts = seurat, min.cells = 3, min.features = 200, project = "Single_cell_skin")
spliced <- Read10X(data.dir = "~/s_97/Solo.out/Velocyto/raw/spliced", gene.column = 2)
unspliced<- Read10X(data.dir = "~/s_97/Solo.out/Velocyto/raw/unspliced", gene.column = 2)
spliced <- spliced[rownames(spliced) %in% rownames(seurat),
colnames(spliced) %in% colnames(seurat)]
unspliced <- unspliced[rownames(unspliced) %in% rownames(seurat),
colnames(unspliced) %in% colnames(seurat)]
seurat[[’spliced’]] <- CreateAssayObject(spliced)
seurat[[’unspliced’]] <- CreateAssayObject(unspliced)
seurat$Sample <- factor(c(’5w1’))

#3. (optional) If work with more than a single sample, unificate them using 'merge()'
