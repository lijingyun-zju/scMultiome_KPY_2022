# R 4.1.1 required
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)

set.seed(1)

# load RNA and ATAC data
data.KPY = Read10X_h5('./Data/KPY/filtered_feature_bc_matrix.h5')

frag.KPY <- "./Data/KPY/atac_fragments.tsv.gz"

#### processing RNA-seq
RNA.KPY = data.KPY$`Gene Expression`

# first assign sample name to column names
colnames(RNA.KPY) = paste0("KPY:",colnames(RNA.KPY))

# get gene annotations
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

RNA.KPY.mat = as.matrix(RNA.KPY)
# Save for DEG testing
saveRDS(RNA.KPY.mat,'KPY.raw.mat.rds')






