# R 4.1.1 required
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)

set.seed(1)

# load RNA and ATAC data
data.KPY = Read10X_h5('./Data/KPY/filtered_feature_bc_matrix.h5')
data.YFP = Read10X_h5('./Data/YFP/filtered_feature_bc_matrix.h5')

frag.KPY <- "./Data/KPY/atac_fragments.tsv.gz"
frag.YFP <- "./Data/YFP/atac_fragments.tsv.gz"

#### processing RNA-seq
RNA.KPY = data.KPY$`Gene Expression`
RNA.YFP = data.YFP$`Gene Expression`

# first assign sample name to column names
colnames(RNA.KPY) = paste0("KPY:",colnames(RNA.KPY))
colnames(RNA.YFP) = paste0("YFP:",colnames(RNA.YFP))

data.merge <- cbind(RNA.KPY,RNA.YFP)

# get gene annotations
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
Kras <- CreateSeuratObject(
  counts = data.merge,
  assay = "RNA"
)

# add sample name info
library(stringr)
samp <- str_split_fixed(rownames((Kras@meta.data)),":",2)
Kras$lib = samp[,1]

# QC
Kras[["percent.mt"]] <- PercentageFeatureSet(Kras, pattern = "^mt-")

pdf("MultiOmic.RNA.QC.metric.20220430.pdf",useDingbats=F)
VlnPlot(Kras, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,group.by="lib" ,flip =TRUE,pt.size=0)
VlnPlot(Kras, features = c("percent.mt"), ncol = 1,group.by="lib" ,flip =TRUE,pt.size=0)

plot1 <- FeatureScatter(Kras, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by="lib")
plot2 <- FeatureScatter(Kras, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="lib")
plot1
plot2
dev.off()


# filter cells with low quality, requring num of genes > 200, num of reads > 2500, Mitochondrial% < 10%
# mt% is set to 10%, based on the mt vs nCount
Kras <- subset(Kras, subset = nFeature_RNA > 200 & nFeature_RNA > 2500 & percent.mt < 15)

# Start processing

# normalization
Kras <- NormalizeData(Kras, normalization.method = "LogNormalize", scale.factor = 10000)

# find variation genes
Kras <- FindVariableFeatures(Kras, selection.method = "vst", nfeatures = 2000)    # use 2k here, can be optimizedd?
top10 <- head(VariableFeatures(Kras), 10)
top10_plusKras = c(top10,"Kras")
pdf("Kras_P53.variable.gene.pdf")
plot1 <- VariableFeaturePlot(Kras)
plot2 <- LabelPoints(plot = plot1, points = top10_plusKras, repel = TRUE)
plot1
plot2
dev.off()

# scaling
all.genes <- rownames(Kras)
Kras <- ScaleData(Kras, features = all.genes)

# run PCA
Kras <- RunPCA(Kras,npcs=100, features = VariableFeatures(object = Kras))
Kras <- JackStraw(Kras, num.replicate = 100,dims = 100)
Kras <- ScoreJackStraw(Kras, dims = 1:100)

pdf("Kras_P53.PCs.pdf",width=12)
#JackStrawPlot(Kras, dims = 1:75)
ElbowPlot(Kras,ndims=100)

DimHeatmap(Kras, dims = c(1:3, 70:75), cells = 500, balanced = TRUE)
dev.off()

#decide to use top 60 PCs
# graph-based clustering
Kras <- FindNeighbors(Kras, dims = 1:60)
Kras <- FindClusters(Kras, resolution = 0.5)

# UMAP/tSNE visualization
Kras <- RunUMAP(Kras, dims = 1:60)
Kras <- RunTSNE(Kras, dims = 1:60, tsne.method = "Rtsne", nthreads = 4, max_iter = 2000,reduction.name = "Rtsne")

pdf("Kras_P53.2D_reduction.pdf")
DimPlot(Kras, reduction = "umap",shuffle=T)
DimPlot(Kras,group.by="lib",reduction = "umap",shuffle=T,label=T)
DimPlot(Kras, reduction = "pca",group.by="lib",shuffle=T,label=T)
DimPlot(Kras, reduction = "Rtsne")
dev.off()

# individally show each sample
meta <- Kras@meta.data

p1 <- DimPlot(Kras,group.by="lib",reduction = "umap",shuffle=T,label=T,cells=rownames(meta[meta$lib=="KPY",]))
p2 <- DimPlot(Kras,group.by="lib",reduction = "umap",shuffle=T,label=T,cells=rownames(meta[meta$lib=="YFP",]))

p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)

pdf("Kras_P53.2D_reduction.bysamp.pdf",width=14)
p1 + p2
dev.off()

pdf("CellType.marker.gene.plot.Kras_P53.pdf")
DimPlot(Kras,group.by = "seurat_clusters",reduction = "umap",shuffle=T,label=T)
FeaturePlot(Kras, features = c("Nkx2-1","Sox9","Ly6a","Hmga2","Cd24a","Epcam"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Pdpn","Hopx","Cav1","Aqp5","Vegfa","Spock2","Lmo7"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Sftpa1","Sftpc","Soat1","Lpcat1","Etv5","Abca3","Cebpa","Acly"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Cpm","Nr4a1","Hmgcr","Hmgcs1","Lyz2","Sftpb","Sftpd"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Scgba3a2","H2-K1","Scgb1a1","Ccsp","Cc10"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Krt5"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Foxj1"),cols = c("lightgrey","firebrick3"))
dev.off()

library(Hmisc)
pdf("CellType.marker.gene.plot.Kras_P53.gastric.pdf",width=10,height=12)
gastric = capitalize(c('cdx2','hnf4a','pdx1','muc2','muc6','tff2','prss1','cps1','spdef','muc5ac'))
FeaturePlot(Kras,features=gastric,cols = c("lightgrey","firebrick3"),ncol=3)
FeaturePlot(Kras,features=gastric,cols = c("lightgrey","firebrick3"),ncol=3,order=T)
dev.off()


## compiled gene score
# read in gene lists
Kras.gene <- read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/Kras.gene.mouse.txt",head=T)
AT2.gene <- read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/AT2marker.gene.mouse.txt",head=T)
Prolif.gene <- read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/Proliferate.gene.txt",head=F)
colnames(Prolif.gene) = c("mouse")
G12C.induced <- read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/G12C.induced.mouse.txt",head=T)
G12C.suppress <- read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/G12C.suppress.mouse.txt",head=T)
KRAS.inhib.Resist <- read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/KRAS.inhib.Resist.mouse.txt",head=T)
KRAS.inhib.Sensit <- read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/KRAS.inhib.Sensit.mouse.txt",head=T)

Kras.gene.list <- unique(Kras.gene$mouse)

AT2.gene.list <- unique(AT2.gene$mouse)
Prolif.gene.list <- unique(Prolif.gene$mouse)
G12C.induced.list <- unique(G12C.induced$mouse)
G12C.suppress.list <- unique(G12C.suppress$mouse)
KRAS.inhib.Resist.list <- unique(KRAS.inhib.Resist$mouse)
KRAS.inhib.Sensit.list <- unique(KRAS.inhib.Sensit$mouse)

Kras.score = FetchData(Kras,Kras.gene.list)
AT2.score = FetchData(Kras,AT2.gene.list)
Prolif.score = FetchData(Kras,Prolif.gene.list)
G12C.induced.score = FetchData(Kras,G12C.induced.list)
G12C.suppress.score = FetchData(Kras,G12C.suppress.list)
KRAS.inhib.Resist.score = FetchData(Kras,KRAS.inhib.Resist.list)
KRAS.inhib.Sensit.score = FetchData(Kras,KRAS.inhib.Sensit.list)

meta <- Kras@meta.data
meta$Kras.score = rowMeans(Kras.score)
meta$AT2.score = rowMeans(AT2.score)
meta$Prolif.score = rowMeans(Prolif.score)
meta$G12C.induced.score = rowMeans(G12C.induced.score)
meta$G12C.suppress.score = rowMeans(G12C.suppress.score)
meta$KRAS.inhib.Resist.score = rowMeans(KRAS.inhib.Resist.score)
meta$KRAS.inhib.Sensit.score = rowMeans(KRAS.inhib.Sensit.score)

Kras$Kras.score = meta$Kras.score
Kras$AT2.score = meta$AT2.score
Kras$Prolif.score = meta$Prolif.score
Kras$G12C.induced.score = meta$G12C.induced.score
Kras$G12C.suppress.score = meta$G12C.suppress.score
Kras$KRAS.inhib.Resist.score = meta$KRAS.inhib.Resist.score
Kras$KRAS.inhib.Sensit.score = meta$KRAS.inhib.Sensit.score

pdf("KP.organoid.gene.scores.pdf")
FeaturePlot(Kras, features = c("Kras.score","AT2.score","Prolif.score"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("G12C.induced.score","G12C.suppress.score","KRAS.inhib.Resist.score","KRAS.inhib.Sensit.score"),cols = c("lightgrey","firebrick3"))
dev.off()

meta <- Kras@meta.data

## plot by clsuter; using violin plot
library(ggplot2)
mytheme10 <- theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1,axis.text = element_text(size=rel(1.2)),axis.title = element_text(size=rel(1.2)),legend.text=element_text(size=12),plot.title=element_text(size=rel(1.2)),legend.title =element_text(size=rel(1.2)))

colname.list = colnames(meta)
pdf("Aggregated.scores.by.clusters.pdf")
for(col in 8:14){
 p = ggplot(meta,aes(x=seurat_clusters,y=meta[,col])) + geom_violin() + mytheme10 + ylab(as.character(colname.list[col])) + ggtitle(as.character(colname.list[col]))
 print(p)
 p = ggplot(meta,aes(x=seurat_clusters,y=meta[,col])) + geom_boxplot() + mytheme10 + ylab(as.character(colname.list[col])) + ggtitle(as.character(colname.list[col]))
 print(p)
}
dev.off()

saveRDS(Kras,"Kras_P53.2Kgene.RNA.rds")







