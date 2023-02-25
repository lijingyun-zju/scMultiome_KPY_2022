# R 4.1.1 required
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(argparser, quietly=TRUE)

### first remove the wierd cluster, to enable better separation of the main clusters

# parameters to use
nfeat = 2000 # number of variable feature
nPC = 75 # number of PCs

print(nfeat)
print(nPC)
set.seed(1)

outdir = paste0("Output.nfeat",nfeat,".nPC",nPC,'.rm_cellcycle')
system(paste0("mkdir -p ",outdir))
Kras = readRDS('Kras_P53.2Kgene.RNA.KPY_only.rds')

# remove the outlier 
Kras <- subset(Kras, subset = RNA_snn_res.0.5 != 7)

# re-start processing

# normalization
Kras <- NormalizeData(Kras, normalization.method = "LogNormalize", scale.factor = 10000)

# get cell cycle genes
cc = read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Reference/mouse/cell.cycle/Mus_musculus.parse.csv",head=T)
g2m_genes = as.character(cc[cc$phase == "G2M",]$genename)
s_genes = as.character(cc[cc$phase == "S",]$genename)

Kras <- CellCycleScoring(Kras,
                                g2m.features = g2m_genes,
                                 s.features = s_genes)
# find variation genes
Kras <- FindVariableFeatures(Kras, selection.method = "vst", nfeatures = nfeat)    # use 2k here, can be optimizedd?

# scaling
all.genes <- rownames(Kras)
Kras <- ScaleData(Kras, features = all.genes)

Kras <- RunPCA(Kras,npcs=100, features = VariableFeatures(object = Kras))
pdf(paste0(outdir,"/PCA.by.CellCycle.before_regress.pdf"))
DimPlot(Kras,
        reduction = "pca",
        group.by= "Phase")
dev.off()

# run PCA
## Regress out cell cycle scores during data scaling
Kras <- ScaleData(Kras, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Kras))
Kras <- RunPCA(Kras,npcs=100, features = VariableFeatures(object = Kras))
pdf(paste0(outdir,"/PCA.by.CellCycle.after_regress.pdf"))
DimPlot(Kras,
        reduction = "pca",
        group.by= "Phase")
dev.off()

# graph-based clustering
Kras <- FindNeighbors(Kras, dims = 1:nPC)
Kras <- FindClusters(Kras, resolution = 1)

# UMAP/tSNE visualization
Kras <- RunUMAP(Kras, dims = 1:nPC)
Kras <- RunTSNE(Kras, dims = 1:nPC, tsne.method = "Rtsne", nthreads = 4, max_iter = 2000,reduction.name = "Rtsne")

pdf(paste0(outdir,'/Kras_P53.2D_reduction.KPY_only.pdf'))
DimPlot(Kras, reduction = "umap",shuffle=T)
DimPlot(Kras, reduction = "Rtsne")
dev.off()

# individally show each sample
meta <- Kras@meta.data

pdf(paste0(outdir,"/CellType.marker.gene.plot.Kras_P53.KPY_only.pdf"))
DimPlot(Kras,group.by = "seurat_clusters",reduction = "umap",shuffle=T,label=T)
FeaturePlot(Kras, features = c("Nkx2-1","Sox9","Ly6a","Hmga2","Cd24a","Epcam","Cd44"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Pdpn","Hopx","Cav1","Aqp5","Vegfa","Spock2","Lmo7"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Sftpa1","Sftpc","Soat1","Lpcat1","Etv5","Abca3","Cebpa","Acly"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Cpm","Nr4a1","Hmgcr","Hmgcs1","Lyz2","Sftpb","Sftpd"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Scgba3a2","H2-K1","Scgb1a1","Ccsp","Cc10"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Krt5"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("Foxj1"),cols = c("lightgrey","firebrick3"))
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

pdf(paste0(outdir,"/KP.organoid.gene.scores.KPY_only.pdf"))
FeaturePlot(Kras, features = c("Kras.score","AT2.score","Prolif.score"),cols = c("lightgrey","firebrick3"))
FeaturePlot(Kras, features = c("G12C.induced.score","G12C.suppress.score","KRAS.inhib.Resist.score","KRAS.inhib.Sensit.score"),cols = c("lightgrey","firebrick3"))
dev.off()

## plot by clsuter; using violin plot
library(ggplot2)
mytheme10 <- theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1,axis.text = element_text(size=rel(1.2)),axis.title = element_text(size=rel(1.2)),legend.text=element_text(size=12),plot.title=element_text(size=rel(1.2)),legend.title =element_text(size=rel(1.2)))

Marks = FetchData(Kras,c('Cd44','Procr','F3'))
meta <- Kras@meta.data
all(rownames(Marks) == rownames(meta))
meta[,c('Cd44','Procr','F3')] = Marks[,c('Cd44','Procr','F3')]

colname.list = colnames(meta)
pdf(paste0(outdir,"/Aggregated.scores.by.clusters.KPY_only.pdf"))
for(col in c("Kras.score","AT2.score","Prolif.score","G12C.induced.score","G12C.suppress.score","KRAS.inhib.Resist.score","KRAS.inhib.Sensit.score",'Cd44','Procr','F3')){
 p = ggplot(meta,aes(x=seurat_clusters,y=meta[,col])) + geom_violin() + mytheme10 + ylab(as.character(col)) + ggtitle(as.character(col))
 print(p)
 p = ggplot(meta,aes(x=seurat_clusters,y=meta[,col])) + geom_boxplot() + mytheme10 + ylab(as.character(col)) + ggtitle(as.character(col))
 print(p)
}
dev.off()

meta = Kras@meta.data
umapCoord <- as.data.frame(Embeddings(object = Kras[["umap"]]))
all(rownames(umapCoord) == rownames(meta))
meta$UMAP1.rna = umapCoord$UMAP_1
meta$UMAP2.rna = umapCoord$UMAP_2

saveRDS(meta,paste0('Kras_P53.2Kgene.RNA.KPY_only.',outdir,'.meta.rds'))
saveRDS(Kras,paste0('Kras_P53.2Kgene.RNA.KPY_only.',outdir,'.rds'))




