# R 4.1.1 required
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(argparser, quietly=TRUE)

Kras = readRDS('Kras_P53.2Kgene.RNA.KPY_only.Output.nfeat2000.nPC75.rm_cellcycle.rds')

meta = Kras@meta.data
meta$Population1[meta$RNA_snn_res.1 %in% c(5)] = 'Hmga2_A'
meta$Population1[meta$RNA_snn_res.1 %in% c(3,8,1,0,10)] = 'Hmga2_B'
meta$Population1[meta$RNA_snn_res.1 %in% c(2,9,6,7)] = 'Nkx2-1_B'
meta$Population1[meta$RNA_snn_res.1 %in% c(4)] = 'Nkx2-1_A'

meta$Population2[meta$Population1 %in% c('Hmga2_A','Hmga2_B')] = 'Hmga2'
meta$Population2[meta$Population1 %in% c('Nkx2-1_A','Nkx2-1_B')] = 'Nkx2-1'

Kras$Population1 = meta$Population1
Kras$Population2 = meta$Population2

### cancer state genes from cancer cell paper
state.gene <- read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/TableS2.CancaerCell.txt",head=T,sep="\t")
state.name <- colnames(state.gene)
meta <- Kras@meta.data

for(i in 2:ncol(state.gene)){
  df = state.gene[order(state.gene[,i],decreasing=T),]
  df.genes = as.character(df[1:100,1])
  df.score = FetchData(Kras,df.genes)
  meta[,state.name[i]] = rowMeans(df.score)
}

meta.cut = meta[,c(paste0("NMF_0",1:9),"NMF_10","NMF_11")]
Kras <- AddMetaData(Kras,meta.cut )  # add column data

pdf("Kras.cancerState.gene.scores.tumor_cls.20220518.pdf")
for(i in 2:ncol(state.gene)){
  p = FeaturePlot(Kras, features = state.name[i],cols = c("dodgerblue3","firebrick1"))
  print(p)
}
dev.off()

library(ggplot2)
mytheme10 <- theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1,axis.text = element_text(size=rel(1.2)),axis.title = element_text(size=rel(1.2)),legend.text=element_text(size=12),plot.title=element_text(size=rel(1.2)),legend.title =element_text(size=rel(1.2)))

pdf("Aggregated.cancerState.gene.scores.tumor_cls.boxplot.20220518.pdf",useDingbats=F)
for(state in state.name){
 p = ggplot(meta,aes(x=seurat_clusters,y=meta[,state])) + geom_violin() + mytheme10 + ylab(state) + ggtitle(as.character(state))
 print(p)
 p = ggplot(meta,aes(x=seurat_clusters,y=meta[,state])) + geom_boxplot() + mytheme10 + ylab(state) + ggtitle(as.character(state))
  print(p)
}
dev.off()

# calculate the median of expression across cells in clusters
state.exp.mid = aggregate(meta[,c(paste0("NMF_0",1:9),"NMF_10","NMF_11")], list(meta$seurat_clusters), FUN=median)

library(pheatmap)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
# name of cancer scores
cvt.name = read.table('/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/TableS2.CancaerCell.names.txt',head=F,sep='\t')
colnames(state.exp.mid) = cvt.name$V2
pdf("Kras.cancerState.gene.scores.tumor_cls.20220518.heatmap.pdf")# this one not used; look at the local R studio
pheatmap(t(state.exp.mid),scale = "none",cluster_row = F,cluster_cols = T,main="cancer state no_norm")
pheatmap(t(state.exp.mid),scale = "row",cluster_row = F,cluster_cols = T,main="cancer state by row")
saveRDS(state.exp.mid,'Aggregated.cancerState.gene.scores.tumor_cls.state.exp.mid.20220518.rds')

# calculate the median of expression across cells in populations
state.exp.mid = aggregate(meta[,c(paste0("NMF_0",1:9),"NMF_10","NMF_11")], list(meta$Population1), FUN=median)

library(pheatmap)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
colnames(state.exp.mid) = cvt.name$V2
pheatmap(t(state.exp.mid),scale = "none",cluster_row = F,cluster_cols = T,main="cancer state no_norm")
pheatmap(t(state.exp.mid),scale = "row",cluster_row = F,cluster_cols = T,main="cancer state by row")
dev.off()
saveRDS(state.exp.mid,'Aggregated.cancerState.gene.scores.tumor_cls.state.exp.mid_cross_population1.20220518.rds')

# calculate the mean of expression across cells in clusters
state.exp.mid = aggregate(meta[,c(paste0("NMF_0",1:9),"NMF_10","NMF_11")], list(meta$seurat_clusters), FUN=mean)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
# name of cancer scores
cvt.name = read.table('/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/TableS2.CancaerCell.names.txt',head=F,sep='\t')
colnames(state.exp.mid) = cvt.name$V2
saveRDS(state.exp.mid,'Aggregated.cancerState.gene.scores.tumor_cls.state.exp.mean.20220518.rds')

# calculate the mean of expression across cells in populations
state.exp.mid = aggregate(meta[,c(paste0("NMF_0",1:9),"NMF_10","NMF_11")], list(meta$Population1), FUN=mean)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
# name of cancer scores
cvt.name = read.table('/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/TableS2.CancaerCell.names.txt',head=F,sep='\t')
colnames(state.exp.mid) = cvt.name$V2
saveRDS(state.exp.mid,'Aggregated.cancerState.gene.scores.tumor_cls.state.exp.mean_cross_population1.20220518.rds')

### using top 10 genes
### cancer state genes from cancer cell paper
state.gene <- read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/TableS2.CancaerCell.txt",head=T,sep="\t")
state.name <- colnames(state.gene)
meta <- Kras@meta.data

for(i in 2:ncol(state.gene)){
  df = state.gene[order(state.gene[,i],decreasing=T),]
  df.genes = as.character(df[1:10,1])
  df.score = FetchData(Kras,df.genes)
  meta[,state.name[i]] = rowMeans(df.score)
}

meta.cut = meta[,c(paste0("NMF_0",1:9),"NMF_10","NMF_11")]
Kras <- AddMetaData(Kras,meta.cut )  # add column data

# calculate the median of expression across cells in clusters
state.exp.mid = aggregate(meta[,c(paste0("NMF_0",1:9),"NMF_10","NMF_11")], list(meta$seurat_clusters), FUN=median)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
# name of cancer scores
cvt.name = read.table('/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/TableS2.CancaerCell.names.txt',head=F,sep='\t')
colnames(state.exp.mid) = cvt.name$V2
saveRDS(state.exp.mid,'Aggregated.cancerState.gene.scores.tumor_cls.state.exp.mid.20220518.top10.rds')

# calculate the median of expression across cells in populations
state.exp.mid = aggregate(meta[,c(paste0("NMF_0",1:9),"NMF_10","NMF_11")], list(meta$Population1), FUN=median)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
# name of cancer scores
cvt.name = read.table('/n/data2/bch/hemonc/ckim/JINGYUN/Mouse.KRAS.scRNA.20201208/Refs/TableS2.CancaerCell.names.txt',head=F,sep='\t')
colnames(state.exp.mid) = cvt.name$V2
saveRDS(state.exp.mid,'Aggregated.cancerState.gene.scores.tumor_cls.state.exp.mid_cross_population1.20220518.top10.rds')

########### Marker genes overlapped TFs and surface markers

####### find marker genes for each cluster
Idents(Kras) <- "Population1"
cls.list = unique(Kras@meta.data$Population1)
for(cls in cls.list){
 cls.marker <- FindMarkers(Kras,ident.1=cls,only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
 write.table(data.frame(gene=rownames(cls.marker),cls.marker),paste0("MarkerGene.cluster",cls,".KPY_only.0518.tumor_cls.txt"),quote=F,row.names=F,sep="\t")
}

## overlap and show matrix
tf = read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Reference/mouse/TF.list.txt",head=T)
tf.genes = unique(as.character(tf$TFs))

cls.list = unique(Kras@meta.data$Population1)
markers.merge = data.frame()
for(cls in cls.list){
 cls.marker <- read.table(paste0("MarkerGene.cluster",cls,".KPY_only.0518.tumor_cls.txt"),head=T)
 cls.marker = cls.marker[cls.marker$p_val_adj<0.05 & cls.marker$avg_log2FC>log2(1.5),]
 cls.marker = cls.marker[cls.marker$gene %in% tf.genes,]
 if(nrow(cls.marker)==0){next}
# cls.marker = cls.marker[1:5,]
 cls.marker$cls = paste0("c",cls)
 markers.merge = rbind(markers.merge,cls.marker)
}
markers.genes = unique(as.character(markers.merge$gene))
tf.markers = c(markers.genes,tf.genes)
tf.markers =  tf.markers[duplicated(c(markers.genes,tf.genes))]

tf.markers.exp = FetchData(Kras,tf.markers)
meta = Kras@meta.data
all(rownames(meta) == rownames(tf.markers.exp))  ## check cells
meta = cbind(meta,tf.markers.exp)
## expression pseudobulk
# mean for each cluster
marker.df = aggregate(meta[,tf.markers], list(meta$seurat_clusters), FUN=mean)
rownames(marker.df) = marker.df[,1]
marker.df = marker.df[,-1]
saveRDS(marker.df,"TF.Marker.expression.mean_cross_cluster.20220518.rds")

# mean for each population
marker.df = aggregate(meta[,tf.markers], list(meta$Population1), FUN=mean)
rownames(marker.df) = marker.df[,1]
marker.df = marker.df[,-1]
saveRDS(marker.df,"TF.Marker.expression.mean_cross_population1.20220518.rds")

## overlap with surface markers
surface = read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Kras.P53.scRNA.20210805/2.integration.with.KRAS/Surface.Marker.txt",head=F)
surface.genes = as.character(surface$V1)
surface.genes = surface.genes[order(surface.genes)]
surface.genes = unique(surface.genes)

cls.list = unique(Kras@meta.data$Population1)
markers.merge = data.frame()
for(cls in cls.list){
 cls.marker <- read.table(paste0("MarkerGene.cluster",cls,".KPY_only.0518.tumor_cls.txt"),head=T)
 cls.marker = cls.marker[cls.marker$p_val_adj<0.05 & cls.marker$avg_log2FC>log2(1.5),]
 cls.marker = cls.marker[cls.marker$gene %in% surface.genes,]
 if(nrow(cls.marker)==0){next}
# cls.marker = cls.marker[1:6,]
 cls.marker$cls = paste0(cls)
 markers.merge = rbind(markers.merge,cls.marker)
}
markers.genes = unique(as.character(markers.merge$gene))
surface.markers = c(markers.genes,surface.genes)
surface.markers =  surface.markers[duplicated(c(markers.genes,surface.genes))]

# add two genes Procr F3: based on Wechat Jingyun Sent on Nov 17, 2021
surface.markers = c(surface.markers,"Procr","F3")

surface.markers.exp = FetchData(Kras,surface.markers)
meta = Kras@meta.data
all(rownames(meta) == rownames(surface.markers.exp))  ## check cells
meta = cbind(meta,surface.markers.exp)
## expression pseudobulk
# mean for each cluster
marker.df = aggregate(meta[,surface.markers], list(meta$seurat_clusters), FUN=mean)
rownames(marker.df) = marker.df[,1]
marker.df = marker.df[,-1]
saveRDS(marker.df,"Surface.Marker.expression.mean_cross_cluster.20220518.rds")

# mean for each population
marker.df = aggregate(meta[,surface.markers], list(meta$Population1), FUN=mean)
rownames(marker.df) = marker.df[,1]
marker.df = marker.df[,-1]
saveRDS(marker.df,"Surface.Marker.expression.mean_cross_population1.20220518.rds")

######## find marker for cluster B
Idents(Kras) <- "Population1"

Hmga2_B.vs.A <- FindMarkers(Kras,ident.1='Hmga2_B',ident.2 = 'Hmga2_A',only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
Hmga2_B.vs.Nkx2 <- FindMarkers(Kras,ident.1='Hmga2_B',ident.2 = c('Nkx2-1_A','Nkx2-1_B'),only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
Hmga2_B.vs.A$gene = rownames(Hmga2_B.vs.A)
Hmga2_B.vs.Nkx2$gene = rownames(Hmga2_B.vs.Nkx2)
Hmga2_B.marker = merge(Hmga2_B.vs.A,Hmga2_B.vs.Nkx2,suffix=c('Hmga2_B.vs.A','Hmga2_B.vs.Nkx2'),by='gene')

Nkx2_B.vs.A = FindMarkers(Kras,ident.1='Nkx2-1_B',ident.2 = 'Nkx2-1_A',only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
Nkx2_B.vs.Hmga2 <- FindMarkers(Kras,ident.1='Nkx2-1_B',ident.2 = c('Hmga2_A','Hmga2_B'),only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
Nkx2_B.vs.A$gene = rownames(Nkx2_B.vs.A)
Nkx2_B.vs.Hmga2$gene = rownames(Nkx2_B.vs.Hmga2)
Nkx2_B.marker = merge(Nkx2_B.vs.A,Nkx2_B.vs.Hmga2,suffix=c('Nkx2_B.vs.A','Nkx2_B.vs.Hmga2'),by='gene')

surface = read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Kras.P53.scRNA.20210805/2.integration.with.KRAS/Surface.Marker.txt",head=F)
surface.genes = as.character(surface$V1)
surface.genes = surface.genes[order(surface.genes)]
surface.genes = unique(surface.genes)

Hmga2_B.marker.surf = Hmga2_B.marker[Hmga2_B.marker$gene %in% surface.genes,]
Nkx2_B.marker.surf = Nkx2_B.marker[Nkx2_B.marker$gene %in% surface.genes,]

# plot 
# Hmga2
pdf('KPY_only.Hmga2_B.marker.surface.pdf')
FeaturePlot(Kras, features = Hmga2_B.marker.surf$gene,cols = c("lightgrey","firebrick3"))
dev.off()

pdf('KPY_only.Nkx2_B.marker.surface.pdf',width=7,height=10)
FeaturePlot(Kras, features = Nkx2_B.marker.surf$gene,cols = c("lightgrey","firebrick3"))
dev.off()




