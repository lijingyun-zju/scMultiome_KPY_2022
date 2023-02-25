# R 4.1.1 required
library(Signac)
library(Seurat)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(stringr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(SCENIC)
library(SummarizedExperiment)

set.seed(1)

## 1_7.scenic.longtime.R was run first; 1_7.SCENIC_Kras.re was copied to 1_7.SCENIC_Kras

setwd('1_7.SCENIC_Kras')

# read in the output of 1_7.scenic.longtime.R 
load('./run.scenic.tmp2.RData')

# re-Initialize SCENIC settings
org <- "mgi" # mouse
dbDir <- "/n/data2/bch/hemonc/ckim/JINGYUN/software/cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC on Mouse Kras" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs['500bp'] = 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'
dbs['10kb'] = 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather'

scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 1

# real run 
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions ) #
# Optional: log expression (for TF expression plot, it does not affect any other calculation)
exprMat_log <- log2(exprMat_filtered+1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status

###### for plot
# plot regulon activity (AUC) on UMAP! 
# first read in AUC regulon
setwd('1_7.SCENIC_Kras')
auc = readRDS('/n/data2/bch/hemonc/ckim/JINGYUN/MultiOmic.20220429/1_7.SCENIC_Kras/int/3.4_regulonAUC.Rds')
auc.df = data.frame(t(data.frame(assays(auc)$AUC)))
saveRDS(auc.df,'/n/data2/bch/hemonc/ckim/JINGYUN/MultiOmic.20220429/1_7.SCENIC_Kras/int/3.4_regulonAUC.df.rds')

# read in seurat for UMAP and annotation
Kras = readRDS("../Kras_P53.2Kgene.RNA.KPY_only.Output.nfeat2000.nPC75.rm_cellcycle.assigned.rds")
meta = Kras@meta.data
umapCoord <- as.data.frame(Embeddings(object = Kras[["umap"]]))
all(rownames(umapCoord) == rownames(meta))
meta$UMAP1.rna = umapCoord$UMAP_1
meta$UMAP2.rna = umapCoord$UMAP_2

rownames(meta) = gsub(':','.',rownames(meta))
rownames(meta) = gsub('-','.',rownames(meta))

# check name order
all(rownames(meta) == rownames(auc.df))
# assign umap coordinate and population
all.tfs = colnames(auc.df)
auc.df[,c('Population1','Population2','UMAP1.rna','UMAP2.rna')] = meta[,c('Population1','Population2','UMAP1.rna','UMAP2.rna')]

# plot Nkx2-1
auc.df$Population1 = factor(auc.df$Population1,levels= rev(c('Nkx2-1_A','Nkx2-1_B','Hmga2_B','Hmga2_A')))
pop1.cols = rev(c('aquamarine4','aquamarine1','pink','pink4'))

library(ggplot2)
mytheme2 <-theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),aspect.ratio = 1,axis.text = element_text(face="bold", size=rel(1)),axis.title = element_text(face="bold", size=rel(1)),legend.text=element_text(face="bold",size=10),plot.title=element_text(face="bold",size=rel(1)),legend.title =element_text(face="bold", size=rel(1)))
p.list = list()
cols = c("grey","yellow", "darkred")
png(paste0("Plot/Target_gene.exp.",'Nkx2.1',".umap.png"),width=7, height=7,units='in', res=200)
ggplot(auc.df,aes(x=UMAP1.rna,y=UMAP2.rna,color=Nkx2.1..17g.)) + geom_point(shape=16,size=1) + mytheme2 + ggtitle(paste0('Nkx2-1',' target gene.exp')) + scale_colour_gradientn(colours = cols)
dev.off()

pdf(paste0("Plot/Target_gene.exp.",'Nkx2.1',".pdf"),width=7, height=7,useDingbats=F)
ggplot(auc.df,aes(x=Population1,y=Nkx2.1..17g.,color=Population1)) + geom_boxplot() + mytheme2 + theme(aspect.ratio=2) + ggtitle(paste0('Nkx2-1',' target gene.exp boxplot')) + ylab(paste0('Nkx2-1',' target gene.exp')) + scale_colour_manual(values=pop1.cols)
dev.off()

## plot other TFs
AT2.TFs = c('Dbp','Foxa2','Cebpa','Rfx2')
# boxplot
pdf(paste0("Plot/Target_gene.exp.AT2_TFs.pdf"),width=9, height=7,useDingbats=F)
for(TF in AT2.TFs){
	TF2 = all.tfs[grepl(TF,all.tfs)]
	for(k in 1:length(TF2)){
	 p = ggplot(auc.df,aes(x=Population1,y=auc.df[,TF2[[k]]],color=Population1)) + geom_boxplot() + mytheme2 + theme(aspect.ratio=2) + ggtitle(paste0(TF,' target gene.exp boxplot')) + ylab(paste0(TF2[[k]],' target gene.exp')) + scale_colour_manual(values=pop1.cols)
	 print(p)
	}
}
dev.off()
# umaps
for(TF in AT2.TFs){
        TF2 = all.tfs[grepl(TF,all.tfs)]
        print(TF2)
        for(k in 1:length(TF2)){
	 png(paste0("Plot/Target_gene.exp.",TF2[[k]],".umap.png"),width=7, height=7,units='in', res=200)
	 p = ggplot(auc.df,aes(x=UMAP1.rna,y=UMAP2.rna,color=auc.df[,TF2[[k]]])) + geom_point(shape=16,size=1) + mytheme2 + ggtitle(paste0(TF2[[k]],' target gene.exp')) + scale_colour_gradientn(colours = cols)
         print(p)
	 dev.off()
        }
}

EMT.TFs = c('Ets2','Tcf12','Myc','Hivep2','Nfkb1','Rel','Tcf4','Sox4')
pdf(paste0("Plot/Target_gene.exp.EMT_TFs.pdf"),width=7, height=7,useDingbats=F)
for(TF in EMT.TFs){
        TF2 = all.tfs[grepl(TF,all.tfs)]
        for(k in 1:length(TF2)){
         p = ggplot(auc.df,aes(x=Population1,y=auc.df[,TF2[[k]]],color=Population1)) + geom_boxplot() + mytheme2 + theme(aspect.ratio=2) + ggtitle(paste0(TF,' target gene.exp boxplot')) + ylab(paste0(TF2[[k]],' target gene.exp')) + scale_colour_manual(values=pop1.cols)
         print(p)
        }
}
dev.off()
# umaps
for(TF in EMT.TFs){
        TF2 = all.tfs[grepl(TF,all.tfs)]
	print(TF2)
        for(k in 1:length(TF2)){
         png(paste0("Plot/Target_gene.exp.",TF2[[k]],".umap.png"),width=7, height=7,units='in', res=200)
         p = ggplot(auc.df,aes(x=UMAP1.rna,y=UMAP2.rna,color=auc.df[,TF2[[k]]])) + geom_point(shape=16,size=1) + mytheme2 + ggtitle(paste0(TF2[[k]],' target gene.exp')) + scale_colour_gradientn(colours = cols)
         print(p)
         dev.off()
        }
}

## 2022/10/01; plot Zeb1
pdf(paste0("Plot/Target_gene.exp.Zeb1.pdf"),width=7, height=7,useDingbats=F)
TF = 'Zeb1'
        TF2 = all.tfs[grepl(TF,all.tfs)]
        for(k in 1:length(TF2)){
         p = ggplot(auc.df,aes(x=Population1,y=auc.df[,TF2[[k]]],color=Population1)) + geom_boxplot() + mytheme2 + theme(aspect.ratio=2) + ggtitle(paste0(TF,' target gene.exp boxplot')) + ylab(paste0(TF2[[k]],' target gene.exp')) + scale_colour_manual(values=pop1.cols)
         print(p)
        }
dev.off()
# umaps
        TF2 = all.tfs[grepl(TF,all.tfs)]
        print(TF2)
        for(k in 1:length(TF2)){
         png(paste0("Plot/Target_gene.exp.",TF2[[k]],".umap.png"),width=7, height=7,units='in', res=200)
         p = ggplot(auc.df,aes(x=UMAP1.rna,y=UMAP2.rna,color=auc.df[,TF2[[k]]])) + geom_point(shape=16,size=1) + mytheme2 + ggtitle(paste0(TF2[[k]],' target gene.exp')) + scale_colour_gradientn(colours = cols)
         print(p)
         dev.off()
        }

### boxplot
# calculate the median of target-gene expression across cells in clusters
TFs.merged = c('Nkx2.1..17g.','Dbp_extended..26g.','Foxa2_extended..86g.','Cebpa_extended..287g.','Rfx2_extended..50g.','Ets2_extended..8265g.','Tcf12_extended..15g.','Myc_extended..20g.','Hivep2_extended..74g.','Nfkb1_extended..50g.','Rel_extended..94g.','Tcf4_extended..85g.','Sox4_extended..471g.','Zeb1_extended..27g.')
state.exp.mid = aggregate(auc.df[,TFs.merged], list(auc.df$Population1), FUN=median)
library(pheatmap)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
colnames(state.exp.mid) = c('Nkx2.1','Dbp','Foxa2','Cebpa','Rfx2','Ets2','Tcf12','Myc','Hivep2','Nfkb1','Rel','Tcf4','Sox4')
state.exp.mid.t = t(state.exp.mid)
state.exp.mid.t = state.exp.mid.t[,c('Nkx2-1_A','Nkx2-1_B','Hmga2_B','Hmga2_A')]
pdf('Target_gene.exp.heatmap.pdf',useDingbats=F,width=5)
#pheatmap(t(state.exp.mid),scale = "none",cluster_row = F,cluster_cols = T,main="Target_gene exp no_norm",border_color=NA)
pheatmap(state.exp.mid.t,scale = "row",cluster_row = F,cluster_cols = T,main="Target_gene exp by row",border_color=NA)
pheatmap(state.exp.mid.t,scale = "row",cluster_row = F,cluster_cols = F,main="Target_gene exp by row",border_color=NA)
dev.off()

# calculate the median of the TF expression 
TFs.merged = c('Nkx2-1','Dbp','Foxa2','Cebpa','Rfx2','Ets2','Tcf12','Myc','Hivep2','Nfkb1','Rel','Tcf4','Sox4')
meta <- Kras@meta.data
for(tf in TFs.merged){
	meta[,tf] = FetchData(Kras,tf)
}

state.exp.mid = aggregate(meta[,TFs.merged], list(meta$Population1), FUN=mean)
library(pheatmap)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
state.exp.mid.t = t(state.exp.mid)
state.exp.mid.t = state.exp.mid.t[,c('Nkx2-1_A','Nkx2-1_B','Hmga2_B','Hmga2_A')]
pdf('TF.exp.heatmap.pdf',useDingbats=F,width=5)
#pheatmap(t(state.exp.mid),scale = "none",cluster_row = F,cluster_cols = T,main="Target_gene exp no_norm",border_color=NA)
pheatmap(state.exp.mid.t,scale = "row",cluster_row = F,cluster_cols = T,main="TF exp by row",border_color=NA)
pheatmap(state.exp.mid.t,scale = "row",cluster_row = F,cluster_cols = F,main="TF exp by row",border_color=NA)
dev.off()

# re-draw heatmap using the intersected TFs
TFs = readRDS('1_10.TF.overlaps/TFs.merged.new.rds')
auc.df = readRDS('/n/data2/bch/hemonc/ckim/JINGYUN/MultiOmic.20220429/1_7.SCENIC_Kras/int/3.4_regulonAUC.df.rds')

Kras = readRDS("./Kras_P53.2Kgene.RNA.KPY_only.Output.nfeat2000.nPC75.rm_cellcycle.assigned.rds")
meta = Kras@meta.data
umapCoord <- as.data.frame(Embeddings(object = Kras[["umap"]]))
all(rownames(umapCoord) == rownames(meta))
meta$UMAP1.rna = umapCoord$UMAP_1
meta$UMAP2.rna = umapCoord$UMAP_2

rownames(meta) = gsub(':','.',rownames(meta))
rownames(meta) = gsub('-','.',rownames(meta))

# check name order
all(rownames(meta) == rownames(auc.df))
# assign umap coordinate and population
all.tfs = colnames(auc.df)
auc.df[,c('Population1','Population2','UMAP1.rna','UMAP2.rna')] = meta[,c('Population1','Population2','UMAP1.rna','UMAP2.rna')]

TFs.merged = c("Cebpa_extended..287g.","Foxa2_extended..86g.","Nkx2.1..17g.","Stat3_extended..70g.","Tcf12_extended..15g.","Tcf4_extended..85g.","Rel..70g.","Nfkb1_extended..50g.","Sox4_extended..471g.","Hivep2_extended..74g.")

state.exp.mid = aggregate(auc.df[,TFs.merged], list(auc.df$Population1), FUN=median)
library(pheatmap)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
colnames(state.exp.mid) = TFs

state.exp.mid.t = t(state.exp.mid)
state.exp.mid.t = state.exp.mid.t[,c('Nkx2-1_A','Nkx2-1_B','Hmga2_B','Hmga2_A')]
pdf('Target_gene.exp.heatmap.0716_re.pdf',useDingbats=F,width=5)
#pheatmap(t(state.exp.mid),scale = "none",cluster_row = F,cluster_cols = T,main="Target_gene exp no_norm",border_color=NA)
pheatmap(state.exp.mid.t,scale = "row",cluster_row = F,cluster_cols = T,main="Target_gene exp by row",border_color=NA)
pheatmap(state.exp.mid.t,scale = "row",cluster_row = F,cluster_cols = F,main="Target_gene exp by row",border_color=NA)
dev.off()

# calculate the median of the TF expression
TFs.merged = TFs
meta <- Kras@meta.data
for(tf in TFs.merged){
        meta[,tf] = FetchData(Kras,tf)
}

state.exp.mid = aggregate(meta[,TFs.merged], list(meta$Population1), FUN=mean)
library(pheatmap)
rownames(state.exp.mid) = state.exp.mid[,1]
state.exp.mid = state.exp.mid[,-1]
state.exp.mid.t = t(state.exp.mid)
state.exp.mid.t = state.exp.mid.t[,c('Nkx2-1_A','Nkx2-1_B','Hmga2_B','Hmga2_A')]
pdf('TF.exp.heatmap.0716_re.pdf',useDingbats=F,width=5)
#pheatmap(t(state.exp.mid),scale = "none",cluster_row = F,cluster_cols = T,main="Target_gene exp no_norm",border_color=NA)
pheatmap(state.exp.mid.t,scale = "row",cluster_row = F,cluster_cols = T,main="TF exp by row",border_color=NA)
pheatmap(state.exp.mid.t,scale = "row",cluster_row = F,cluster_cols = F,main="TF exp by row",border_color=NA)
dev.off()






