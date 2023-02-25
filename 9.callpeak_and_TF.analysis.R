# for ArchR, use R3.6.1
library(ArchR)
set.seed(1)
addArchRThreads(threads = 6)
addArchRGenome("mm10")   ## use mm10
pathToMacs2 <- findMacs2()

# read in ATAC object
proj = loadArchRProject(path = "Save-Proj2.RNA_subseted")

## peak calling
proj <- addGroupCoverages(ArchRProj = proj,groupBy = "Population1")
pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    peakMethod = "Macs2",
    groupBy = "Population1",
    force = TRUE,
    maxPeaks=300000,
    cutOff = 0.01,
    pathToMacs2 = pathToMacs2
)
getPeakSet(proj)

peaks = getPeakSet(proj)
saveRDS(peaks,'Peaks/Peaks.ArchR.by_population1.rds')
peaks = data.frame(peaks)
write.table(peaks,'Peaks/Peaks.ArchR.by_population1.txt',row.names=F,quote=F,sep='\t')

# add peak matrix
proj <- addPeakMatrix(proj)

# DEP in 2 groups
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Population2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
# get a glimps
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

# do TF enrichment test
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

Hmga2 = data.frame(TF = rownames(enrichMotifs), mlog10Padj = assays(enrichMotifs)[['mlog10Padj']][,'Hmga2'],mlog10p=assays(enrichMotifs)[['mlog10p']][,'Hmga2'],Enrichment=assays(enrichMotifs)[['Enrichment']][,'Hmga2'])
Nkx2 = data.frame(TF = rownames(enrichMotifs), mlog10Padj = assays(enrichMotifs)[['mlog10Padj']][,'Nkx2-1'],mlog10p=assays(enrichMotifs)[['mlog10p']][,'Nkx2-1'],Enrichment=assays(enrichMotifs)[['Enrichment']][,'Nkx2-1'])

tf.merge = merge(Hmga2,Nkx2,by='TF',suffix=c('.Hmga2','.Nkx2'))
library(stringr)
tf.merge$TF.name <- str_split_fixed(tf.merge$TF,"_",2)[,1]
tf.merge$group = 'nc'
tf.merge$group[tf.merge$mlog10p.Hmga2>2 & tf.merge$Enrichment.Hmga2>1] = 'Hmga2.enriched'
tf.merge$group[tf.merge$mlog10p.Nkx2>2 & tf.merge$Enrichment.Nkx2>1] = 'Nkx2.enriched'
tf.merge$group[tf.merge$mlog10p.Hmga2>2 & tf.merge$Enrichment.Hmga2>1 & tf.merge$mlog10p.Nkx2>2 & tf.merge$Enrichment.Nkx2>1] = 'Both.enriched'

tf.merge.Hmga = head(tf.merge[order(tf.merge$Enrichment.Hmga,decreasing=T),],50)
tf.merge.Nkx2 = head(tf.merge[order(tf.merge$Enrichment.Nkx2,decreasing=T),],50)
merge.top = rbind(tf.merge.Hmga,tf.merge.Nkx2)
merge.top = merge.top[! merge.top$group %in% 'nc',]

# ggplot
library(ggplot2)
library(ggpubr)
library(ggrepel)

pdf('TF.enrichment.2group.pdf',useDingbats=F)
ggplot(tf.merge, aes(log2(Enrichment.Hmga2), log2(Enrichment.Nkx2),color=group)) +
        geom_point(cex=.5) +
        geom_text_repel(data=merge.top, aes(log2(Enrichment.Hmga2), log2(Enrichment.Nkx2), label=TF.name,color=group), size=2, max.overlaps=50) +
        scale_y_continuous(expand=c(0,0)) + xlim(-2,2) + ylim(-2,2)+
        theme_pubr() + theme(legend.position='none') + ggtitle('TF enrichments in diff. peaks') +
        geom_hline(yintercept = 0,color="grey1",linetype="dashed") + geom_vline(xintercept = 0,color="grey1",linetype="dashed") 

write.table(tf.merge,'TF.enrichment.2group.txt',row.names=F,quote=F,sep='\t')
dev.off()

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj3.RNA_subseted", load = FALSE)

# look at the gene expression
# R 4.1.1 required (reopen R)
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

Idents(Kras) <- "Population2"
method = 'wilcox'

Kras.DEGs <- FindAllMarkers(Kras, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0,test.use=method)

Kras.DEGs.Hmga2 = Kras.DEGs[Kras.DEGs$cluster == 'Hmga2',]
Kras.DEGs.Nkx2 = Kras.DEGs[Kras.DEGs$cluster == 'Nkx2-1',]

# read in TF enrichment results
tf.merge = read.table('TF.enrichment.2group.txt',head=T)
colnames(tf.merge) = paste0(colnames(tf.merge),'.TF')
tf.merge$TF.name.TF[tf.merge$TF.name.TF=='Nkx21'] = 'Nkx2-1'

## work on Hmga2
# first flip DEG to Hmga2
Kras.DEGs.Nkx2_flip = Kras.DEGs.Nkx2
Kras.DEGs.Nkx2_flip$avg_log2FC = -Kras.DEGs.Nkx2$avg_log2FC
Kras.DEGs.Nkx2_flip$pct.1 = Kras.DEGs.Nkx2$pct.2
Kras.DEGs.Nkx2_flip$pct.2 = Kras.DEGs.Nkx2$pct.1
Kras.DEGs.re = rbind(Kras.DEGs.Hmga2,Kras.DEGs.Nkx2_flip)

colnames(Kras.DEGs.re) = paste0(colnames(Kras.DEGs.re),'.DEG')
merge = merge(tf.merge,Kras.DEGs.re,by.x='TF.name.TF',by.y='gene.DEG')
merge.top = rbind(head(merge[order(merge$avg_log2FC.DEG,decreasing=T),],40),head(merge[order(merge$Enrichment.Hmga2,decreasing=T),],40))
merge.top = unique(merge.top)

library(ggplot2)
#library(ggpubr)
library(ggrepel)

pdf('TF.enrichment.vs.DEG.2group.pdf',useDingbats=F)

ggplot(merge, aes(avg_log2FC.DEG, log2(Enrichment.Hmga2.TF),color=pct.1.DEG)) +
        geom_point(cex=.5) +
        geom_text_repel(data=merge.top, aes(avg_log2FC.DEG, log2(Enrichment.Hmga2.TF), label=TF.name.TF), size=2, max.overlaps=50) +
        scale_y_continuous(expand=c(0,0)) + xlim(-2,2) + ylim(-2,2) + #+ theme(legend.position='none') +
        theme_classic() + ggtitle('TF enrich vs. DEG; for Hmga2 cluster') +
        geom_hline(yintercept = 0,color="grey1",linetype="dashed") + geom_vline(xintercept = 0,color="grey1",linetype="dashed")

# use log(Hmga2/NKx2) TF as y axis
merge$log2FC.TF_enrich = log2(merge$Enrichment.Hmga2.TF/merge$Enrichment.Nkx2.TF)
merge.top$log2FC.TF_enrich = log2(merge.top$Enrichment.Hmga2.TF/merge.top$Enrichment.Nkx2.TF)

ggplot(merge, aes(avg_log2FC.DEG, log2FC.TF_enrich,color=pct.1.DEG)) +
        geom_point(cex=.5) +
        geom_text_repel(data=merge.top, aes(avg_log2FC.DEG, log2FC.TF_enrich, label=TF.name.TF), size=2, max.overlaps=50) +
        scale_y_continuous(expand=c(0,0)) + xlim(-2,2) + ylim(-2,2)+ #theme(legend.position='none') +
        theme_classic() + ggtitle('TF enrich Hmga2/Nkx2 vs. DEG; for Hmga2 cluster') +
        geom_hline(yintercept = 0,color="grey1",linetype="dashed") + geom_vline(xintercept = 0,color="grey1",linetype="dashed")

## work on Nkx2
# first flip DEG to Nkx2
Kras.DEGs.Hmga2_flip = Kras.DEGs.Hmga2
Kras.DEGs.Hmga2_flip$avg_log2FC = -Kras.DEGs.Hmga2$avg_log2FC
Kras.DEGs.Hmga2_flip$pct.1 = Kras.DEGs.Hmga2$pct.2
Kras.DEGs.Hmga2_flip$pct.2 = Kras.DEGs.Hmga2$pct.1
Kras.DEGs.re = rbind(Kras.DEGs.Nkx2,Kras.DEGs.Hmga2_flip)

colnames(Kras.DEGs.re) = paste0(colnames(Kras.DEGs.re),'.DEG')
merge = merge(tf.merge,Kras.DEGs.re,by.x='TF.name.TF',by.y='gene.DEG')
merge.top = rbind(head(merge[order(merge$avg_log2FC.DEG,decreasing=T),],40),head(merge[order(merge$Enrichment.Nkx2,decreasing=T),],40))
merge.top = unique(merge.top)

ggplot(merge, aes(avg_log2FC.DEG, log2(Enrichment.Nkx2.TF),color=pct.1.DEG)) +
        geom_point(cex=.5) +
        geom_text_repel(data=merge.top, aes(avg_log2FC.DEG, log2(Enrichment.Nkx2.TF), label=TF.name.TF), size=2, max.overlaps=50) +
        scale_y_continuous(expand=c(0,0)) + xlim(-2,2) + ylim(-2,2) + #+ theme(legend.position='none') +
        theme_classic() + ggtitle('TF enrich vs. DEG; for Nkx2 cluster') +
        geom_hline(yintercept = 0,color="grey1",linetype="dashed") + geom_vline(xintercept = 0,color="grey1",linetype="dashed")

# use log(Hmga2/NKx2) TF as y axis
merge$log2FC.TF_enrich = log2(merge$Enrichment.Nkx2.TF/merge$Enrichment.Hmga2.TF)
merge.top$log2FC.TF_enrich = log2(merge.top$Enrichment.Nkx2.TF/merge.top$Enrichment.Hmga2.TF)
ggplot(merge, aes(avg_log2FC.DEG, log2FC.TF_enrich,color=pct.1.DEG)) +
        geom_point(cex=.5) +
        geom_text_repel(data=merge.top, aes(avg_log2FC.DEG, log2FC.TF_enrich, label=TF.name.TF), size=2, max.overlaps=50) +
        scale_y_continuous(expand=c(0,0)) + xlim(-2,2) + ylim(-2,2)+ #theme(legend.position='none') +
        theme_classic() + ggtitle('TF enrich Nkx2/Hmga2 vs. DEG; for Nkx2 cluster') +
        geom_hline(yintercept = 0,color="grey1",linetype="dashed") + geom_vline(xintercept = 0,color="grey1",linetype="dashed")

dev.off()




