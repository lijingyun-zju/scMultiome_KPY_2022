# for ArchR, use R3.6.1
library(ArchR)
set.seed(1)
addArchRThreads(threads = 6)
addArchRGenome("mm10")   ## use mm10
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
register(MulticoreParam(6, progressbar = TRUE))


# read in ATAC object
proj = loadArchRProject(path = "Save-Proj2.RNA_subseted")

if(0){  # has already done
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
}

# get peak matrix
mat = getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

library(BSgenome.Mmusculus.UCSC.mm10)
mat <- addGCBias(mat, genome = BSgenome.Mmusculus.UCSC.mm10)

# read in cisbp motif from chromVARmotifs! 
# https://github.com/GreenleafLab/chromVARmotifs
load('/n/data2/bch/hemonc/ckim/JINGYUN/software/chromVARmotifs/data/mouse_pwms_v2.rda')
motifs = mouse_pwms_v2

library(motifmatchr)
names(assays(mat)) = 'counts' # rename 
motif_ix <- matchMotifs(motifs, mat, 
                        genome = BSgenome.Mmusculus.UCSC.mm10)

dev <- computeDeviations(object = mat, annotations = motif_ix)

saveRDS(dev,'ChromVar.Kras.rds')
# use background
bg <- getBackgroundPeaks(object = mat)
dev2 <- computeDeviations(object = mat, annotations = motif_ix,
                         background_peaks = bg)
saveRDS(dev2,'ChromVar.Kras.useBgd.rds')
save.image('ChromVar.tmp.Rdata')


# plot on UMAP; 
# genes to plot
tf.list = c('sox4','hivep2','ets2','tcf12','myc','rel','tcf4','nfkb1','dbp','foxa2','rfx2','cebpa','nkx21')  # nkx21 is actually nkx2-1
library(Hmisc)
tf.list = capitalize(tf.list)
mytheme2 <-theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),aspect.ratio = 1,axis.text = element_text(face="bold", size=rel(1)),axis.title = element_text(face="bold", size=rel(1)),legend.text=element_text(face="bold",size=10),plot.title=element_text(face="bold",size=rel(1)),legend.title =element_text(face="bold", size=rel(1)))
# first get RNA	umap coordinates
rna.meta = readRDS('Kras_P53.2Kgene.RNA.KPY_only.Output.nfeat2000.nPC75.rm_cellcycle.meta.rds')
library(stringr)
rna.meta$cell = str_split_fixed(rownames(rna.meta),":",2)[,2]
rownames(rna.meta) = rna.meta$cell

dev.meta = data.frame(colData(dev))
dev.meta$cell = str_split_fixed(rownames(dev.meta),"#",2)[,2]

rna.meta = rna.meta[rna.meta$cell %in% dev.meta$cell,]
rna.meta = rna.meta[dev.meta$cell,]
all(rna.meta$cell == dev.meta$cell)
dev.meta$UMAP1.rna = rna.meta$UMAP1.rna
dev.meta$UMAP2.rna = rna.meta$UMAP2.rna

z.mat = t(assays(dev)$z)

all(rownames(dev.meta) == rownames(z.mat))  # check cell order
# do the work!

p.box.list = list()
k = 1
all.tf = rownames(dev)
dev.meta$Population1 = factor(dev.meta$Population1,levels= rev(c('Nkx2-1_A','Nkx2-1_B','Hmga2_B','Hmga2_A')))
pop1.cols = rev(c('aquamarine4','aquamarine1','pink','pink4'))

for(tf.name in tf.list){
tf.tmp = paste0("_",tf.name,"_")  # to avoid multiple matching
tf = all.tf[grepl(tf.tmp,all.tf)]
print(tf)
dev.meta$tf = z.mat[,tf]
p.list = list()
#cols = c("gray80","steelblue4")
cols = c("cyan3","yellow", "darkred")
p.list[[1]] = ggplot(dev.meta,aes(x= -UMAP1,y=- UMAP2,color=tf)) + geom_point(shape=16,size=1) + mytheme2 + ggtitle(paste0(tf.name,' motif-z ATAC-umap')) + scale_colour_gradientn(colours = cols)
p.list[[2]] = ggplot(dev.meta,aes(x= UMAP1.rna,y= UMAP2.rna,color=tf)) + geom_point(shape=16,size=1) + mytheme2 + ggtitle(paste0(tf.name,' motif-z RNA-umap')) + scale_colour_gradientn(colours = cols)
p.list[[3]] = ggplot(dev.meta,aes(x=Population1,y=tf,color=Population1)) + geom_boxplot() + mytheme2 + theme(aspect.ratio=2) + ggtitle(paste0(tf.name,' motif-z boxplot')) + ylab('TF chromVar') + scale_colour_manual(values=pop1.cols)

library(gridExtra)
png(paste0("Motif.ChromVar.out/ChromVar.no_bgd.",tf,"z.png"),width=16, height=7,units='in', res=200)
do.call("grid.arrange", c(p.list, nrow=1))
dev.off()
}
#pdf('Motif.ChromVar.out/ChromVar.no_bgd.boxplot.pdf',useDingbats=F)
#print(p.box.list)
#dev.off()


# use the one with bgd
# first get RNA	umap coordinates
rna.meta = readRDS('Kras_P53.2Kgene.RNA.KPY_only.Output.nfeat2000.nPC75.rm_cellcycle.meta.rds')
library(stringr)
rna.meta$cell = str_split_fixed(rownames(rna.meta),":",2)[,2]
rownames(rna.meta) = rna.meta$cell

dev2.meta = data.frame(colData(dev2))
dev2.meta$cell = str_split_fixed(rownames(dev2.meta),"#",2)[,2]

rna.meta = rna.meta[rna.meta$cell %in% dev2.meta$cell,]
rna.meta = rna.meta[dev2.meta$cell,]
all(rna.meta$cell == dev2.meta$cell)
dev2.meta$UMAP1.rna = rna.meta$UMAP1.rna
dev2.meta$UMAP2.rna = rna.meta$UMAP2.rna

z.mat = t(assays(dev2)$z)

all(rownames(dev2.meta) == rownames(z.mat))  # check cell order
# do the work!

#tf = 'ENSMUSG00000056758_LINE31_Hmga2_D'
tf = 'ENSMUSG00000001496_LINE1147_Nkx21_D'
dev2.meta$tf = z.mat[,tf]
p.list = list()
#cols = c("gray80","steelblue4")
cols = c("cyan3","yellow", "darkred")
p.list[[1]] = ggplot(dev2.meta,aes(x= -UMAP1,y=- UMAP2,color=tf)) + geom_point(shape=16,size=1) + mytheme2 + ggtitle(paste0(tf,' motif-z bgd ATAC-umap')) + scale_colour_gradientn(colours = cols)
p.list[[2]] = ggplot(dev2.meta,aes(x= UMAP1.rna,y= UMAP2.rna,color=tf)) + geom_point(shape=16,size=1) + mytheme2 + ggtitle(paste0(tf,' motif-z bgd RNA-umap')) + scale_colour_gradientn(colours = cols)

library(gridExtra)
png(paste0("Motif.ChromVar.out/ChromVar.with_bgd.",tf,"z.png"),width=14, height=7,units='in', res=200)
do.call("grid.arrange", c(p.list, nrow=1))
dev.off()



