# R4.1.1
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

ldat <- ReadVelocity(file = "output_KPY/gex_possorted_bam_OL7OP.loom")

# use the RNA cells only
Kras.rna = readRDS('../Kras_P53.2Kgene.RNA.KPY_only.Output.nfeat2000.nPC75.rm_cellcycle.assigned.rds')
rna.meta = Kras.rna@meta.data
rna.meta$cell = rownames(rna.meta)
rna.meta$cell = gsub('KPY','gex_possorted_bam_OL7OP',rna.meta$cell)   # to match velocity cell id
rna.meta$cell = gsub('-1','x',rna.meta$cell)   # to match velocity cell id

# remove the cell cycle genes!
cc = read.table("/n/data2/bch/hemonc/ckim/JINGYUN/Reference/mouse/cell.cycle/Mus_musculus.parse.csv",head=T)
ldat$spliced = ldat$spliced[!(rownames(ldat$spliced) %in% cc$genename),as.character(rna.meta$cell)]
ldat$unspliced = ldat$unspliced[!(rownames(ldat$unspliced) %in% cc$genename),as.character(rna.meta$cell)]
ldat$ambiguous = ldat$ambiguous[!(rownames(ldat$ambiguous) %in% cc$genename),as.character(rna.meta$cell)]

# do the work!
kras <- as.Seurat(x = ldat)
kras <- SCTransform(object = kras, assay = "spliced")
kras <- RunPCA(object = kras, verbose = FALSE)
kras <- FindNeighbors(object = kras, dims = 1:30)
kras <- FindClusters(object = kras)
kras <- RunUMAP(object = kras, dims = 1:30)
kras <- RunVelocity(object = kras, deltaT = 1, kCells = 25, fit.quantile = 0.02)
# saveRDS(kras,'Kras.velocity.rds')
save.image('Kras.velocity.SeuratWrappers.rm_cc.RData')

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = kras)))
names(x = ident.colors) <- levels(x = kras)
cell.colors <- ident.colors[Idents(object = kras)]
names(x = cell.colors) <- colnames(x = kras)

pdf('Kras.velocity.SeuratWrappers.rm_cc.pdf')
show.velocity.on.embedding.cor(emb = Embeddings(object = kras, reduction = "umap"), vel = Tool(object = kras, 
    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1,n.cores = 1)
dev.off()

# export umap coordinates
meta <- kras@meta.data
umapCoord <- as.data.frame(Embeddings(object = kras[["umap"]]))
all(rownames(umapCoord) == rownames(meta))
meta$UMAP1.velocity = umapCoord$UMAP_1
meta$UMAP2.velocity = umapCoord$UMAP_2
saveRDS(meta,'Kras.velocity.SeuratWrappers.rm_cc.UMAP.rds')

# plot cell clusters on the velocity umap
Kras.rna = readRDS('../Kras_P53.2Kgene.RNA.KPY_only.Output.nfeat2000.nPC75.rm_cellcycle.assigned.rds')
rna.meta = Kras.rna@meta.data

library(stringr)
rna.meta$cell = str_split_fixed(rownames(rna.meta),":",2)[,2]
rna.meta$cell = gsub('-1','',rna.meta$cell)   # to match velocity cell id

meta$cell = str_split_fixed(rownames(meta),":",2)[,2]
meta$cell = gsub('x','',meta$cell)   # to match velocity cell id
meta.cut = meta[,c('cell','UMAP1.velocity','UMAP2.velocity')]

meta.merge = merge(rna.meta,meta.cut,by='cell')

pop1.cols = c('aquamarine4','aquamarine1','pink','pink4')
library(BuenColors)
tmp<-shuf(meta.merge)
coord = c()
types<-unique(tmp$Population1)
for (i in 1:length(types)){
        ind=which(tmp$Population1 == types[i])
        loc.x = median(tmp[ind,]$UMAP1.velocity)
        loc.y = median(tmp[ind,]$UMAP2.velocity)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)

p.list = list()
tmp$Population1 = factor(tmp$Population1,levels=c('Nkx2-1_A','Nkx2-1_B','Hmga2_B','Hmga2_A'))
p.list[[1]] = ggplot(tmp,aes(x=UMAP1.velocity, y=UMAP2.velocity,colour=Population1)) + geom_point(shape=16,size=1.2) + theme_classic() + ggtitle("velocity-UMAP; RNA_cluster")+
        geom_text(data=coord, mapping= aes(x=x,y=y,label=clust), color="black", size=4,fontface = "bold")+ scale_colour_manual(values=pop1.cols)
p.list[[2]] = ggplot(tmp,aes(x=UMAP1.velocity, y=UMAP2.velocity,colour=Population1)) + geom_point(shape=16,size=1.2) + theme_classic() + ggtitle("velocity-UMAP; RNA_cluster")+
	scale_colour_manual(values=pop1.cols)

library(gridExtra)
png("Velocity.UMAP.by.RNA_clusters.rm_cc.png",width = 14, height = 7,units='in', res=400)
do.call("grid.arrange", c(p.list, nrow=1))
dev.off()







