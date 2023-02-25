# for ArchR, use R3.6.1
library(ArchR)
set.seed(1)
addArchRThreads(threads = 6)
addArchRGenome("mm10")   ## use mm10
pathToMacs2 <- findMacs2()

# read in RNA meta data
rna.meta = readRDS('Kras_P53.2Kgene.RNA.KPY_only.Output.nfeat2000.nPC75.rm_cellcycle.meta.rds')
library(stringr)
rna.meta$cell = str_split_fixed(rownames(rna.meta),":",2)[,2]

# read in ATAC object
proj = readRDS('ArchR.iter3.sampCellAll.sampFeat50000.Seurat.step3.rds')
proj$UMAP1 = proj@embeddings$UMAP$df[,1]
proj$UMAP2 = proj@embeddings$UMAP$df[,2]
atac.meta = proj@cellColData
atac.meta$cell = str_split_fixed(rownames(atac.meta),"#",2)[,2]

# find overlap between ATAC and RNA cells
meta.merge = merge(atac.meta,rna.meta,by='cell',suffix=c('atac','rna'))

meta.merge$Population1[meta.merge$RNA_snn_res.1 %in% c(5)] = 'Hmga2_A'
meta.merge$Population1[meta.merge$RNA_snn_res.1 %in% c(3,8,1,0,10)] = 'Hmga2_B'
meta.merge$Population1[meta.merge$RNA_snn_res.1 %in% c(2,9,6,7)] = 'Nkx2-1_B'
meta.merge$Population1[meta.merge$RNA_snn_res.1 %in% c(4)] = 'Nkx2-1_A'

meta.merge$Population2[meta.merge$Population1 %in% c('Hmga2_A','Hmga2_B')] = 'Hmga2'
meta.merge$Population2[meta.merge$Population1 %in% c('Nkx2-1_A','Nkx2-1_B')] = 'Nkx2-1'

library(ggplot2)
library(BuenColors)
snap.cols = c("grey", "#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
              "#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744",
              "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
              "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
              "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
              "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
              "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
              "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
              "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
              "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
              "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
              "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")
mytheme <- theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),aspect.ratio=1)

j = 20
cls.cols = snap.cols[j:(j+40)]

meta.merge = data.frame(meta.merge)
meta.merge$UMAP1.atac = meta.merge$UMAP1
meta.merge$UMAP2.atac = meta.merge$UMAP2

p.list = list()
p.list2 = list()     # without names on the umap

### ATAC UMAP
# by RNA cluster
pop1.cols = c('aquamarine4','aquamarine1','pink','pink4')
pop2.cols = c('aquamarine1','pink')
tmp<-shuf(meta.merge)
coord = c()
types<-unique(tmp$Population2)
for (i in 1:length(types)){
        ind=which(tmp$Population2 == types[i])
        loc.x = median(tmp[ind,]$UMAP1.atac)
        loc.y = median(tmp[ind,]$UMAP2.atac)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)

tmp$Population2 = factor(tmp$Population2,levels=c('Nkx2-1','Hmga2'))
p.list[[1]] = ggplot(tmp,aes(x= -UMAP1.atac, y=UMAP2.atac,colour=Population2)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("ATAC-UMAP; RNA_cluster")+
        geom_text(data=coord, mapping= aes(x= -x,y=y,label=clust), color="black", size=4,fontface = "bold")+ scale_colour_manual(values=pop2.cols)
p.list2[[1]] = ggplot(tmp,aes(x= -UMAP1.atac, y=UMAP2.atac,colour=Population2)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("ATAC-UMAP; RNA_cluster")+
	scale_colour_manual(values=pop2.cols)

tmp<-shuf(meta.merge)
coord = c()
types<-unique(tmp$Population1)
for (i in 1:length(types)){
        ind=which(tmp$Population1 == types[i])
        loc.x = median(tmp[ind,]$UMAP1.atac)
        loc.y = median(tmp[ind,]$UMAP2.atac)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)

tmp$Population1 = factor(tmp$Population1,levels=c('Nkx2-1_A','Nkx2-1_B','Hmga2_B','Hmga2_A'))
p.list[[2]] = ggplot(tmp,aes(x= -UMAP1.atac, y=UMAP2.atac,colour=Population1)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("ATAC-UMAP; RNA_cluster")+
        geom_text(data=coord, mapping= aes(x= -x,y=y,label=clust), color="black", size=4,fontface = "bold")+ scale_colour_manual(values=pop1.cols)
p.list2[[2]] = ggplot(tmp,aes(x= -UMAP1.atac, y=UMAP2.atac,colour=Population1)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("ATAC-UMAP; RNA_cluster")+
	scale_colour_manual(values=pop1.cols)

# by ATAC cluster
j = 1
cls.cols = snap.cols[j:(j+40)]

tmp<-shuf(meta.merge)
coord = c()
types<-unique(tmp$Clusters)
for (i in 1:length(types)){
        ind=which(tmp$Clusters == types[i])
        loc.x = median(tmp[ind,]$UMAP1.atac)
        loc.y = median(tmp[ind,]$UMAP2.atac)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)

p.list[[3]] = ggplot(tmp,aes(x= -UMAP1.atac, y=UMAP2.atac,colour=Clusters)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("ATAC-UMAP; ATAC_cluster")+
        geom_text(data=coord, mapping= aes(x= -x,y=y,label=clust), color="black", size=4,fontface = "bold")+ scale_colour_manual(values=cls.cols)
p.list2[[3]] = ggplot(tmp,aes(x= -UMAP1.atac, y=UMAP2.atac,colour=Clusters)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("ATAC-UMAP; ATAC_cluster")+
	scale_colour_manual(values=cls.cols)

### RNA UMAP
# by RNA cluster
tmp<-shuf(meta.merge)
coord = c()
types<-unique(tmp$Population2)
for (i in 1:length(types)){
        ind=which(tmp$Population2 == types[i])
        loc.x = median(tmp[ind,]$UMAP1.rna)
        loc.y = median(tmp[ind,]$UMAP2.rna)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)
tmp$Population2 = factor(tmp$Population2,levels=c('Nkx2-1','Hmga2'))
p.list[[4]] = ggplot(tmp,aes(x=UMAP1.rna, y=UMAP2.rna,colour=Population2)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("RNA-UMAP; RNA_cluster")+
        geom_text(data=coord, mapping= aes(x=x,y=y,label=clust), color="black", size=4,fontface = "bold")+ scale_colour_manual(values=pop2.cols)
p.list2[[4]] = ggplot(tmp,aes(x=UMAP1.rna, y=UMAP2.rna,colour=Population2)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("RNA-UMAP; RNA_cluster")+
	scale_colour_manual(values=pop2.cols)

tmp<-shuf(meta.merge)
coord = c()
types<-unique(tmp$Population1)
for (i in 1:length(types)){
        ind=which(tmp$Population1 == types[i])
        loc.x = median(tmp[ind,]$UMAP1.rna)
        loc.y = median(tmp[ind,]$UMAP2.rna)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)
tmp$Population1 = factor(tmp$Population1,levels=c('Nkx2-1_A','Nkx2-1_B','Hmga2_B','Hmga2_A'))
p.list[[5]] = ggplot(tmp,aes(x=UMAP1.rna, y=UMAP2.rna,colour=Population1)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("RNA-UMAP; RNA_cluster")+
        geom_text(data=coord, mapping= aes(x=x,y=y,label=clust), color="black", size=4,fontface = "bold")+ scale_colour_manual(values=pop1.cols)
p.list2[[5]] = ggplot(tmp,aes(x=UMAP1.rna, y=UMAP2.rna,colour=Population1)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("RNA-UMAP; RNA_cluster")+
	scale_colour_manual(values=pop1.cols)

# by ATAC cluster
j = 1
cls.cols = snap.cols[j:(j+40)]

tmp<-shuf(meta.merge)
coord = c()
types<-unique(tmp$Clusters)
for (i in 1:length(types)){
        ind=which(tmp$Clusters == types[i])
        loc.x = median(tmp[ind,]$UMAP1.rna)
        loc.y = median(tmp[ind,]$UMAP2.rna)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)

p.list[[6]] = ggplot(tmp,aes(x=UMAP1.rna, y=UMAP2.rna,colour=Clusters)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("RNA-UMAP; ATAC_cluster")+
        geom_text(data=coord, mapping= aes(x=x,y=y,label=clust), color="black", size=4,fontface = "bold")+ scale_colour_manual(values=cls.cols)
p.list2[[6]] = ggplot(tmp,aes(x=UMAP1.rna, y=UMAP2.rna,colour=Clusters)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("RNA-UMAP; ATAC_cluster")+
	scale_colour_manual(values=cls.cols)

library(gridExtra)
png("ArchR.iter3.sampCellAll.sampFeat50000.Seurat.KPY_only.color_by_RNA_cluster.filtered.rm_cc.grouped.png",width = 14, height = 10,units='in', res=450)
do.call("grid.arrange", c(p.list, nrow=2))
dev.off()

png("ArchR.iter3.sampCellAll.sampFeat50000.Seurat.KPY_only.color_by_RNA_cluster.filtered.rm_cc.grouped.noname.png",width = 14, height = 10,units='in', res=450)
do.call("grid.arrange", c(p.list2, nrow=2))
dev.off()

### RNA UMAP
# by RNA cluster
p.list = list()
j = 1
cls.cols = snap.cols[j:(j+40)]
tmp<-shuf(meta.merge)
coord = c()
types<-unique(tmp$RNA_snn_res.1)
for (i in 1:length(types)){
        ind=which(tmp$RNA_snn_res.1 == types[i])
        loc.x = median(tmp[ind,]$UMAP1.rna)
        loc.y = median(tmp[ind,]$UMAP2.rna)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)

p.list[[1]] = ggplot(tmp,aes(x=UMAP1.rna, y=UMAP2.rna,colour=RNA_snn_res.1)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("RNA-UMAP; RNA_cluster")+
        geom_text(data=coord, mapping= aes(x=x,y=y,label=clust), color="black", size=4,fontface = "bold")+ scale_colour_manual(values=cls.cols)

### ATAC UMAP
# by RNA cluster
j = 1
cls.cols = snap.cols[j:(j+40)]
tmp<-shuf(meta.merge)
coord = c()
types<-unique(tmp$RNA_snn_res.1)
for (i in 1:length(types)){
        ind=which(tmp$RNA_snn_res.1 == types[i])
        loc.x = median(tmp[ind,]$UMAP1.atac)
        loc.y = median(tmp[ind,]$UMAP2.atac)
        coord = rbind(coord,data.frame(clust=types[i],x=loc.x,y=loc.y))
}
coord$clust<-as.factor(coord$clust)

p.list[[2]] = ggplot(tmp,aes(x= -UMAP1.atac, y=UMAP2.atac,colour=RNA_snn_res.1)) + geom_point(shape=16,size=1.2) + mytheme + ggtitle("ATAC-UMAP; RNA_cluster")+
        geom_text(data=coord, mapping= aes(x=-x,y=y,label=clust), color="black", size=4,fontface = "bold")+ scale_colour_manual(values=cls.cols)

library(gridExtra)
png("ArchR.iter3.sampCellAll.sampFeat50000.Seurat.KPY_only.color_by_RNA_cluster.filtered.rm_cc.png",width = 14, height = 7,units='in', res=450)
do.call("grid.arrange", c(p.list, nrow=1))
dev.off()

# separately drawing each RNA cluster
p.list = list()
cls.list = sort(unique(tmp$RNA_snn_res.1))

df = tmp
for(i in 1:length(cls.list)){
 cls = cls.list[[i]]
 df$cls = "Others"
 df$cls[df$RNA_snn_res.1 == cls] = cls
 df = df[order(df$cls,decreasing=T),]
 p.list[[i]] = ggplot(df,aes(x=UMAP1.atac, y=UMAP2.atac,colour=cls)) + geom_point(shape=16,size=1.5) + mytheme + ggtitle(cls) + scale_colour_manual(values=c("darkred","grey"))
}

library(gridExtra)
png("ArchR.iter3.sampCellAll.sampFeat50000.Seurat.KPY_only.color_by_RNA_cluster.filtered.rm_cc.sep_cls.png",width=28, height=21,units='in', res=300)
do.call("grid.arrange", c(p.list, nrow=3))
dev.off()

# flip it
p.list = list()
cls.list = sort(unique(tmp$RNA_snn_res.1))

df = tmp
for(i in 1:length(cls.list)){
 cls = cls.list[[i]]
 df$cls = "Others"
 df$cls[df$RNA_snn_res.1 == cls] = cls
 df = df[order(df$cls,decreasing=T),]
 p.list[[i]] = ggplot(df,aes(x= -UMAP1.atac, y=UMAP2.atac,colour=cls)) + geom_point(shape=16,size=1.5) + mytheme + ggtitle(cls) + scale_colour_manual(values=c("darkred","grey"))
}

library(gridExtra)
png("ArchR.iter3.sampCellAll.sampFeat50000.Seurat.KPY_only.color_by_RNA_cluster.filtered.rm_cc.sep_cls.flip.png",width=28, height=21,units='in', res=300)
do.call("grid.arrange", c(p.list, nrow=3))
dev.off()

## ATAC-RNA matching plot
# this is the tutorial https://r-graph-gallery.com/322-custom-colours-in-sankey-diagram.html
links = data.frame(table(meta.merge[,c('RNA_snn_res.1','Clusters')]))
colnames(links) = c('Cls.RNA','Cls.ATAC','Freq')
links = links[links$Freq>0,]
links$Cls.RNA = paste0('RNA_C',links$Cls.RNA)
links$Cls.ATAC = paste0('ATAC_',links$Cls.ATAC)

nodes <- data.frame(
  name=c(as.character(links$Cls.RNA), as.character(links$Cls.ATAC)) %>% 
    unique()
)

links$IDsource <- match(links$Cls.RNA, nodes$name)-1 
links$IDtarget <- match(links$Cls.ATAC, nodes$name)-1

p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
              Value = "Freq", NodeID = "name")

print(p)



