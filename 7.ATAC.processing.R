library(ArchR)
set.seed(1)
addArchRThreads(threads = 6)
addArchRGenome("mm10")   ## use mm10

## read in data
#frag.KPY <- "./Data/KPY/atac_fragments.tsv.gz"
frag.YFP <- './Data/YFP/atac_fragments.tsv.gz'
#list.file = c(frag.KPY,frag.YFP)
list.file = c(frag.YFP)
names(list.file) <- c('YFP')

ArrowFiles <- createArrowFiles(
  inputFiles = list.file,
  sampleNames = names(list.file),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles[[2]] = "KPY.arrow"
names(ArrowFiles) = c('YFP','KPY')

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Save-Proj1.all",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)

proj <- proj[proj@cellColData$TSSEnrichment>=6]

ncell=length(proj$cellNames)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", iterations = 3, clusterParams = list(resolution = c(0.2), sampleCells = ncell, n.start = 10), name = "IterativeLSI",varFeatures = 50000,force = TRUE)

proj <- addClusters(input = proj, reducedDims = "IterativeLSI",force = TRUE ,method = "Seurat")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force = TRUE)
proj <- addTSNE(ArchRProj = proj, reducedDims = "IterativeLSI",perplexity = 30, name = "TSNE",force = TRUE)

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
pdf("ArchR.iter3.sampCellAll.sampFeat50000.Seurat.pdf")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
dev.off()

meta<-cbind(proj@embeddings$UMAP$df,proj@embeddings$TSNE$df,proj@cellColData)
write.table(meta,file = "meta_exported.iter3.sampCellAll.sampFeat50000.Seurat.xls",append=FALSE,
		quote= FALSE,sep="\t", eol = "\n",  dec = ".",
		row.names = TRUE, col.names = TRUE, fileEncoding = "")

saveRDS(proj,"ArchR.iter3.sampCellAll.sampFeat50000.Seurat.step3.rds")
