## R4.1.2, for Monocle
library(monocle)
exprs = readRDS('./1_9.monocle.tmp/exprs.rds')
phenoData = readRDS('./1_9.monocle.tmp/phenoData.rds')
featureData = readRDS('./1_9.monocle.tmp/featureData.rds')

# build object
pd <- new("AnnotatedDataFrame", data = phenoData)
fd <- new("AnnotatedDataFrame", data = featureData)
HSMM <- newCellDataSet(exprs,
    phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())

# Estimate size factors and dispersions 
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# Filtering low-quality cells/low-expressed genes
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
    num_cells_expressed >= 10))
# skip cell QC since it's been done in seurat

###### do trajectory
# step 1: feature selection
if(0){
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
              fullModelFormulaStr = "~ Population2 + Phase")
}
# using Nebula DEGs
# from Luria: /net/bmc-lab5/data/kellis/users/xiongxs/Proj_mix/JY/Kras.multiOmic/DEG.0508/dereg.re
# covariates used 
deg = read.table('1_9.monocle.tmp/nebula_ruv.4group.w_cov.pct0.05.tsv.gz',head=T)   
deg.sign = deg[deg$padj<0.01,]
deg.sign = deg.sign[order(abs(deg.sign$logFC),decreasing=T),]
deg.sign = deg.sign[deg.sign$gene %in% expressed_genes,]
rownames(deg.sign) = deg.sign$gene

# ordering_genes <- row.names(subset(diff_test_res, qval < 0.001))  
ordering_genes = rownames(head(deg.sign,2000)) # use the top 2000

HSMM <- setOrderingFilter(HSMM, ordering_genes)

# Trajectory step 2: reduce data dimensionality
HSMM <- reduceDimension(HSMM, max_components = 2,
    method = 'DDRTree')

# do it!
HSMM <- orderCells(HSMM)
pdf('./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DDRTree.pdf')
plot_cell_trajectory(HSMM, color_by = "Population1")
plot_cell_trajectory(HSMM, color_by = "Phase")
dev.off()

saveRDS(HSMM,'./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DDRTree.rds')

##### use all DEGs
ordering_genes = row.names(deg.sign)
HSMM <- setOrderingFilter(HSMM, ordering_genes)

# Trajectory step 2: reduce data dimensionality
HSMM <- reduceDimension(HSMM, max_components = 2,
    method = 'DDRTree')

# do it!
HSMM <- orderCells(HSMM)
pdf('./1_9.monocle.tmp/Monocle.all_nebula_deg.DDRTree.pdf')
plot_cell_trajectory(HSMM, color_by = "Population1")
plot_cell_trajectory(HSMM, color_by = "Phase")
dev.off()
saveRDS(HSMM,'./1_9.monocle.tmp/Monocle.all_nebula_deg.DDRTree.rds')

###### using all DEGs doesn't work well; use 2k instead
HSMM = readRDS('./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DDRTree.rds')

pdf('./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DDRTree.pseudotime.pdf')
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "State") +
    facet_wrap(~State, nrow = 1)
plot_cell_trajectory(HSMM, color_by = "Pseudotime") +
    facet_wrap(~State, nrow = 1)
plot_cell_trajectory(HSMM, color_by = "Pseudotime") +
    facet_wrap(~ Population1, nrow = 1)
dev.off()

HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
          gene_short_name %in% c("Sftpc","Nkx2-1","Hmga2","Cd44")))
cds_subset <- HSMM_filtered[my_genes,]
pdf('./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DDRTree.pseudotime.marker.pdf')
plot_genes_in_pseudotime(cds_subset, color_by = "Population1")
dev.off()

###### do Differential Expression Analysis

# against pseudotime
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
              fullModelFormulaStr = "~ sm.ns(Pseudotime)")
saveRDS(diff_test_res,'./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DEG_pseudotime.with_cov.rds')

sig_gene_names = row.names(head(diff_test_res[order(diff_test_res$qval),],100))# use the top 8000
sig_gene_names = c(sig_gene_names,"Sftpc","Nkx2-1","Hmga2","Cd44")
sig_gene_names = unique(sig_gene_names)

pdf('./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DEG_pseudotime.with_cov.clusters.pdf')
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                num_clusters = 4,
                cores = 1,
                show_rownames = T)
dev.off()

# itersect TF
DEG = row.names(head(diff_test_res[order(diff_test_res$qval),],5000))
tf.merge = read.table('TF.enrichment.2group.txt',head=T)
colnames(tf.merge) = paste0(colnames(tf.merge),'.TF')
tf.merge$TF.name.TF[tf.merge$TF.name.TF=='Nkx21'] = 'Nkx2-1'
tf.merge = tf.merge[tf.merge$group.TF != 'nc',]
tf.merge = tf.merge[tf.merge$TF.name.TF %in% DEG,]
sig_gene_names = tf.merge$TF.name.TF
pdf('./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DEG_pseudotime.with_cov.TF.pdf')
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                num_clusters = 4,
                cores = 1,
                show_rownames = T)
dev.off()

# marker genes; for each group, select top 25
deg = readRDS('1_5.DEG.seurat/Markergene.4groups.wilcox.rds')
groups = unique(deg$cluster)
deg.top = data.frame()
for(group in groups){
 deg.group = head(deg[deg$cluster == group,],25)
 deg.top = rbind(deg.top,deg.group)
}
sig_gene_names = deg.top$gene
sig_gene_names = unique(sig_gene_names[sig_gene_names %in% expressed_genes])
pdf('./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DEG_pseudotime.with_cov.marker.pdf')
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                num_clusters = 4,
                cores = 1,
                show_rownames = T)
dev.off()


# add covariates
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
              fullModelFormulaStr = "~ sm.ns(Pseudotime) + S.Score + G2M.Score")

saveRDS(diff_test_res,'./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DEG_pseudotime.cellstate_cov.rds')

sig_gene_names = row.names(head(diff_test_res[order(diff_test_res$qval),],100))# use the top 8000
sig_gene_names = c(sig_gene_names,"Sftpc","Nkx2-1","Hmga2","Cd44")
sig_gene_names = unique(sig_gene_names)

pdf('./1_9.monocle.tmp/Monocle.top2k.nebula_deg.DEG_pseudotime_cov.with_cov.clusters.pdf')
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                num_clusters = 4,
                cores = 1,
                show_rownames = T)
dev.off()





