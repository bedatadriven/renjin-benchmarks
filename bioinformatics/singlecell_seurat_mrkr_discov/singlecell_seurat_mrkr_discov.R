START_WORKFLOW <- as.numeric(Sys.time())
# Copyright (c) 2015 Rahul Satija
# based on code from http://www.satijalab.org/seurat-intro.html
# Copyright (c) 2015-2016 BeDataDriven B.V.
# License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
#
# Single cell transcriptome analysis with Seurat (v1)
#

# reproducibility
set.seed(8008)
## packages
rVersion <- R.Version()
if ('engine' %in% names(rVersion) && rVersion$engine == "Renjin") {
	library(org.renjin.github.satijalab:Seurat)
} else {
	library(Seurat)
}


# Load in dataset from Pollen et al. 2014
nbt.data <- read.table("HiSeq301-RSEM-linear_values.txt", sep = "\t", header = TRUE, row.names = 1)

# transform data to log scale
nbt.data <- log(nbt.data + 1)

# Look at the transformed data matrix
corner(nbt.data)

dim(nbt.data)

nbt <- new("seurat", raw.data = nbt.data)

# Take all genes in > 3 cells, all cells with > 1k genes, use an expression threshold of 1
# Cell type is encoded in the second _ field, will be stored in nbt@ident and also placed in the "orig.ident" field of object@data.info
nbt <- setup(nbt, project = "NBT", min.cells = 3, names.field = 2, names.delim = "_", min.genes = 1000, is.expr =  1, )

nbt

# Look at some canonical marker genes and metrics
vlnPlot(nbt, c("DPPA4","GATA1","BMP3","nGene"))

# Plot two cells against each other
# Set do.ident=TRUE to identify outliers by clicking on points (ESC to exit after)

par(mfrow = c(2, 2))
cellPlot(nbt, nbt@cell.names[1], nbt@cell.names[2], do.ident = FALSE)
cellPlot(nbt, nbt@cell.names[3], nbt@cell.names[4], do.ident = FALSE)

# Plot two genes against each other, can do this in limited groups of cells
genePlot(nbt, "DLX1", "DLX2", cex.use = 1)
genePlot(nbt, "DLX1", "DLX2", cell.ids = which.cells(nbt,"GW16"), cex.use = 1)

# Identify variable genes across the single cells
# Genes placed into 20 bins based on X-axis (average expression). Y-axis is within-bin z-score of log(Variance/mean).
# Running this sets object@var.genes by default
nbt <- mean.var.plot(nbt, y.cutoff = 2, x.low.cutoff = 2, fxn.x = expMean, fxn.y = logVarDivMean)

length(nbt@var.genes)

# Perform linear dimensional reduction (PCA)
# Run a PCA using the variable genes as input (to change the input gene set, use the pc.genes argument)
# For example, try running PCA on all genes - i.e. , nbt=pca(nbt,pc.genes=rownames(nbt@data)) - which does not perform as well
nbt <- pca(nbt, do.print = FALSE)
pca.plot(nbt, 1, 2, pt.size = 2)

# Examine  and visualize PCA results a few different ways
print.pca(nbt, 1)

viz.pca(nbt, 1:2)

# Draw a heatmap where both cells and genes are ordered by PCA score
# Options to explore include do.balanced (show equal # of genes with +/- PC scores), and use.full (after projection, see below)
pcHeatmap(nbt, pc.use = 1, do.balanced = FALSE)

# Determine statistically significant principal components
# Do 200 random samplings to find significant genes, each time randomly permute 1% of genes
# This returns a 'p-value' for each gene in each PC, based on how likely the gene/PC score woud have been observed by chance
# Note that in this case we get the same result with 200 or 1000 samplings, so we do 200 here for expediency
nbt <- jackStraw(nbt, num.replicate = 200, do.print = FALSE)

# The jackStraw plot compares the distribution of P-values for each PC with a uniform distribution (dashed line)
# 'Significant' PCs will have a strong enrichment of genes with low p-values (solid curve above dashed line)
# In this case PC1-9 are strongly significant
jackStrawPlot(nbt, PCs = 1:12)

# Optional step, grow gene list through PCA projection
# Previous analysis was performed on < 400 variable genes. To identify a larger gene set that may drive
# biological differences, but did not pass our mean/variability thresholds, we first calculate PCA scores
# for all genes (PCA projection)
nbt <- project.pca(nbt, do.print = FALSE)

# Visualize the full projected PCA, which now includes new genes which were not previously included (use.full=TRUE)
pcHeatmap(nbt, pc.use = 1, use.full = TRUE, do.balanced = TRUE, remove.key = TRUE)

# Choose significant genes for PC1-9, allow each PC to contribute a max of 200 genes (to avoid one PC swamping the analysis)
nbt.sig.genes <- pca.sig.genes(nbt, 1:9, pval.cut = 1e-5, max.per.pc = 200)
length(nbt.sig.genes)

# Now redo the PCA analysis with the new gene list
nbt <- pca(nbt, pc.genes = nbt.sig.genes, do.print = FALSE)

# Redo random sampling, PCs 1-11 are significant (additional PCs are borderline, but signal disappears with 1k replicates, though we do 200 here for expediency)
nbt <- jackStraw(nbt, num.replicate = 200, do.print = FALSE)
jackStrawPlot(nbt, PCs = 1:15)

# When we run tSNE using our 11 significant PCs as input (spectral tSNE), we get distinct point clouds
# All cell types cluster togther, with the exception of GW16/21 (gestational week 16/21 neurons)
# These cluster separately by maturation state and subtype (see part 2)
nbt <- run_tsne(nbt, dims.use = 1:11, max_iter = 2000)

tsne.plot(nbt, pt.size = 1)

# You can also run tSNE at the gene level (set genes.use instead of dims.use), but this fails to acheive clean separation between cell types
# Note that, depending on the perplexity parameter, this approach may be sub-optimal for finding very rare subtypes (~5 cells or less)
# in this case, Seurat could be complementary with other clustering approaches.

# Cluster the cells
# Density cluster the tSNE map - note that the G.use parameter is the density parameter for the clustering - lower G.use to get finer settings
# Cells which are 'unassigned' are put in cluster 1 - though in this case there are none
# Assigned cluster will be placed in the 'DBClust.ident' field of nbt@data.info. Putting set.ident=TRUE means that the assigned clusters will also be stored in nbt@ident
nbt <- DBclust_dimension(nbt, 1, 2, reduction.use = "tsne", G.use = 8, set.ident = TRUE)
tsne.plot(nbt, pt.size = 1)

# Build a phylogenetic tree, and rename/reorder cluster names according to their position on the tree
# See help for details on tree building strategy
# This gives closely related clusters similar cluster IDs, which is occasionally useful for visualization later on
# Assigned cluster will be placed in the 'tree.ident' field of nbt@data.info, and also stored in nbt@ident
nbt <- buildClusterTree(nbt, do.reorder = TRUE, reorder.numeric = TRUE, pcs.use = 1:11)

# View the t-SNE plot with the new labels, in a slightly different format
tsne.plot(nbt, do.label = TRUE, label.pt.size = 0.5)

# Pulling data from a Seurat object
# First, we introduce the fetch.data function, a very useful way to pull information from the dataset.
# Essentially it is a wrapper to pull from nbt@data, nbt@ident, nbt@pca.rot, nbt@data.info, etc...

# Pulls the identity class (cluster ID), PC1 scores, # of genes, original identity class (parsed from the cell name), gene expression levels for SOX1 and ACTB
# Info is pulled for all cells, but displayed for the last 5
# Note that the GW16 neurons fall into three clusters (1-3), with variable expression of PAX6 and DLX2 (see below for more)
my.data <- fetch.data(nbt, c("ident", "PC1", "nGene", "orig.ident", "PAX6", "DLX2", "ACTB"))
tail(my.data, 5)

# Can use the which.cells function to pull the names of cells that were experimentally marked as iPS cells, and find that they have been placed in cluster 7
ips.cells <- which.cells(nbt, "iPS", id = "orig.ident")
my.data <- fetch.data(nbt, "ident", cells.use = ips.cells)
head(my.data, 5)

# You can also switch easily switch the cell's identity (for example, going back to the original annotation)
nbt <- set.all.ident(nbt, "orig.ident")
# And switch back - to the cluster ID defined by the tree building
nbt <- set.all.ident(nbt, "tree.ident")
# Also see set.ident for changing identity classes for only a subset of cells

# Finding differentially expressed genes (cluster biomarkers)
#find all markers of cluster 8
#thresh.use speeds things up (increase value to increase speed) by only testing genes whose average expression is > thresh.use between cluster
#Note that Seurat finds both positive and negative markers (avg_diff either >0 or <0)
ips.markers <- find.markers(nbt, 7, thresh.use = 2)
print(head(ips.markers, 5))

# note that Seurat has four tests for differential expression:
# ROC test ("roc"), t-test ("t"), LRT test based on zero-inflated data ("bimod", default), LRT test based on tobit-censoring models ("tobit")
# The ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect). Though not a statistical test, it is often very useful for finding clean markers.
ips.markers <- find.markers(nbt, 7, thresh.use = 2, test.use = "roc")

#Visualize these new markers with a violin plot
vlnPlot(nbt, c("CRABP1", "LINC-ROR"))

#plots a correlation analysis of gene/gene (ie. 'FACS' plot - cells colored by cluster number)
genePlot(nbt, "CRABP1", "LINC-ROR")

# Neuronal cells in the dataset (GW represents gestational week) cluster into three groups (1-3) on the phylogenetic tree, let's explore these grouos
plotClusterTree(nbt)

# Names of cells in C1
which.cells(nbt, 1)

# Find markers of C1 (compared to a background of all other cells)
c1.markers <- find.markers(nbt, 1, thresh.use = 2, test.use = "roc")

# Find markers of C1 (compared to a background of clusters 2 and 3)
# As noted in Pollen et al., Many of these markers (DLX6, DLX2, GAD2, etc.) are known markers of inhibitory interneurons
c1.markers <- find.markers(nbt, 1, c(2, 3), thresh.use = 2,test.use = "roc")
print(head(c1.markers, 10))

# Find markers that distinguish clusters 2 and 3 (markers with +avg_diff distinguish C2, markers with -avg_diff characterize C3)
# As noted in Pollen et al., these markers are consistent with different neuronal maturation states
c2.markers <- find.markers(nbt, 2, 3, thresh.use = 2,test.use = "roc")
print(head(c2.markers, 10))

# Visualizing markers of different clusters
par(mfrow = c(1, 2))

genes.viz <- c("DLX6-AS1", "NEUROD2", "GAD2", "TFAP2C")
vlnPlot(nbt, genes.viz)

# Note that you can also view the violin plot grouping the cells by the original annotation - where these markers appear highly heterogeneous
# vlnPlot(nbt, genes.viz,group.by="orig.ident")

# View which cells express a given gene (red is high expression) on either a tSNE plot
feature.plot(nbt, genes.viz, pt.size = 1)

# Can also visualize on a PCA plot using:
# feature.plot(nbt,genes.viz,pt.size = 1,reduction.use = "pca")

# A different visualization, where the size of the dot is the percentage of cells expressing, and the color is the average expression level (green is high)
# Visualize markers that define cluster 2 vs cluster 3
dot.plot(nbt, genes.plot = rownames(c2.markers)[1:10], cex.use = 4 )

# Find markers for all clusters, set do.print=TRUE to output progress
markers.all <- find_all_markers(nbt, thresh.test = 3, test.use = "roc", do.print = TRUE)

head(markers.all)

# Select markers for plotting on a Heatmap (positive markers with high discriminatory power)
markers.use <- subset(markers.all, avg_diff > 0 & power > 0.8)$gene
markers.use.neuronal <- subset(markers.all, avg_diff > 0 & power > 0.8 & (cluster < 4 | cluster == 8))$gene

# Draw a heatmap of all cells for these marker genes
# Gene names are a bit too small to see in the tutorial, but you can blow up the heatmap in R
doHeatMap(nbt, genes.use = markers.use, slim.col.label = TRUE, remove.key = TRUE, cexRow = 0.1)

# Just for the neuronal cells
doHeatMap(nbt, genes.use = markers.use.neuronal, slim.col.label = TRUE, remove.key = TRUE, cells.use = which.cells(nbt, c(1,2,3,8)))

# A few more tricks
# Since we have placed cells into clusters, we can look at and compare the average expression within a cluster
# For example, we can compare cluster 7 (ES cells), and cluster 8 (NPCs)
nbt.avg <- average.expression(nbt)
colnames(nbt.avg) <- paste("c", colnames(nbt.avg), sep = "")
plot(nbt.avg[ , "c7"], nbt.avg[ , "c8"], pch = 16, cex = 0.8)

# uncomment the next line to identify individual genes
# identify(nbt.avg[,"c7"],nbt.avg[,"c8"],labels = rownames(nbt.avg))

# What if we want to define markers of a clade, instead of an individual cluster? Re-examine the cluster tree above
# Find markers of the largest tree split (defined by node 12)
# Markers with a avg_diff>0 mark the left branch, avg_diff<0 mark the right branch
node12.markers <- find.markers.node(nbt, 12, thresh.use = 2, test.use = "roc")
head(node12.markers, 5)

# These are markers which are shared between clusters (1:3).
vlnPlot(nbt, c("NNAT", "TUBA1A", "DCX", "FXYD5"))

# Find markers of a middle split (node 14). These separate clusters 4:6 from clusters 7:11
node14.markers <- find.markers.node(nbt, 14, thresh.use = 2, test.use = "roc")
vlnPlot(nbt, c("LCP1", "NGFRAP1", "PTPRF", "CNN3"))

# Rename a cluster - for example let's rename cluster 1 to be Interneurons
# Now 1 will be replaced with Interneurons for future plots
nbt <- rename.ident(nbt, 1, "Interneurons")

# Final visualization! Splits a 'feature plot' into clusters, very useful for seeing lots of info across many clusters
# Applied here to PC scores. You can see that cluster 3 is uniquely marked by PC8, and cluster 8 is uniquely marked by PC11, etc.
pcs.plot <- paste("PC", 1:11, sep = "")
feature.heatmap(nbt, pcs.plot, cols.use = heat.colors(10), pt.size = 2)

END_WORKFLOW <- as.numeric(Sys.time())
TOTAL_TIME <- END_WORKFLOW - START_WORKFLOW
print(TOTAL_TIME)
write(TOTAL_TIME, file = "TIMINGS", append = TRUE)
