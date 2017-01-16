START_WORKFLOW <- as.numeric(Sys.time())
# Copyright (c) 2015 Rahul Satija
# based on code from http://www.satijalab.org/seurat-intro.html
# Copyright (c) 2015-2016 BeDataDriven B.V.
# License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
#
# Single cell transcriptome analysis with Seurat (Spatial inference)
#

# reproducibility
set.seed(8008)
## packages
library(useful)
library(methods)
library(utils)
rVersion <- R.Version()
if ('engine' %in% names(rVersion) && rVersion$engine == "Renjin") {
	library(org.renjin.github.satijalab:Seurat)
} else {
	library(Seurat)
}
library(XLConnect)
library(rgl)
library(knitr)
knit_hooks$set(rgl = hook_rgl)

### Setup Seurat Object
# Read in log-space expression matrix. Has been pre-computed and normalized
# (see manuscript for exact details). The data was run in three batches (zf1,
# zf2, zf3), which is denoted in the column name
zfish.data <- read.table("zdata.matrix.txt", sep = "\t", header = TRUE)

# Due to the high dynamic range in RNA-seq, transform data to log-space. Not
# required for Seurat, but highly recommended
zfish.log <- log(zfish.data + 1)

# Create and setup the Seurat object. Include all genes detected in > 3 cells
# (expression >0.01), and all cells with > 2k genes. Cells will be initially
# assigned to an identity class (grouping) based on the first field (underscore-delimited)
zf.all <- new("seurat", raw.data = zfish.log)
zf.all <- Setup(zf.all, project = "zfish", min.cells = 3, min.genes = 2000, is.expr = 0.01, names.field = 1, names.delim = "_")

zf.all

### Basic exploration of data

# View the expression of genes and basic metrics across samples
VlnPlot(zf.all, c("ACTB2", "VENT", "CHD", "nGene"))

# Plot two cells against each other
#Set do.ident=TRUE to identify outliers by clicking on points (ESC to exit after)

par(mfrow = c(2,2))
CellPlot(zf.all, zf.all@cell.names[1], zf.all@cell.names[2], do.ident = FALSE)
CellPlot(zf.all ,zf.all@cell.names[3], zf.all@cell.names[4], do.ident = FALSE)

#Plot two genes against each other, can do this in limited groups of cells
GenePlot(zf.all, "MIXL1", "OSR1", cex.use = 1)
GenePlot(zf.all, "MIXL1", "SOX3", cell.ids = WhichCells(zf.all, "zf1"), cex.use = 1)

###Identify variable genes across the single cells
#Genes placed into 20 bins based on X-axis (average expression). Y-axis is within-bin z-score of log(Variance/mean).
zf.all <- MeanVarPlot(zf.all, y.cutoff = 2, do.plot = TRUE, x.low.cutoff = 0.25, x.high.cutoff = 7, fxn.x = expMean, fxn.y = logVarDivMean)

#remove markers which primarily define batch
markers.remove <- BatchGene(zf.all, idents.use = c("zf1", "zf2", "zf3"), genes.use = zf.all@var.genes, auc.cutoff = 0.7)
zf.all@var.genes <- zf.all@var.genes[!(zf.all@var.genes %in% markers.remove)]
length(zf.all@var.genes)

###Run PCA, examine results
#To demonstrate, run PCA on all genes (note that the result isn't great, PC1/2 contain little known biology)
zf.all <- PCA(zf.all, pc.genes = rownames(zf.all@data), pcs.print = 3, genes.print = 5)

#Instead, run PCA on zf.all@var.genes (default value for pc.genes).
zf.all <- PCA(zf.all, pcs.print = 3, genes.print = 5)
PCAPlot(zf.all, pt.size = 2)

#If desired, project PCA structure across entire dataset (so all genes are included)
zf.all <- ProjectPCA(zf.all, do.print = FALSE)
#visualize top genes for each PC, use.full selects the projected PCA. See also print.pca
#Explicit thanks to Olga Botvinnik, who showed me this visualization in her Flotilla package
VizPCA(zf.all, pcs.use = 1:3, num.genes = 10, use.full = TRUE, nCol = 3)

###Create new object without EVL cells
#Identify EVL cells (see Manuscript for further detail)
plot(zf.all@pca.rot[ , 1], zf.all@pca.rot[ , 2], pch = 16, xlab = "PC1", ylab = "PC2")

x <- seq(-0.2, 0.2, .01)
lines(x, -x * 0.5 - 0.04, lwd = 2, lty = 2, col = "red")
evl.quant <- zf.all@pca.rot[ , 1] + 2 * zf.all@pca.rot[ , 2] + 0.08
names(evl.quant) <- colnames(zf.all@data)
not.evl <- names(evl.quant[evl.quant > 0])
is.evl <- names(evl.quant[evl.quant < 0])

points(zf.all@pca.rot[is.evl, 1], zf.all@pca.rot[is.evl, 2], pch = 16, col = "red")

#subsetData allows you to create a new Seurat object with a subset of the data
zf <- SubsetData(zf.all, cells.use = not.evl)

###Identify genes to use for building gene expression models
#recalculate a set of variable genes, now that EVL are removed
zf <- MeanVarPlot(zf, y.cutoff = 2, do.plot = FALSE, x.low.cutoff = 1, x.high.cutoff = 7, fxn.x = expMean, fxn.y = logVarDivMean, set.var.genes = TRUE)
markers.remove <- BatchGene(zf, idents.use = c("zf1", "zf2", "zf3"), genes.use = zf@var.genes)
zf@var.genes <- zf@var.genes[!(zf@var.genes %in% markers.remove)]

#redo the PCA on the variable genes.
zf <- PCA(zf, do.print = FALSE)

#Run a 'random' PCA 1,000 times - scrambling a random 2.5% of the data each time
#This enables us to identify statistically significant PCs (in this case, 1:3), and genes with significant PC scores
zf <- JackStraw(zf, num.replicate = 1000, prop.freq = 0.025)
JackStrawPlot(zf)

zf <- ProjectPCA(zf, do.print = FALSE, do.center = FALSE)
genes.sig <- PCASigGenes(zf, pcs.use = c(1, 2, 3), pval.cut = 1e-2, use.full = TRUE)

plot.1 <- PCAPlot(zf, do.return = TRUE)
plot.2 <- PCAPlot(zf, 1, 3, do.return = TRUE)
MultiPlotList(list(plot.1, plot.2), cols = 2)
VizPCA(zf, pcs.use = 1:3, num.genes = 10, nCol = 3)

### Build models of gene expression
# Matrices of gene expression were generated from published in situ stainings, and saved in an Excel file (which eases data entry). So, we import this data and add it to the Seurat object.
# Load in the Excel file.
wb <- loadWorkbook("Spatial_ReferenceMap.xlsx", create = FALSE)
insitu.genes <- getSheets(wb)
insitu.matrix <- data.frame(sapply(1:length(insitu.genes), function(x) as.numeric(as.matrix(wb[x][2:9, 2:9]))))
insitu.genes <- toupper(insitu.genes)
colnames(insitu.matrix) <- (insitu.genes)

# Then, we store this information in the Seurat object.
zf@insitu.matrix <- insitu.matrix[ , insitu.genes]

# Now build models for these insitu genes, and predict robust values
lasso.genes.use <- unique(c(genes.sig,zf@var.genes))

# we will fit models for the landmark genes using the 'structured' genes (with significant PCA scores), and variable genes
#zf <- AddImputedScore(zf, genes.use = lasso.genes.use, genes.fit = insitu.genes, do.print = FALSE, s.use = 40, gram = FALSE)

# Demonstrate the benefit of imputation
#before imputation - MIXL1 and OSR1 should be tightly co-expressed (on t he left)
par(mfrow = c(1, 2))
GenePlot(zf, "MIXL1", "OSR1", col = "black", cex.use = 1)
#after imputation (on the right)
#GenePlot(zf, "MIXL1", "OSR1", use.imputed = TRUE, col = "black", cex.use = 1)


### Build Mixture models of Gene Expression

# The in situ patterns that we use to provide geographical information are scored in a binary on/off format. In order to translate the continuous RNAseq data into this form, we model it as mixtures of 2 normal distributions that represent the on state and off state. We then use this to estimate whether each cell should be considered on or off for each gene.
insitu.genes <- colnames(zf@insitu.matrix)
for(i in rev(insitu.genes)) zf <- FitGeneK(zf, i, do.plot = FALSE, do.k = 2, start.pct = mean(zf@insitu.matrix[ , i]), num.iter = 1)

#show an example mixture model
par(mfrow = c(2, 2))
zf_temp <- FitGeneK(zf, "SOX3", do.plot = TRUE, do.k = 2, start.pct = mean(zf@insitu.matrix[ , "SOX3"]))
zf_temp <- FitGeneK(zf, "OSR1", do.plot = TRUE, do.k = 2, start.pct = mean(zf@insitu.matrix[ , "OSR1"]))
zf_temp <- FitGeneK(zf, "BAMBIA", do.plot = TRUE, do.k = 2, start.pct = mean(zf@insitu.matrix[ , "BAMBIA"]))
zf_temp <- FitGeneK(zf, "SEBOX", do.plot = TRUE, do.k = 2, start.pct = mean(zf@insitu.matrix[ , "SEBOX"]))


### Project each cell into its proper location
#Perform an initial mapping of cells to bins
zf <- InitialMapping(zf)

#Now, perform the quantitative refinement procedure (see manuscript for details)

#first identify the genes to use for the refinement (top 3 genes in both directions for PC1-3)
genes.use <- PCTopGenes(zf, pc.use = 1:3, num.genes = 6, use.full = TRUE, do.balanced = TRUE)

#impute values for these genes if needed
new.imputed <- genes.use[!genes.use %in% rownames(zf@imputed)]
lasso.genes.use <- unique(c(zf@var.genes, PCASigGenes(zf, pcs.use = c(1, 2, 3), pval.cut = 1e-2, use.full = TRUE, )))
zf <- AddImputedScore(zf, genes.use = lasso.genes.use, genes.fit = new.imputed, do.print = FALSE, s.use = 40, gram = FALSE)

#refine the mapping with quantitative models that also consider gene covariance
zf <- RefinedMapping(zf, genes.use)


### Brief analysis of the mapped cells
#view the raw mapping probabilities
corner(zf@final.prob)

#calculate centroids for each cell
zf.centroids <- GetCentroids(zf)
colnames(zf.centroids) <- c("bin.tier", "bin.dv")
corner(zf.centroids)

#add this info into the Seurat object
#note that you can do this for any metaData info (alignment rate, external measurements for each cell, etc.)
zf <- addMetaData(zf, zf.centroids)

#We can use genePlot to examine relationships between different types of variables
#Note that individual genes (or PCs) display strong, but noisy and non-linear, relationships with Seurat's mapping positions along both axes
par(mfrow = c(2, 2))
GenePlot(zf, "OSR1", "bin.tier", cex.use = 1)
GenePlot(zf, "BAMBIA", "bin.dv", cex.use = 1)
GenePlot(zf, "PC1", "bin.tier", cex.use = 1)
GenePlot(zf, "PC2", "bin.dv", cex.use = 1)


##### Beautiful Zebrafish visaulizations courtesy of Jeff Farrell (jfarrell@g.harvard.edu)
##### Before running these functions (those starting with zf.), or loading the rgl library, please make sure an X11 client (i.e. XQuartz) is installed
###Draw inferred in situ patterns for a few known genes

zf.insitu.lateral(zf, "GSC", label = FALSE)
zf.insitu.lateral(zf, "SOX3", label = FALSE)
zf.insitu.lateral(zf, "VED", label = FALSE)


###Draw inferred in situ patterns for a few new genes
zf.insitu.lateral(zf, "RIPPLY1", label = FALSE)
zf.insitu.lateral(zf, "DUSP4", label = FALSE)
zf.insitu.lateral(zf, "ETS2", label = FALSE)


###Focused analysis of cells near the margin
zf.centroids <- get.centroids(zf)
colnames(zf.centroids) <- c("bin.tier", "bin.dv")
margin.cells <- rownames(subset(zf.centroids, bin.tier > 5))
zf.margin <- subsetData(zf, cells.use = margin.cells)
#note, could replace the last three lines with: zf.margin=subsetData(zf,subset.name = "bin.tier",accept.low = 5)

zf.margin <- PCA(zf.margin, do.print = FALSE)
zf.margin <- JackStraw(zf.margin, num.replicate = 1000, prop.freq = 0.025)
zf.margin <- ProjectPCA(zf.margin, pcs.print = 3, genes.print = 8, do.center = FALSE)

kmeans.genes <- PCASigGenes(zf.margin, pcs.use = 1:3, pval.cut = 1e-3)
zf.margin <- doKMeans(zf.margin, genes.use = kmeans.genes, k.genes = 7, k.cells = 8, pc.row.order = 3, pc.col.order = 3, rev.pc.order = TRUE, cexRow = 0.5, cexCol = 0.5)
#We can see high Gsc in cluster 5, and high Sox32 in cluster 1.
#Note that clusters (cells and genes) are ordered by average PC3 score, just for visualization.

endo.cells <- WhichCells(zf.margin, 1)
plate.cells <- WhichCells(zf.margin, 5)
zf.margin <- set.ident(zf.margin, zf.margin@cell.names, "All")
zf.margin <- set.ident(zf.margin, endo.cells, "Endoderm")
zf.margin <- set.ident(zf.margin, plate.cells, "PCP")

pop.cols <- c("#999999","#1B75BB","#37B34A")

#plot PCA, coloring cells by their status/ID (stored in zf@ident)
PCAPlot(zf.margin, 1, 2, pt.size = 3.5, cols.use = pop.cols)
pc.13 <- PCAPlot(zf.margin, 1, 3, pt.size = 3.5, cols.use = pop.cols)

#gene-gene plot, again, coloring cells by their stat/ID
genePlot(zf.margin, "GSC", "FRZB", col.use = pop.cols)


###Define markers for Endodermal and PCP
#Uses DE test from McDavid et al, Bioinformatics, 2013
#since both the PCP and the Endoderm are on the lowest band of the margin, we will use Tier 8 cells as a control population
tier8.cells <- rownames(subset(zf.centroids, bin.tier >= 7))

#Note that find.markers accepts either a list of cells (plate.cells, endo.cells) or an identity class (e.g. "PCP, Endoderm")
#If the second argument is left blank (currently, tier8.cells), Seurat will test against all other cells in the object.
pcp.markers <- find.markers(zf.margin, plate.cells, tier8.cells, test.use = "bimod", thresh.use = 2)
endoderm.markers <- find.markers(zf.margin, "Endoderm", tier8.cells, test.use = "bimod", thresh.use = 1)

subset(pcp.markers, avg_diff > 0 & p_val < 0.01)[1:15, ]

subset(endoderm.markers, avg_diff > 0 & p_val < 0.01)[1:15, ]

# Draw violin plots of known and new markers
VlnPlot(zf.margin, c("GSC", "SOX32", "RIPPLY1", "FRZB"), cols.use = pop.cols, nCol = 2)


###Localize cellular populati ons
#prechordal plate
zf.cells.render(zf, plate.cells, do.rotate = FALSE, radius.use = 0.0625, col.use = pop.cols[3], do.new = TRUE)

#endoderm progenitors
zf.cells.render(zf, endo.cells, do.rotate = FALSE, radius.use = 0.0625, col.use = pop.cols[2], do.new = TRUE)

END_WORKFLOW <- as.numeric(Sys.time())
TOTAL_TIME <- END_WORKFLOW - START_WORKFLOW
print(TOTAL_TIME)
write(TOTAL_TIME,file="TIMINGS",append=TRUE)
