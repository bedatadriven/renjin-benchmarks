################################################################
Single Cell (Seurat, Spatial Inference)
################################################################

All the functions that take place within a cell are performed through proteins.
These proteins are coded within the DNA (Deoxyribonucleic acid) of the cell.
A gene is a sequence of DNA that encodes for a particular protein. In order to
make the necessary proteins, the transcriptional machinary of a cell makes
special copies of the respective genes that can be translated to protein
sequences. These special copies are called messenger RNAs (Ribonucleic acids).

The amount of mRNA produced by a specific gene is used as a surrogate marker for
quantification of the gene activity. With Current highthroughput sequencing
technologies (such as in this case RNAseq), all (fragmented) RNA molecules
from a biological sample are sequenced. These sequences are then matched to
the annotated genome sequence to identify to which gene each sequenced fragment
belongs. More sequenced fragments mapping to a gene in one samples as compared
to another, means that the specific gene had a higher activity. In addition,
single cell sequencing technologies have allowed us to identify the heterogeneous
nature of phenotypically pure populations and identity new sub populations of
cells based on difference in gene expression.

In this workflow developed by `Satija lab <http://www.satijalab.org>`_ imaging data is
used together with single cell transcrptomic data to resolve the spatial position
of the sequenced cells based on their expression profile.


******************************
Packages and Dependencies
******************************
There are 4 packages used in this workflow, which depend
on 68 additional packages (dependencies)

+++++++++++++++
Used packages:
+++++++++++++++

- **Github**: Seurat (satijalab/seurat)

- **CRAN**: XLConnect, rgl, knitr

++++++++++++++++++++++
Package dependencies:
++++++++++++++++++++++

Depends on:
- **CRAN**: ggplot2, reshape2, useful, gridExtra, gplots,
            gdata, XLConnectJars, ROCR, stringr, mixtools,
            lars, fastICA, tsne, Rtsne, fpc,
            ape, VGAM, jackstraw, XLConnect, methods,
            rJava, plyr, scales, gtable, digest,
            MASS, proto, Rcpp, dichromat, labeling,
            RColorBrewer, munsell, colorspace, magrittr, stringi,
            dplyr, assertthat, R6, DBI, lazyeval,
            gtools, KernSmooth, caTools, bitops, boot,
            segmented, diptest, mvtnorm, flexmix, mclust,
            trimcluster, kernlab, robustbase, class, cluster,
            prabclus, nnet, lattice, modeltools, DEoptimR,
            nlme, corpcor, evaluate, markdown, yaml,
            highr, formatR, mime

+++++++++++++++++++++++
Data
+++++++++++++++++++++++

Single cell transcriptome data from `Satija*, Farrell* et al., 2015 <http://doi.org/10.1038/nbt.3192>`_

********************
License
********************
Copyright (c) 2015 Rahul Satija
based on code from `Tutorial - Spatial Inference of single cell data <http://www.satijalab.org/seurat-intro.html>`_
Copyright (c) 2015-2016 BeDataDriven B.V.
License: `GPL version 2 or higher <http://www.gnu.org/licenses/gpl.html>`_
