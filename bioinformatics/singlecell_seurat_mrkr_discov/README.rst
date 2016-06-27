
Single Cell (Seurat, Clustering and marker discovery)
=====================================================
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

In this workflow developed by `Satija lab`_ single
cell transcriptomic data is used in combination of various dimensional reduction
algorithms (such as PCA and t-SNE) to cluster cells. The leading genes in these
clusterings can be used as genetic markers.


Packages and Dependencies
-------------------------
There is 1 package used in this workflow, which depend
on 60 additional packages from CRAN (dependencies)

**Used packages:**

- **Github**: Seurat (satijalab/seurat)

**Package dependencies:**

* *CRAN*: ggplot2, reshape2, useful, gridExtra, gplots, ROCR, stringr, mixtools, lars, fastICA, tsne, Rtsne, fpc, ape, VGAM, jackstraw, XLConnect, plyr, scales, gtable, digest, MASS, proto, Rcpp, dichromat, labeling, RColorBrewer, munsell, colorspace, magrittr, stringi, dplyr, assertthat, R6, DBI, lazyeval, gtools, KernSmooth, caTools, bitops, boot, segmented, diptest, mvtnorm, flexmix, mclust, trimcluster, kernlab, robustbase, class, cluster, prabclus, nnet, lattice, modeltools, DEoptimR, nlme, corpcor, XLConnectJars, rJava

Data:
-------
Single cell transcriptome data from `Pollen et al. 2014 <http://doi.org/10.1038/nbt.2967>`_

License
-------
| Copyright (c) 2015 Rahul Satija
| based on code from `Tutorial - Unsupervised clustering and marker discovery <http://www.satijalab.org/seurat-intro.html>`_
| Copyright (c) 2015-2016 BeDataDriven B.V.
| License: `GPL version 2 or higher`_

.. _Satija lab: http://www.satijalab.org
.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html