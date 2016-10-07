
RNAseq DESeq2
=============

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
to another, means that the specific gene had a higher activity.

This workflow uses the exploratory analysis of RNAseq data using many well
stablished Bioconductor packages such as DESeq2, Rsamtools, and
GenomicAlignments as described in 
`Love MI et al., 2015 <http://doi.org/10.12688/f1000research.7035.1>`_


Packages and Dependencies
-------------------------
There are 17 packages used in this workflow, which depend
on 79 additional packages (dependencies).

**Used packages:**

* *Bioconductor*: Rsamtools, GenomicFeatures, GenomicAlignments, BiocParallel, DESeq2, genefilter, AnnotationDbi, org.Hs.eg.db, ReportingTools

* *CRAN*: pheatmap, RColorBrewer, PoiClaClu, ggplot2, Gviz, fission, sva, fission

**Package dependencies:**

* *Bioconductor*: BiocGenerics, Biostrings, GenomeInfoDb, XVector, S4Vectors, rtracklayer, Biobase, biomaRt, IRanges, zlibbioc, GenomicRanges, geneplotter, ggbio, PFAM.db, limma, edgeR, GSEABase, GOstats, VariantAnnotation, BSgenome, OrganismDbi, biovizBase, GO.db

* *CRAN*: bitops, RCurl, RSQLite, DBI, XML, snow, futile.logger, Rcpp, RcppArmadillo, Hmisc, locfit, plyr, scales, reshape2, gtable, digest, MASS, proto, Formula, lattice, gridExtra, nnet, acepack, latticeExtra, cluster, rpart, foreign, survival, annotate, dichromat, labeling, munsell, stringr, xtable, colorspace, magrittr, stringi, Category, R.utils, hwriter, knitr, GGally, graph, Matrix, RBGL, AnnotationForge, evaluate, markdown, yaml, highr, formatR, reshape, mime, matrixStats, mgcv, nlme

Data
-------
RNAseq bam-files from `Solaimani Kartalaei P, (2014) <http://www.doi.org/10.1084/jem.20140767>`_


License
-------
| Copyright (c) 2015 Wolfgang Huber
| based on  `DOI: 10.12688/f1000research.7035.1 <http://www.doi.org/10.12688/f1000research.7035.1>`_
| Copyright (c) 2016 BeDataDriven B.V.
| License: Artistic 2.0

