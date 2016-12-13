
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

.. graphviz::
   :caption: Diagram of RNAseq analysis using DESeq2.

   digraph GRAPH_ID {
    fontname="sans-serif";
    penwidth="0.1";
    compound="true";
    edge [comment="Wildcard edge", 
          fontname="sans-serif", 
          fontsize=10, 
          colorscheme="blues3", 
          color=2, 
          fontcolor=3];
    node [fontname="serif", 
          fontsize=10, 
          fillcolor="1", 
          colorscheme="blues4", 
          color="2", 
          fontcolor="4", 
          style="filled"];
    subgraph "cluster0" {
        label="Read raw data";
        edge [comment="subgraph edge wildcard"];
        node [comment="subgraph node wildcard"];
        "bams" [shape="invhouse", 
                label="BAM Files"];
        "bamfilelist" [label="BamFileList \n (yieldsize=2e6)"];
        "gtf" [shape="invhouse", 
               label="GTF"];
        "c0_mtxdb" [label="makeTxDbFromGFF \n (format=GTF) \n\n exonsBy \n (by=gene)"];
        "se" [group=g1, 
              label="summarizeOverlaps()"];
        "dds" [colorscheme="ylorrd3", 
               label="DESeqDataSet \n (design=~cell+dex)"];
        "c0_countdata" [label="assay() \n colData()"];
        "ddsMat" [shape="box", 
                  label="DESeqDataSetFromMatrix( \n countdata, colData, \n design=~cell+dex)"];
        "bams" -> "bamfilelist";
        "gtf" -> "c0_mtxdb";
        "c0_mtxdb" -> "se";
        "bamfilelist" -> "se";
        "se" -> "dds";
        "se" -> "c0_countdata";
        "c0_countdata" -> "ddsMat";
        "dds" -> "ddsMat";
    }

    subgraph "cluster1" {
        label="Exploratory analysis";
        edge [comment="subgraph edge wildcard"];
        node [comment="subgraph node wildcard"];
        "c1_dds" [colorscheme="ylorrd3", 
                  label="DESeqDataSet"];
        "plot_rlog" [shape="box", 
                     label="plot() \n plotPCA() \n MDS plot"];
        "dds_eSF" [label="estimateSizeFactors()"];
        "dist_dds_eSF" [label="dist()"];
        "pheatmap_dds" [shape="box", 
                        label="pheatmap() \n MDS Plot"];
        "dds_rlog" [label="rlog()"];
        "dds_cnt_poisdist" [label="PoissonDistance()"];
        "pheatmap_cnt_pois" [shape="box", 
                             label="pheatmap() \n MDS plot"];
        "c1_dds" -> "dds_eSF";
        "dds_eSF" -> "dist_dds_eSF";
        "dist_dds_eSF" -> "pheatmap_dds";
        "c1_dds" -> "dds_rlog"  [ltail="cluster0"];
        "dds_rlog" -> "plot_rlog";
        "dds_eSF" -> "dds_cnt_poisdist";
        "dds_cnt_poisdist" -> "pheatmap_cnt_pois";
    }

    subgraph "cluster2" {
        label="Differential Expression";
        edge [comment="subgraph edge wildcard"];
        node [comment="subgraph node wildcard"];
        "c2_dds" [label="DESeqDataSet", 
                  colorscheme="ylorrd3"];
        "c2_rlog" [label="rlog()"];
        "c2_assay" [label="↳ assay() \n ↳ rowVars() \n ↳ order() \n ↳ head()"];
        "c2_pheatmap" [shape="box", 
                       label="pheatmap()"];
        "dds_eSF_deseq" [label="DESeq() \n results() \n quantile() \n cut()"];
        "c2_plotCounts" [shape="box", 
                         label="plotCounts()\n hist()\n plotMA() \n barplot()"];
        "c2_res_granges" [label="results(format=GRanges) \n mapIds() \n strand()"];
        "c2_init_tracks" [label="GenomeAxisTrack() \n AnnotationTrack() \n DataTrack()"];
        "c2_plot_track" [shape="box", 
                         label="plotTrack()"];
        "c2_counts" [label="counts(normalized=TRUE) \n svaseq(n.sv=2)"];
        "c2_design" [label="design(~SV1+SV2+dex) \n DESeq()"];
        "c2_stripchart" [shape="box", 
                         label="stripchart()"];
        "c2_dds" -> "dds_eSF_deseq";
        "dds_eSF_deseq" -> "c2_plotCounts";
        "c2_dds" -> "c2_res_granges";
        "c2_res_granges" -> "c2_init_tracks";
        "c2_init_tracks" -> "c2_plot_track";
        "c2_dds" -> "c2_rlog";
        "c2_rlog" -> "c2_assay";
        "c2_assay" -> "c2_pheatmap";
        "c2_dds" -> "c2_counts";
        "c2_counts" -> "c2_design";
        "c2_counts" -> "c2_stripchart";
    }

}

	

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

