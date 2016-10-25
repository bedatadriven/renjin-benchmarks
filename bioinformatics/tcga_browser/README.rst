
TCGA browser
============
All cancers are caused by mutations in DNA. `The Cancer Genome Atlas (TCGA)`_
is a project that begun in 2005 and aims to catalog the mutations that cause
cancer in patients. DNA sequence, gene expression profile, and other relevant
patient information are collected for identification of risk factors and better
understanding the disease mechanism.

First, gene expression profiles (the amount of RNA transcripts produced by each
gene), DNA mutations, and patient information which have been previously
retrieved from TCGA repository are loaded. These data are stored as data frames
using 'data.table' package which is one of the most efficient (memory- and
speed wise) packages for handeling large data frames. Significant differences
in gene activity (expression) are calculated using the 'limma' package. This
type of analysis result often in large lists of differentially expressed, each
of these genes is in turn is involved in a few and sometimes many biological
processes. In analysis of gene expression data, it is a common practice to
statistically test whether there are biological processes/pathways whose genes
are significantly overrepresented using a gene set analysis (GSA) method. This
gives us a view at the level of biological processes. In this workflow 'gage'
package has been used which uses Generally Applicable Gene-set Enrichment
`GAGE method`_. Gene transcripts are translated to proteins and proteins can
interact and inhibit or activate each others function or expression level.
Information about protein-protein interactions are continuously added to
`STRING database`_ which is accessible through 'STRINGdb' package and used to
show the interaction between the genes with significant differential
expression.

.. graphviz::
   :caption: TCGA browser shiny app modules.

   digraph TCGA_browser {
      TCGA_mrna         [group = g1;  shape = invhouse, label = "Data (TCGA)\nRNAseq"];
      TCGA_patient      [group = g1;  shape = invhouse, label = "Data (TCGA)\nPhenotype"];
      TCGA_exome        [group = g1;  shape = invhouse, label = "Data (TCGA)\nExome"];
      Data_string       [group = g1;  shape = invhouse; lable = "STRINGdb"];
      Data_kegg         [group = g1;  shape = invhouse; label = "KEGG"];
      pp                [shape = box; label = "preprocess    \nRNAseq"];
      pp_sub            [shape = box; label = "restructure data     \nsubset()\nmerge()\nmelt()"];
      an_mrna_deg       [shape = box; label = "Differential Expression     \n(limma)"];
      an_mrna_hm        [shape = box; label = "Plot\nheatmap.2()"];
      an_mrna_gage      [shape = box; label = "Geneset Analysis     \ngage()"];
      an_mrna_string    [shape = box; label = "Protein interaction network       \nSTRINGdb"];
      an_string_plot    [shape = box; label = "Protein Network Plot     "];
      an_surv           [shape = box; label = "Simple survival analysis     \nSurv()\nsurvdiff()\nsurvfit()"];
      an_surv_plot      [shape = box; label = "Kaplan-Meier plot     \nplot()"];
      pp_exome          [shape = box; label = "Integrate     \nExome & Phenotype\ndata.table()\nmelt()"];
      an_exome          [shape = box; label = "Plot     \nggvis\nggplot2"];
      TCGA_mrna     ->      pp          ->   pp_sub;
      TCGA_exome    ->      pp;
      TCGA_patient  ->      pp;
      pp_sub        ->      an_mrna_deg ->   an_mrna_hm;
                            an_mrna_deg ->   an_mrna_gage;
                            Data_kegg   ->   an_mrna_gage;
                            an_mrna_deg ->   an_mrna_string -> an_string_plot;
                            Data_string ->   an_mrna_string;
      pp_sub        ->      an_surv     ->   an_surv_plot;
      pp_sub        ->      pp_exome    ->   an_exome;
   }

.. _The Cancer Genome Atlas (TCGA): http://cancergenome.nih.gov/
.. _GAGE method: http://doi.org/10.1186/1471-2105-10-161
.. _STRING database: http://string-db.org/

Packages and Dependencies
-------------------------
There are 24 packages used in this workflow, which depend
on 71 additional packages (dependencies)

**Used packages:**

* *Bioconductor*: limma, edgeR, gage, STRINGdb,

* *CRAN*: ncvreg, boot, lars, lasso2, mda, leaps, data.table, reshape2, ggplot2, magrittr, survival, googleVis, plyr, grid, d3heatmap, ggvis, RColorBrewer, DT, jsonlite

* *Github*: rCharts

**Package dependencies:**

* *Bioconductor*: BiocGenerics, GenomeInfoDb, S4Vectors, Biobase, IRanges, Biostrings, AnnotationDbi, KEGGREST, XVector

* *CRAN*: irlba, Matrix, NMF, lattice, foreach, gridBase, pkgmaker, stringr, colorspace, doParallel, digest, rngtools, cluster, registry, codetools, iterators, xtable, Rcpp, stringi, scales, gtable, MASS, proto, dichromat, labeling, munsell, RSQLite, gsubfn, chron, DBI, graph, png, httr, zlibbioc, R6, jsonlite, curl, mime, RCurl, igraph, hash, gplots, plotrix, sqldf, bitops, gtools, KernSmooth, caTools, gdata, class, RJSONIO, htmlwidgets, dendextend, base64enc, yaml, htmltools, whisker, dplyr, assertthat, lazyeval, shiny, httpuv

Data
------
- Exome, Patient, RNAseq data from TCGA consortium and processed by Phil Cheng
- Protein data downlaoded from STRINGdb package for tax_id 9606 (Human)

License
-------
| Copyright (c) 2015 by Phil Cheng
| Copyright (c) 2015-2016 BeDataDriven B.V.
| License: `GPL version 2 or higher`_

.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html
