################################
Integration: iGraph
################################

Genes are sequences of DNA that encode functional molecules (such
as RNAs or Proteins). These functional molecules most of time interact
and affect eachothers function. Therefor, to understand the role of a
specific gene in a disease, scientist need to know about all the published
functions.

Given many published scientific manuscripts about any given gene, the
`National Center for Biotechnology Information (NCBI) <.http://www.ncbi.nlm.nih.gov/>`_ has introduced
`GeneRIF <http://www.ncbi.nlm.nih.gov/gene/about-generif>`_. GeneRIF provides
scientist a simple mechanism to add functional annotation to known genes based
on scientific publications. These annotations are than processed and approved
by GeneRIF staff. All these functional annotations are freely accessible and
can be used by anyone to extract functions of a gene or genes with specific
function.

This workflow imports human Phophatases and Kinases and the publications. And
uses 'igraph' to store this information and to plot the relation between these
publications as well as interaction network between molecules.


******************************
Packages and Dependencies
******************************
There are 7 packages used in this workflow, which depend
on 35 additional packages from CRAN (dependencies)

+++++++++++++++
Used packages:
+++++++++++++++

- **CRAN**: igraph, XML, R.utils, plyr, reshape, utils, sqldf

++++++++++++++++++++++
Package dependencies:
++++++++++++++++++++++

- **CRAN**: magrittr, irlba, Matrix, NMF, lattice, foreach,
gridBase, pkgmaker, reshape2, stringr, colorspace,
doParallel, digest, rngtools, ggplot2, RColorBrewer,
cluster, registry, codetools, iterators, xtable, plyr,
Rcpp, stringi, scales, gtable, MASS, proto, dichromat,
labeling, munsell, RSQLite, gsubfn, chron, DBI

+++++++++++++++++++++++
Data
+++++++++++++++++++++++
Basic GeneRIF data:
  - ftp://ftp.ncbi.nih.gov/gene/GeneRIF/generifs_basic.gz
SQL script written by Ieuan Clay:
  - get_all_RIF.sql
Human kinases and phosphatases protein:
  - GeneOntology.org list of human (Taxonomy id: 9606) genes using terms GO:0050222
(protein kinase activity) and GO:0004721 (phosphoprotein phosphatase activity).


********************
License
********************
Copyright (c) 2015 Ieuan Clay
based on code from https://github.com/biolion/genbench
Copyright (c) 2015-2016 BeDataDriven B.V.
License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
