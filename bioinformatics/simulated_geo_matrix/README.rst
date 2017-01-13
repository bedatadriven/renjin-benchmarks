
Simulated GEO matrix
====================

Mutations in DNA that is coding for proteins (genes) or regulatory elements are 
often the cause of most diseases. Affymatirx arrays allow quantification of 
specific sequences of DNA or RNA in biological samples such as blood, tissue, or 
tumors. These arrays give the relative levels of specific transcript (between 
two samples) so that these can be compared (between eg. healty vs diseased 
tissues or before vs after treatment).

To do so, DNA sequences unique to specific genes/mutations are printed 
as small spots on a glass surface (arrays), with each spot containing thousends 
of copies of a single sequence. The DNA from each sample is than tagged with 
different colors of flourecent molecules (red or green) and hybredized with the array. 
This way, DNA from sample can bind to sequences on the surface that are antisense of 
its own sequence (specific binding).  The array undergoes several washing steps 
to remove the none-specifically bound DNA sequences, which bind loosly. The 
plate is than scanned and the flourecence intensity for each color at each spot 
is recorded. This intensity is an indicator of the ammount of DNA bound and 
probably the levels of that specific DNA present in the sample. The information 
regarding the specific sequence printed on each spot, control spots, and other 
array specific information is stored in a CDF format/file. By comparing the 
intensity of red to green at each spot you can know which sample contained more 
of that specific transcript. 

For this workflow a simulated dataset was used containing patient 
meta-information, gene expression levels (microarray), and gene ontology 
information for 500 patients and 500 genes. This workflow performs Standard 
statistical methods such as linear regression model (lm.fit from stats package), 
covariance ('cov' from 'stats' package), biclustering ('biclust' package), SVM 
('irlba' package), and differential gene expression using two sample Wilcoxon 
test (also know as 'Mann-Whitney' test) are performed. 

.. graphviz::
   :caption: Diagram of Simulated GEO Matrix analysis workflow.

   digraph SIMULATED_GEO_MATRIX {
    fontname="sans-serif";
    penwidth="0.1";
    edge [comment="Wildcard edge",
          fontname="sans-serif",
          fontsize=10,
          colorscheme="blues3",
          color=2,
          fontcolor=3];
    node [fontname="serif",
          fontsize=13,
          fillcolor="1",
          colorscheme="blues4",
          color="2",
          fontcolor="4",
          style="filled"];
    subgraph cluster0 {
        label="Load data";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c0_geo [shape="invhouse",
                label=<GEO Expression<BR/><I>Simulated, 500x500</I>>];
        c0_gen [shape="invhouse",
                label=<Gene meta info<BR/><I>Simulated, 500x500</I>>];
        c0_pat [shape="invhouse",
                label=<Patient meta info<BR/><I>Simulated, 500x500</I>>];
        c0_gos [shape="invhouse",
                label=<Gene Ontology info<BR/><I>Simulated, 500x500</I>>];
        c0_rds [label=<Read and store<BR/><I>readRDS(), list()</I>>];
        c0_dat [label=<Data>];
        c0_geo -> c0_rds;
        c0_gen -> c0_rds;
        c0_pat -> c0_rds;
        c0_gos -> c0_rds;
        c0_rds -> c0_dat;
    }

    subgraph cluster1 {
        label="Perform linear regression analysis";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c1_dat [label="Data"];
        c0_dat -> c1_dat  [ltail=cluster0,
                           lhead=cluster1];
        c1_gen [label=<Subset Genes<BR/>based on <I>function</I>>];
        c1_res [label=<Select <I>drug.response</I>>];
        c1_aaa [label=<Create expression matrix<BR/>for each gene in each <BR/>patient using <I>merge()</I>>];
        c1_cst [label=<cast to matrix<BR/><I>patientid x geneid</I>>];
        c1_lmf [label=<fit linear regression<BR/>using <I>lm.fit()</I>>];
        c1_prt [label=<<I>print()</I>>];
        c1_dat -> c1_gen  [label="$gen"];
        c1_dat -> c1_res  [label="$pat"];
        c1_dat -> c1_aaa  [label="$geo"];
        c1_gen -> c1_aaa;
        c1_aaa -> c1_cst;
        c1_cst -> c1_lmf  [label="x"];
        c1_res -> c1_lmf  [label="y"];
        c1_lmf -> c1_prt;
        {
            rank=same;
            edge [comment="Wildcard node added automatic in EG."];
            node [comment="Wildcard node added automatic in EG."];
            c1_cst;
            c1_res;
        }

    }

    subgraph cluster2 {
        label="Calculate covariance";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c2_dat [label="Data"];
        c0_dat -> c2_dat  [ltail=cluster0,
                           lhead=cluster2];
        c2_pat [label=<Subset patients<BR/>based on <I>disease</I>>];
        c2_aaa [label=<Create expression matrix<BR/>for each gene in each <BR/>patient using <I>merge()</I>>];
        c2_cst [label=<cast to matrix<BR/><I>patientid x geneid</I>>];
        c2_cov [label=<compute covariance<BR/><I>stats::cov()</I>>];
        c2_scv [label=<select top 25 &#x25; genes<BR/><I>which( &#x3E; .75 &#x2A; max())</I>>];
        c2_prt [label=<Print complete cases<BR/><I>complete.cases(), print()</I>>];
        c2_dat -> c2_pat  [label="$pat"];
        c2_dat -> c2_aaa  [label="$geo"];
        c2_pat -> c2_aaa;
        c2_aaa -> c2_cst;
        c2_cst -> c2_cov;
        c2_cov -> c2_scv;
        c2_scv -> c2_prt;
    }

    subgraph cluster3 {
        label="Perform bi-clustering";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c3_dat [label="Data"];
        c3_pat [label=<Subset paients<BR/>based on <I>gender &#x26; age</I>>];
        c3_aaa [label=<Create expression matrix<BR/>for each gene in each <BR/>patient using <I>merge()</I>>];
        c3_cst [label=<cast to matrix<BR/><I>patientid x geneid</I>>];
        c3_bic [label=<Perform biclustering<BR/><I>biclust(method=BCssvd, K=5)</I>>];
        c3_cdf [label=<Convert to data.frame<BR/><I>biclust::writeclust()</I>>];
        c3_prt [label=<Print complete cases<BR/><I>complete.cases(), print()</I>>];
        c0_dat -> c3_dat  [ltail=cluster0,
                           lhead=cluster3];
        c3_dat -> c3_pat  [label="$pat"];
        c3_dat -> c3_aaa  [label="$geo"];
        c3_pat -> c3_aaa;
        c3_aaa -> c3_cst;
        c3_cst -> c3_bic;
        c3_bic -> c3_cdf;
        c3_cdf -> c3_prt;
    }

    subgraph cluster4 {
        label="Compute largest singular values";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c4_dat [label="Data"];
        c0_dat -> c4_dat  [ltail=cluster0,
                           lhead=cluster4];
        c4_gen [label=<Subset Genes<BR/>based on <I>function</I>>];
        c4_aaa [label=<Create expression matrix<BR/>for each gene in each <BR/>patient using <I>merge()</I>>];
        c4_cst [label=<cast to matrix<BR/><I>patientid x geneid</I>>];
        c4_irl [label=<Perform SVD<BR/><I>irlba(bu=50, nv=50)</I>>];
        c4_prt [label=<Print complete cases<BR/><I>complete.cases(), print()</I>>];
        c4_dat -> c4_gen  [label="$gen"];
        c4_dat -> c4_aaa  [label="$geo"];
        c4_gen -> c4_aaa;
        c4_aaa -> c4_cst;
        c4_cst -> c4_irl;
        c4_irl -> c4_prt;
    }

    subgraph cluster5 {
        label="Statistical test";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c5_dat [label="Data"];
        c0_dat -> c5_dat  [ltail=cluster0,
                           lhead=cluster5];
        c5_cst_ge [label=<cast to matrix<BR/><I>Gene ID x Patient ID</I>>];
        c5_cst_go [label=<cast to matrix<BR/><I>Gene ID x GO ID</I>>];
        c5_st1 [label=<Expression levels<BR/>genes in GO term>];
        c5_st2 [label=<Expression levels<BR/>genes NOT in GO term>];
        c5_wlc [label=<Perform Mann–Whitney test<BR/><I>lapply(), wilcox.test()</I>>];
        c5_sbs [label=<Subset significant GOs<BR/><I>subset(p &#x3C; 1e-3)</I>>];
        c5_prt [label=<Print complete cases<BR/><I>complete.cases(), print()</I>>];
        c5_dat -> c5_cst_ge  [label="$geo"];
        c5_cst_ge -> c5_st1;
        c5_cst_ge -> c5_st2;
        c5_dat -> c5_cst_go  [label="$go"];
        c5_cst_go -> c5_st1  [label="1"];
        c5_cst_go -> c5_st2  [label="0"];
        c5_st1 -> c5_wlc  [label="x"];
        c5_st2 -> c5_wlc  [label="y"];
        c5_wlc -> c5_sbs;
        c5_sbs -> c5_prt;
    }
   }

Packages and Dependencies
-------------------------

There are 3 packages used in this workflow, which depend
on 7 additional packages from CRAN (dependencies)

**Used packages:**

* *CRAN*: biclust, s4vd, irlba

**Package dependencies:**

* *CRAN*: lattice, colorspace, MASS, flexclust, modeltools, biclust, Matrix

Data
------

This workflow uses simulated data from GenBase `data_generator.py`_ with 500 as 
size of columns and rows. Here you can read about the original GenBase study by 
`Taft R et al., 2014`_.

.. _data_generator.py: https://github.com/mitdbg/genbase/blob/master/data/data_generator.py
.. _Taft R et al., 2014: http://dx.doi.org/10.1145/2588555.2595633

License
-------

- Copyright (c) 2015 MIT DB Group based on code from `GenBase <https://github.com/mitdbg/genbase/blob/master/code/R_benchmark/vanilla_R_benchmark.R>`_
- Copyright (c) 2015 Hannes Mühleisen based on code from `GenBase (fork) <https://github.com/hannesmuehleisen/genbase/blob/master/code/R_benchmark/vanilla_R_benchmark.R>`_
- Copyright (c) 2015 Ieuan Clay based on code from `genbench <https://github.com/biolion/genbench>`_
- Copyright (c) 2015-2016 BeDataDriven B.V.  License: `GPL version 2 or higher`_

.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html

