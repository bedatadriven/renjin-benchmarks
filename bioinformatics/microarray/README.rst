
Microarray
==========

Mutations in DNA that is coding for proteins (genes) or regulatory elements are 
often the cause of most diseases. Affymatirx arrays allow quantification of specific 
sequences of DNA or RNA in biological samples such as blood, tissue, or tumors. 
These arrays give the relative levels of specific transcript (between two samples) 
so that these can be compared (between eg. healty vs diseased tissues or before vs 
after treatment).

To do so, DNA sequences unique to specific genes/mutations are printed as small spots 
on a glass surface (arrays), with each spot containing thousends of copies of a single 
sequence. The DNA from each sample is than tagged with different colors of flourecent 
molecules (red or green) and hybredized with the array. This way, DNA from sample can 
bind to sequences on the surface that are antisense of its own sequence (specific 
binding).  The array undergoes several washing steps to remove the none-specifically 
bound DNA sequences, which bind loosly. The plate is than scanned and the flourecence 
intensity for each color at each spot is recorded. This intensity is an indicator of 
the ammount of DNA bound and probably the levels of that specific DNA present in the 
sample. The information regarding the specific sequence printed on each spot, control 
spots, and other array specific information is stored in a CDF format/file. By 
comparing the intensity of red to green at each spot you can know which sample 
contained more of that specific transcript.

This workflow uses `limma`_ package for analysis of microarray data to analyse the 
data from `Ramsey et al 2013`_. The raw data is normalized and differential gene 
expression analysis and a simple gene set testing is performed as described in limma 
vignet.


.. graphviz::
   :caption: Diagram for the microarray benchmark.

   digraph MICROARRAY {
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
          fontsize=13, 
          fillcolor="1", 
          colorscheme="blues4", 
          color="2", 
          fontcolor="4", 
          style="filled"];
    subgraph "cluster0" {
        label="(down)load data";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c0_data [shape="invhouse", 
                 label="GSE45417_RAW.tar"];
        c0_anno_df [shape="box", 
                    label="create an AnnotatedDataFrame() \n data.frame(Name, FileName, Group, Treatment, Replicate)"];
        c0_read [label="read.affybatch()"];
        c0_data -> c0_anno_df;
        c0_anno_df -> c0_read;
    }

    subgraph "cluster1" {
        label="Perform QC";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c1_qc [label="affy::AffyRNAdeq()"];
        c1_drop [shape="box", 
                 label="Drop samples with abnormal \n RNA degradattion slope"];
        c0_read -> c1_qc  [ltail="cluster0", 
                           lhead="cluster1"];
        c1_filt [label="filtered samples"];
        c1_qc -> c1_drop;
        c1_drop -> c1_filt;
    }

    subgraph "cluster2" {
        label="Normalize and scale";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c2_filt [label="filtered samples"];
        c1_filt -> c2_filt  [ltail="cluster1", 
                             lhead="cluster2"];
        c2_rma [label="Robust multi-array average \n rma(normalize=TRUE)"];
        c2_mas [label="MAS5 normalization \n mas5(normalize=TRUE)"];
        c2_mas_log [label="log2(exprs())"];
        c2_filter [label="filter out probes \n< 50 in 25% samples"];
        c2_filt -> c2_rma;
        c2_filt -> c2_mas;
        c2_mas -> c2_mas_log;
        c2_mas_log -> c2_filter;
        c2_mas_filt [label="filtered data"];
        c2_filter -> c2_mas_filt;
        {
            rank=same;
            edge [comment="Wildcard node added automatic in EG."];
            node [comment="Wildcard node added automatic in EG."];
            c2_rma;
            c2_mas;
        }

    }

    subgraph "cluster3" {
        label="Limma analysis";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c3_mas_filt [label="filtered data"];
        c2_mas_filt -> c3_mas_filt  [ltail="cluster2", 
                                     lhead="cluster3"];
        c3_grp [label="pData()$Group"];
        c3_trt [label="pData()$Treatment"];
        c3_design [label="model.matirx(~group*treatment)"];
        c3_fit [label="lmFit()"];
        c3_cntm [label="cbind(0010,0011,0001)"];
        c3_fit2 [label="contrast.fit()"];
        c3_fit3 [label="eBayes()"];
        c3_top [label="topTable(adjust=BH)"];
        c3_mas_filt -> c3_grp;
        c3_grp -> c3_design;
        c3_mas_filt -> c3_trt;
        c3_trt -> c3_design;
        c3_design -> c3_fit;
        c3_mas_filt -> c3_fit;
        c3_fit -> c3_fit2;
        c3_fit2 -> c3_fit3;
        c3_fit3 -> c3_top;
        c3_cntm -> c3_fit2;
        {
            rank=same;
            edge [comment="Wildcard node added automatic in EG."];
            node [comment="Wildcard node added automatic in EG."];
            c3_grp;
            c3_trt;
        }

    }

    subgraph "cluster4" {
        label="Gene Set Analysis";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        c4_deg [label="simulate"];
        c3_top -> c4_deg  [ltail="cluster3", 
                           lhead="cluster4"];
        c4_simd [label="data \n 1e4 x 20"];
        c4_gset [label="gene sets"];
        c4_mod [label="model.matrix"];
        c4_rmr [label="limma::romer()"];
        c4_rmrtop [label="topRomer()"];
        c4_rst [label="limma::mroast()"];
        c4_rsttop [label="Get top hits \n order(), subset()"];
        c4_rbnd [label="rbind()"];
        c4_deg -> c4_mod;
        c4_mod -> c4_rmr;
        c4_mod -> c4_rst;
        c4_deg -> c4_simd;
        c4_simd -> c4_rmr;
        c4_rmr -> c4_rmrtop;
        c4_rmrtop -> c4_rbnd;
        c4_deg -> c4_gset;
        c4_gset -> c4_rmr;
        c4_simd -> c4_rst;
        c4_rst -> c4_rsttop;
        c4_rsttop -> c4_rbnd;
        c4_gset -> c4_rst;
        {
            rank=same;
            edge [comment="Wildcard node added automatic in EG."];
            node [comment="Wildcard node added automatic in EG."];
        }

    }

}


Packages and Dependencies
-------------------------

There are 4 packages used in this workflow, which depend
on 5 additional packages (dependencies).

**Used packages:**

* *Bioconductor*: Biobase, affy, hgu133plus2cdf, limma

**Package dependencies:**

* *Bioconductor*: BiocGenerics, BiocInstaller, zlibbioc, preprocessCore, affyio

Data
--------

From Gene Expression Omnibus repository accession ID `GSE45417`_.

.. _limma: http://www.bioconductor.org/packages/release/bioc/html/limma.html
.. _Ramsey et al 2013: http://doi.org/10.1016/j.molimm.2013.07.001
.. _GSE45417: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45417

License
-------

* Copyright (c) 2005 Gordon Smyth based on `Limma package <http://www.bioconductor.org/packages/release/bioc/html/limma.html>`_
* Copyright (c) 2015 Ieuan Clay based on code from `genbench <https://github.com/biolion/genbench>`_
* Copyright (c) 2015-2016 BeDataDriven B.V.  License: `GPL version 2 or higher`_

.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html

