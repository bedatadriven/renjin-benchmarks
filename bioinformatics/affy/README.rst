
Affy
====

Mutations in DNA that is coding for proteins (genes) or regulatory elements are often the cause of 
most diseases. Affymatirx arrays allow quantification of specific sequences of DNA or RNA in biological 
samples such as blood, tissue, or tumors. These arrays give the relative levels of specific transcript 
(between two samples) so that these can be compared (between eg. healty vs diseased tissues or before 
vs after treatment).

To do so, DNA sequences unique to specific genes/mutations are printed as small spots on a glass 
surface (arrays), with each spot containing thousends of copies of a single sequence. The DNA from 
each sample is than tagged with different colors of flourecent molecules (red or green) and hybredized 
with the array. This way, DNA from sample can bind to sequences on the surface that are antisense of 
its own sequence (specific binding). The array undergoes several washing steps to remove the 
none-specifically bound DNA sequences, which bind loosly. The plate is than scanned and the flourecence 
intensity for each color at each spot is recorded. This intensity is an indicator of the ammount of DNA 
bound and probably the levels of that specific DNA present in the sample. The information regarding the 
specific sequence printed on each spot, control spots, and other array specific information is stored in 
a CDF format/file. By comparing the intensity of red to green at each spot you can know which sample 
contained more of that specific transcript.

This workflow is provided by `ArrayAnalysis.org`_ and its code has been 
merged into single R script. It performs array normalization and differential expression analysis, and 
plots the important Quality Control plots. Dataset used in this workflow is from a study by `Ramsey JE 
et al 2013`_ in which the effects of presence/absence of ZXDC gene 
before and during differentiation of a white blood cell type is studied. The corresponding data can be 
downloaded from Gene Expression Omnibus repository using accession number GSE45417. For more information 
about this workflow please visit `ArrayAnalysis.org`_.

.. graphviz::
   :caption: Diagram of affy analysis workflow.

   digraph AFFY_diagram {
    ranksep=1;
    fontname="sans-serif";
    compound="true";
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
    subgraph "cluster0" {
        label="Raw data QC graphs";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        "c0_rawData" [shape="invhouse", 
                      label="rawData", 
                      colorscheme="ylorrd3"];
        "c0_plot" [shape="box", 
                   label="plot"];
        "c0_rawData" -> "samplePrepPlot()";
        "samplePrepPlot()" -> "c0_plot";
        "c0_rawData" -> "ratioPlot()";
        "ratioPlot()" -> "c0_plot";
        "c0_rawData" -> "RNAdegPlot()";
        "RNAdegPlot()" -> "c0_plot";
    }

    subgraph "cluster1" {
        label="Spike-in Controls";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        "c0_plot" -> "c1_rawData"  [style="bold", 
                                    ltail="cluster0", 
                                    lhead="cluster1", 
                                    penwidth="5.0"];
        "c1_rawData" [shape="invhouse", 
                      label="rawData", 
                      colorscheme="ylorrd3"];
        "c1_plot" [shape="box", 
                   label="plot"];
        "c1_wt" [label="write.table()"];
        "c1_pmacall" [label="deduceSpecies() \n computePMAtable()"];
        "c1_rawData" -> "hybridPlot()";
        "hybridPlot()" -> "c1_plot";
        "c1_rawData" -> "backgroundPlot()";
        "backgroundPlot()" -> "c1_plot";
        "c1_rawData" -> "percPresPlot()";
        "percPresPlot()" -> "c1_plot";
        "c1_rawData" -> "PNdistrPlot()";
        "PNdistrPlot()" -> "c1_plot";
        "c1_rawData" -> "controlPlots()";
        "controlPlots()" -> "c1_plot";
        "c1_rawData" -> "c1_pmacall";
        "c1_pmacall" -> "c1_plot";
        "c1_pmacall" -> "c1_wt";
    }

    subgraph "cluster2" {
        label="Scale factor";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        "c1_plot" -> "c2_rawData"  [ltail="cluster1", 
                                    penwidth="5.0", 
                                    lhead="cluster2"];
        "c2_rawData" [shape="invhouse", 
                      label="rawData", 
                      colorscheme="ylorrd3"];
        "c2_plot" [shape="box", 
                   label="plot"];
        "c2_rawData.pset" [label="rawData.pset"];
        "c2_rawData" -> "scaleFactPlot()";
        "scaleFactPlot()" -> "c2_plot";
        "c2_rawData" -> "boxplotFun()";
        "boxplotFun()" -> "c2_plot";
        "c2_rawData" -> "densityFun()";
        "densityFun()" -> "c2_plot";
        "c2_rawData" -> "maFun()";
        "maFun()" -> "c2_plot";
        "c2_rawData" -> "plotArrayLayout()";
        "plotArrayLayout()" -> "c2_plot";
        "c2_rawData" -> "PNposPlot()";
        "PNposPlot()" -> "c2_plot";
        "c2_rawData" -> "fitPLM()";
        "fitPLM()" -> "c2_rawData.pset";
        "c2_rawData.pset" -> "spatialImages()";
        "c2_rawData" -> "spatialImages()";
        "spatialImages()" -> "c2_plot";
        "c2_rawData" -> "array.image()";
        "array.image()" -> "c2_plot";
        "c2_rawData" -> "nuseFun()";
        "nuseFun()" -> "c2_plot";
        "c2_rawData" -> "rleFun()";
        "rleFun()" -> "c2_plot";
    }

    subgraph "cluster3" {
        label="Correlation Plot";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        "c2_plot" -> "c3_rawData"  [ltail="cluster2", 
                                    penwidth="5.0", 
                                    lhead="cluster3"];
        "c3_rawData" [shape="invhouse", 
                      label="rawData", 
                      colorscheme="ylorrd3"];
        "c3_plot" [shape="box", 
                   label="plot"];
        "c3_rawData" -> "c3_correlFun()";
        "c3_correlFun()" -> "c3_plot";
        "c3_pcaFun()" [label="pcaFun()"];
        "c3_rawData" -> "c3_pcaFun()";
        "c3_pcaFun()" -> "c3_plot";
        "c3_clusterFun()" [label="clusterFun()"];
        "c3_rawData" -> "c3_clusterFun()";
        "c3_clusterFun()" -> "c3_plot";
    }

    subgraph "cluster4" {
        label="Preprocessing";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        "c3_plot" -> "c4_rawData"  [ltail="cluster3", 
                                    penwidth="5.0", 
                                    lhead="cluster4"];
        "c4_rawData" [shape="invhouse", 
                      label="rawData", 
                      colorscheme="ylorrd3"];
        "c4_plot" [shape="box", 
                   label="plot"];
        "c4_png" [label="png"];
        "c4_rawData" -> "deduceSpecies()";
        "deduceSpecies()" -> "c4_plot";
        "c4_rawData" -> "normalizeData()";
        "normalizeData()" -> "c4_plot";
    }

    subgraph "cluster5" {
        label="Preprocessing Evaluation";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        "c4_plot" -> "c5_rawData"  [ltail="cluster4", 
                                    penwidth="5.0", 
                                    lhead="cluster5"];
        "c5_rawData" [shape="invhouse", 
                      label="rawData", 
                      colorscheme="ylorrd3"];
        "c5_plot" [shape="box", 
                   label="plot"];
        "c5_boxplotFun()" [label="boxplotFun()"];
        "c5_rawData" -> "c5_boxplotFun()";
        "c5_boxplotFun()" -> "c5_plot";
        "c5_densityFun()" [label="densityFun()"];
        "c5_rawData" -> "c5_densityFun()";
        "c5_densityFun()" -> "c5_plot";
        "c5_maFun()" [label="maFun()"];
        "c5_rawData" -> "c5_maFun()";
        "c5_maFun()" -> "c5_plot";
        "c5_rawData" -> "correlFun()";
        "correlFun()" -> "c5_plot";
        "c5_rawData" -> "pcaFun()";
        "pcaFun()" -> "c5_plot";
        "c5_rawData" -> "clusterFun()";
        "clusterFun()" -> "c5_plot";
    }

    subgraph "cluster6" {
        label="Prepare output table";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        "c5_plot" -> "c6_rawData"  [ltail="cluster5", 
                                    penwidth="5.0", 
                                    lhead="cluster6"];
        "c6_rawData" [shape="invhouse", 
                      label="rawData", 
                      colorscheme="ylorrd3"];
        "c6_wt" [label="write.table()"];
        "c6_rawData" -> "createNormDataTable()";
        "createNormDataTable()" -> "c6_wt";
    }

}


Packages and Dependencies
-------------------------

There are 12 packages (mainly affymatrix array analysis related) used in this workflow, which depend
on 28 additional packages from CRAN and Bioconductor (dependencies)

**Used packages:**

* *Bioconductor*: ArrayTools, affy, affycomp, affyPLM, affypdnn, bioDist, simpleaffy, affyQCReport, plier, yaqcaffy

* *CRAN*: gdata, gplots

**Package dependencies:**

* *Bioconductor*: limma, Biobase, BiocInstaller, BiocGenerics, zlibbioc, preprocessCore, affyio, gcrma, Biostrings, XVector, S4Vectors, IRanges, genefilter, AnnotationDbi, annotate, GenomeInfoDb, affyPLM

* *CRAN*: xtable, KernSmooth, DBI, XML, lattice, RSQLite, RColorBrewer, caTools, bitops, gtools, survival,



.. _ArrayAnalysis.org: http://www.arrayanalysis.org
.. _Ramsey JE et al 2013: doi.org/10.1016/j.molimm.2013.07.001


License
-------

* Copyright (c) 2015 Arrayanalysis.org based on code from http://www.arrayanalysis.org/affyQC/doc_affyQC_R.php
* Copyright (c) 2015-2016 BeDataDriven B.V.  License: `Apache License version 2.0 or higher`_

.. _Apache License version 2.0 or higher: http://www.apache.org/licenses/LICENSE-2.0

