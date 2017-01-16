
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
plots the important Quality Control plots. Dataset used in this workflow is from a study by `Ramsey JE et al 2013`_ in which the effects of presence/absence of ZXDC gene 
before and during differentiation of a white blood cell type is studied. The corresponding data can be 
downloaded from Gene Expression Omnibus repository using accession number GSE45417. For more information 
about this workflow please visit `ArrayAnalysis.org`_.

.. figure:: ../../docs/_static/affy.pdf
   :scale: 120 %
   :alt: Diagram affy
   
   Diagram of 'affy' analysis workflow.


Packages and Dependencies
-------------------------

There are 12 packages (mainly affymatrix array analysis related) used in this workflow, which depend
on 28 additional packages from CRAN and Bioconductor (dependencies)

**Used packages:**

* *Bioconductor*: ArrayTools, affy, affycomp, affyPLM, affypdnn, bioDist, simpleaffy, affyQCReport, plier, yaqcaffy

* *CRAN*: gdata, gplots

**Package dependencies:**

* *Bioconductor*: limma, Biobase, BiocInstaller, BiocGenerics, zlibbioc, preprocessCore, affyio, gcrma, Biostrings, XVector, S4Vectors, IRanges, genefilter, AnnotationDbi, annotate, GenomeInfoDb, affyPLM

* *CRAN*: xtable, KernSmooth, DBI, XML, lattice, RSQLite, RColorBrewer, caTools, bitops, gtools, survival



.. _ArrayAnalysis.org: http://www.arrayanalysis.org
.. _Ramsey JE et al 2013: http://dx.doi.org/10.1016/j.molimm.2013.07.001


License
-------

* Copyright (c) 2015 Arrayanalysis.org based on code from http://www.arrayanalysis.org/affyQC/doc_affyQC_R.php
* Copyright (c) 2015-2016 BeDataDriven B.V.  License: `Apache License version 2.0 or higher`_

.. _Apache License version 2.0 or higher: http://www.apache.org/licenses/LICENSE-2.0


.. raw:: latex

    \clearpage

