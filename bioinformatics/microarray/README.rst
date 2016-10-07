
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

