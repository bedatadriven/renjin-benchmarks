################################
Simulated GEO matrix
################################

Mutations in DNA that is coding for proteins (genes) or regulatory elements
 are often the cause of most diseases. Affymatirx arrays allow quantification of
 specific sequences of DNA or RNA in
biological samples such as blood, tissue, or tumors.
These arrays give the relative levels of specific transcript (between two
samples) so that these can be compared (between eg. healty vs diseased tissues or
before vs after treatment).

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
test (also know as ‘Mann-Whitney’ test) are performed.

******************************
Packages and Dependencies
******************************
There are 3 packages used in this workflow, which depend
on 7 additional packages from CRAN (dependencies)

+++++++++++++++
Used packages:
+++++++++++++++

- **CRAN**: biclust, s4vd, irlba

++++++++++++++++++++++
Package dependencies:
++++++++++++++++++++++

- **CRAN**: lattice, colorspace, MASS, flexclust, modeltools, biclust, Matrix

+++++++++++++++++++++++
Data
+++++++++++++++++++++++

This workflow uses simulated data from GenBase `data_generator.py <https://github.com/mitdbg/genbase/blob/master/data/data_generator.py>`_ with 500 as size of
columns and rows. Here you can read about the original GenBase study by `Taft R et al., 2014 <http://dx.doi.org/10.1145/2588555.2595633>`_.

********************
License
********************
Copyright (c) 2015 MIT DB Group
based on code from `GenBase <https://github.com/mitdbg/genbase/blob/master/code/R_benchmark/vanilla_R_benchmark.R>`_
Copyright (c) 2015 Hannes Muehleisen
based on code from `GenBase (fork) <https://github.com/hannesmuehleisen/genbase/blob/master/code/R_benchmark/vanilla_R_benchmark.R>`_
Copyright (c) 2015 Ieuan Clay
based on code from `GenBench <https://github.com/biolion/genbench>`_
Copyright (c) 2015-2016 BeDataDriven B.V.
License: `GPL version 2 or higher <http://www.gnu.org/licenses/gpl.html>`_
