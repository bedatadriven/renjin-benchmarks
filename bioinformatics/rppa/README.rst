###################################
Reverse phase protein array (rppa)
###################################

Processes within a cell are mainly performed by proteins. Special copies of a
gene, called messenger RNAs (mRNAs), are produced and transported to the cells
protein synthesis machinery to produce the corresponding protein. The more mRNA
molecules of a gene cause the more of its corresponding protein is produced.
However, this is not always a linear relationship. Presence or absence of some
molecules (including other proteins) can inhibit or enhance the production of a
protein.

High-throughput technologies such as micro-arrays and RNA/DNA sequencing
technologies measure the level of mRNA as a surrogate marker for protein level.
There are a few highthroughput technologies that allow direct measurement of the
protein levels. The reverse phase protein array (RPPA) is one such technology.
To perform RPPA, cell/tissue lysate of interest is printed as small dots
(droplets) of the same size on surface of a special plate. Each dot on the plate
is then stained with fluorochrome conjugated antibodies against one specific
protein. The plate is then washed and scanned to detect the fluorescence. Since
most antibodies bind to multiple proteins due to non-specific binding, the
requirement for antibodies which only detect a single protein in one of the
limitation of RPPA.

In this workflow, `RPPA dataset<http://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/RPPA_input.csv>`_ from TCGA consortium is
clustered using R hierarchical clustering  and K-means clustering from ‘stats’
package.


******************************
Packages and Dependencies
******************************
There is ony 1 package used in this workflow, which has no dependencies.

+++++++++++++++
Used packages:
+++++++++++++++

- **CORE**: stats

+++++++++++++++++++++++
Data
+++++++++++++++++++++++

- RPPA data: http://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/RPPA_input.csv

********************
License
********************
Copyright (c) 2015 Ieuan Clay
based on code from `GenBench <https://github.com/biolion/genbench>`_
Copyright (c) 2015-2016 BeDataDriven B.V.
License: `GPL version 2 or higher <http://www.gnu.org/licenses/gpl.html>`_
