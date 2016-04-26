################################
clinical (liver cohort)
################################

Single nucleotide ploymorphisms (SNPs) are single basepair mutations in DNA that
 occur with certain frequency within the population. These mutations might be the
 cause of a disease or due to their close vecinity to genetic abration they might
 be co-transfered with the cause.

In the study by `Schadt EE et al., 2008 <http://dx.doi.org/10.1371/journal.pbio.0060107>`_
more than 780 thousand SNPs were
analyzed in more than 400 human liver samples. Since a mutation could cause a
phenotype that in time can cause a disease, the phenotypic information about the
patients was included as well.

In this workflow, these data will be used as input for machine learning
classification using SVM (from e1071 package). Data is devided into training and
test sets and predictive models are generated for different features, these
models are then tested using test dataset and the results are stored in a list.


******************************
Packages and Dependencies
******************************
There is 1 package used in this workflow, which depends
on 2 additional packages from CRAN (dependencies)

+++++++++++++++
Used packages:
+++++++++++++++

- **CRAN**: e1071

++++++++++++++++++++++
Package dependencies:
++++++++++++++++++++++

- **CRAN**: class, MASS

+++++++++++++++++++++++
Data
+++++++++++++++++++++++
Downloaded from https://www.synapse.org/ in June 2015

********************
License
********************
Copyright (c) 2015 Ieuan Clay
based on code from https://github.com/biolion/genbench
Copyright (c) 2015-2016 BeDataDriven B.V.
License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
