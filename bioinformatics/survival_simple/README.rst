
Survival simple
===============

All cancers are caused by mutations in DNA. `The Cancer Genome Atlas (TCGA)`_
is a project that begun in 2005 and aims to catalogue the mutations that cause
cancer in patients. It uses highthroughput sequencing technologies and
bioinformatics to achive this gaol. Beside genomic data, TCGA also records and
provides detailed anonymize meta information about each patient (such as
disease stage, patient age, sex, healt, etc) which are very valuable for
identification of risk factors.

In this workflow, survival analysis is performed on TCGA patient survival data
with BAZ2A gene mutations as predictor. The `survival package`_ is one of the
most used packages for survival analysis in R.

.. _The Cancer Genome Atlas (TCGA): http://cancergenome.nih.gov/
.. _survival package: https://cran.r-project.org/web/packages/survival/index.html

Packages and Dependencies
-------------------------

There is only 1 package used in this workflow, which has
no additional dependencies.

**Used packages:**

* *CRAN*: survival

Data:
------
- pat.gene.rda: patient survival data from TCGA consortium processed and formatted by Phil Cheng

License
-------
| Copyright (c) 2015 Phil Cheng
| Copyright (c) 2015-2016 BeDataDriven B.V.
| License: `GPL version 2 or higher`_

.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html