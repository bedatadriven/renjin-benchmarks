
Survival TCGA
=============

All cancers are caused by mutations in DNA. `The Cancer Genome Atlas (TCGA)`_
is a project that begun in 2005 and aims to catalogue the mutations that cause
cancer in patients. It uses highthroughput sequencing technologies and
bioinformatics to achive this gaol. Beside genomic data, TCGA also records and
provides detailed anonymize meta information about each patient (such as
disease stage, patient age, sex, healt, etc) which are very valuable for
identification of risk factors.

In this workflow, survival analysis is performed on TCGA meta and survival data
using the Cox regression model. The `glmnet package`_ is one of the most
efficient packages for such an analysis in R.

.. _The Cancer Genome Atlas (TCGA): http://cancergenome.nih.gov/
.. _glmnet package: https://cran.r-project.org/web/packages/glmnet/index.html

Packages and Dependencies
-------------------------

There are 2 packages used in this workflow, which depend
on 5 additional packages from CRAN (dependencies)

**Used packages:**

* *CRAN*: glmnet, survival

**Package dependencies:**

* *CRAN*: foreach, Matrix, codetools, iterators, lattice

Data:
------
- Patient survival data from TCGA consortium provided by Andre Verissimo.

License
-------
| Copyright (c) 2015 Andre Verissimo (andre.verissimo@tecnico.ulisboa.pt)
| Copyright (c) 2015-2016 BeDataDriven B.V.
| License: `GPL version 2 or higher`_

.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html