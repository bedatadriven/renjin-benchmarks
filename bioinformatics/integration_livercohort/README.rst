
Integration: Liver cohort
=========================

Recent high-throughput technologies have caused an explosion in the number of
cohort studies which monitor and record medical (phenotype, genetic), social
(phenotype), and biological (phenotype, genotype, transcriptional) information
of large number of human participants for a long period of time. This allows to
identify the variables that, for example, correlate with occurrence of diseases.
Researchers hope that the causative events that have led to the disease are
among the highest correlating variables. Given the limited number of validation
experiments that a scientist can perform during his/her career, being able to
narrow down the number of candidate variables without loosing the causative
variable is of utmost importance. From population health perpective, earlier
diagnosis of disease or slightly better prediction of disease progress (higher
accuracy and precision) can have significant effects on population health and
its costs. Machine learning algorithms can integrate information from multiple
sources and are, therefor, most widely used.

This workflow uses data acquired in a cohort study by `Schadt EE et al., 2008 <http://doi.org/10.1371/journal.pbio.0060107>`_
on over 400 human liver samples. First gene expression, mutation and phenotype
data from are cleaned up by removing variables containing missing data and
patients which lack expression, mutation, or phenotype data. Multiple machine
leaning algorithms such as Support Vector Machines, Native Bayesian, and Robust
regression are than used to create predictive models.

Support Vector Machine ('e1071' package) is used to train a model
(classification and regression) using random sample of 2/3 of samples
(training set) and tested with the remaining 1/3 of samples (test set). Model
is trained based on independent variables such as age and liver triglyceride
levels, and dependent variables such as activity of nine liver enzymes.

Using enzyme activity information, heatmaps ('stats' package) are generated to
visualize correlation between enzymes and correlation between patients. For
clustering of patients based on enzyme activity data, heatmap are grouped based
on Principal Component Analysis results (prcomp from 'stats' package). Naive
Bayesian classifier ('e1071' package) is used to cluster samples based on gene
expression profile with aldehyde oxydase levels or liver enzyme activity as
class vector (independent variable).

Furthermore, Robust linear model ('MASS' package) is used to train model using
liver triglyceride levels and gene expression levels of genes with highest
variance. A random sample of 25 genes are selected from 1000 genes with the
highest variance and used in combination with triglyceride level phenotype.
This is done in 50 iteration and the iteration with highest correlation in
training and test sets is recorded.


Packages and Dependencies
-------------------------

There are 3 packages used in this workflow, which depend
on 1 additional package from CRAN (dependency)

**Used packages:**

* *CRAN*: stats, e1071, MASS

**Package dependencies:**

* *CRAN*: class

Data
-------

from `Human Liver Cohort (Synapse ID: syn4499) <https://www.synapse.org/#!Synapse:syn4499>`_

License
-------

* Copyright (c) 2015 Ieuan Clay based on code from `genbench <https://github.com/biolion/genbench>`_
* Copyright (c) 2015-2016 BeDataDriven B.V.  License: `GPL version 2 or higher`_

.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html

