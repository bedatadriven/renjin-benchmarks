
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


.. graphviz::
   :caption: Diagram for integration livercohort benchmark.

   digraph INTEGRATION_LIVERCOHORT {
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
    subgraph cluster0 {
        label="Load & prepare data";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        load [shape="box", 
              label="Load \n read.delim()"];
        process [shape="box", 
                 label="Process: \n subset(), merge(), \n complete.cases()"];
        load -> process;
    }

    subgraph cluster1 {
        label="SVM predicting modeling";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        process -> dataset  [ltail="cluster0", 
                             lhead="cluster1"];
        dataset [shape="invhouse", 
                 label="Devide in train and test \n sets using rep(), sample()"];
        dataset -> "trainset"  [label="1/3"];
        "trainset" -> "svm()";
        "svm()" -> "predict()"  [label="model"];
        "predict()" -> "classAgreement";
        "trainset" -> "predict()";
        dataset -> "testset"  [label="2/3"];
        "testset" -> "predict()";
        {
            rank=same;
            edge [comment="Wildcard node added automatic in EG."];
            node [comment="Wildcard node added automatic in EG."];
            "trainset";
            "testset";
        }

    }

    subgraph cluster2 {
        label="NaiveBayesian modeling";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        process -> c2_dataset  [ltail="cluster0", 
                                lhead="cluster2"];
        c2_explr [shape="box", 
                  label="Explore dataset using: \n hclust(), prcomp(), t(), \n dist(), cutree(), cor()"];
        c2_explt [shape="box", 
                  label="Plot exploratoy analysis \n with heatmap() and pairs()"];
        c2_explr -> c2_explt;
        c2_dataset [shape="invhouse", 
                    label="curatedPhen"];
        c2_make_cat [shape=box, 
                     label="create binary categories \n cut(quantile()), \n cutree(hclust())"];
        c2_train [label="trainset"];
        c2_test [label="testset"];
        c2_nb [label="naiveBayes()"];
        c2_pred [label="predict()"];
        c2_clsagr [label="classAgreement()"];
        c2_dataset -> c2_explr;
        c2_explr -> c2_make_cat;
        c2_dataset -> c2_train  [label="1/3"];
        c2_dataset -> c2_test  [label="2/3"];
        c2_train -> c2_nb;
        c2_nb -> c2_pred  [label="model"];
        c2_test -> c2_pred;
        c2_pred -> c2_clsagr;
        c2_make_cat -> c2_nb;
        c2_train -> c2_pred;
        {
            rank=same;
            edge [comment="Wildcard node added automatic in EG."];
            node [comment="Wildcard node added automatic in EG."];
            c2_train;
            c2_test;
        }

    }

    subgraph cluster3 {
        label="Robust Linear Model fitting (RLM)";
        edge [comment="Wildcard node added automatic in EG."];
        node [comment="Wildcard node added automatic in EG."];
        process -> c3_expre  [ltail="cluster0", 
                              lhead="cluster3"];
        c3_pheno [shape="invhouse", 
                  label="curatedPhen"];
        c3_expre [shape="invhouse", 
                  label="curatedExpr"];
        c3_dataset [shape="invhouse", 
                    label="Devide in train and test \n sets using rep(), sample()"];
        c3_train [label="trainset"];
        c3_test [label="testset"];
        c3_expre -> c3_dataset;
        c3_dataset -> c3_train  [label="1/3"];
        c3_dataset -> c3_test  [label="2/3"];
        c3_feats [label="selected features"];
        c3_col_feat [shape="box", 
                     label="Remove low variance columns \n var(), rank()"];
        c3_row_feat [shape="box", 
                     label="Remove high correlation rows \n sum(), abs(), cor()"];
        c3_rlm_tri [label="rlm(triglyc ~ ., data)"];
        c3_pred [label="predict()"];
        c3_cor [label="cor()"];
        c3_train -> c3_col_feat;
        c3_col_feat -> c3_feats;
        c3_row_feat -> c3_feats;
        c3_feats -> c3_rlm_tri;
        c3_feats -> c3_pred;
        c3_feats -> c3_cor;
        c3_train -> c3_rlm_tri;
        c3_rlm_tri -> c3_pred  [label="model"];
        c3_pred -> c3_cor;
        c3_test -> c3_pred;
        c3_pheno -> c3_rlm_tri;
        c3_pheno -> c3_cor;
        {
            rank=same;
            edge [comment="Wildcard node added automatic in EG."];
            node [comment="Wildcard node added automatic in EG."];
            c3_train;
            c3_test;
        }

        {
            rank=same;
            edge [comment="Wildcard node added automatic in EG."];
            node [comment="Wildcard node added automatic in EG."];
            c3_pheno;
            c3_expre;
        }

    }

}


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

