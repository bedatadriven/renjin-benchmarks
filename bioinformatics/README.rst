.. _Bioinformatics:

##############
Bioinformatics
##############

************
Introduction
************

Here we have gathered bioinformatic workflows used in different fields of biology. For each benchmark we have 
included a brief description containing information about its usage, package dependencies, and sources of input 
data. Moreover, each benchmark is stored in a seperate folder with defined structure and includes information 
about the package and dataset versions used.

.. figure:: docs/_static/benchmark_structure.pdf
   :scale: 75 %
   :alt: Benchmark structure
   :figwidth: 50 %

   The folder name is the name of the benchmark and the file containing the R code. Each folder 
   contains a Debian Control File containing Title and Description tags and for each input file a File, URL, and 
   Hash tags. In addition, there is “packrat.R” file containing the detailed information about each package that 
   is used in the workflow. 


So far the benchmarks include analysis of gene expression profiling (array and sequencing based 
methods), protein analysis, gene set enrichment analysis, network analysis, analysis of clinical data, and 
visualization of these data.

List of bioinformatic benchmarks included in this repository.

+---------------------------------+----------------------------------------------------------------------------+
| **Benchmark name**              | **Benchmark description**                                                  |
+---------------------------------+----------------------------------------------------------------------------+
| affy                            | Analysis of MicroArray data                                                |
+---------------------------------+----------------------------------------------------------------------------+
| clinical_eslII                  | Analysis of clinical data (according to book ESL 2nd edition)              |
+---------------------------------+----------------------------------------------------------------------------+
| clinical_livercohort            | Analysis of clinical data from liver cohort study                          |
+---------------------------------+----------------------------------------------------------------------------+
| generate_count                  | Generate count tables from RNAseq data                                     |
+---------------------------------+----------------------------------------------------------------------------+
| integration_igraph              | Integrate different biological data with iGraph                            |
+---------------------------------+----------------------------------------------------------------------------+
| integration_livercohort         | Integrate data from liver cohort study with another data source            |
+---------------------------------+----------------------------------------------------------------------------+
| microarray                      | Analysis of MicroArray datasets                                            |  
+---------------------------------+----------------------------------------------------------------------------+
| mutation                        | Analysis of SNP data                                                       |
+---------------------------------+----------------------------------------------------------------------------+
| rnaseq_deseq2                   | Analysis of RNAseq data using DESeq2                                       |
+---------------------------------+----------------------------------------------------------------------------+
| rnaseq_preprocess               | Preprocessing of RNAseq files to generate count data from FASTA files      |
+---------------------------------+----------------------------------------------------------------------------+
| rppa                            | Analysis of Reverse Phase Protein Array (RPPA) data                        |
+---------------------------------+----------------------------------------------------------------------------+
| simulated_geo_matrix            | Analysis of simulated gene expression and patient meta data                |
+---------------------------------+----------------------------------------------------------------------------+
| survival_simple                 | Simple survival analysis                                                   |
+---------------------------------+----------------------------------------------------------------------------+
| survival_tcga                   | Survival analysis using TCGA data                                          |
+---------------------------------+----------------------------------------------------------------------------+
| tcga_browser                    | Modules used in a shiny app to visualize TCGA data                         |
+---------------------------------+----------------------------------------------------------------------------+

We accept contributions to this list through Pull Requests (PR) after we have validated the requested benchmark for 
correctness and if the workflow is not covered by existing benchmarks.

************
Workflows
************

.. toctree::
    :maxdepth: 1

    affy/README
    clinical_eslII/README
    clinical_livercohort/README
    generate_count/README
    integration_igraph/README
    integration_livercohort/README
    microarray/README
    mutation/README
    rnaseq_deseq2/README
    rppa/README
    simulated_geo_matrix/README
    survival_simple/README
    survival_tcga/README
    tcga_browser/README

.. raw:: latex

    \clearpage

********************
Performance in GNU R
********************

.. toctree::
    :maxdepth: 1

    gnuR_profile.rst

