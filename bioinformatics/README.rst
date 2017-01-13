
##############
Bioinformatics
##############

************
Introduction
************

Bioinformatic workflows used in different fields of biology for educational and benchmarking purposes. This is and will stay a work in progress and we welcome any user to submit PRs to include their workflows. 

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
| integration_igraph              | Integerate different biological data with iGraph                           |
+---------------------------------+----------------------------------------------------------------------------+
| integration_livercohort         | Integrate data from liver cohort study with another data sourc             |
+---------------------------------+----------------------------------------------------------------------------+
| microarray                      | Analysis of microarray datasets                                            |  
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
| singlecell_seurat_mrkr_discov   | Single cell gene expression analysis for marker discovery (using Seurat)   |
+---------------------------------+----------------------------------------------------------------------------+
| singlecell_seurat_spatial_infer | Single cell gene expression analysis for sepatial inference (using Seurat) |
+---------------------------------+----------------------------------------------------------------------------+
| survival_simple                 | Simple survival analysis                                                   |
+---------------------------------+----------------------------------------------------------------------------+
| survival_tcga                   | Survival analysis using TCGA data                                          |
+---------------------------------+----------------------------------------------------------------------------+
| tcga_browser                    | Modules used in a shiny app to visualize TCGA data                         |
+---------------------------------+----------------------------------------------------------------------------+

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
    singlecell_seurat_mrkr_discov/README
    singlecell_seurat_spatial_infer/README
    survival_simple/README
    survival_tcga/README
    tcga_browser/README

