
Mutation
========

Two main types of genetic studies are population and familial/pedigree studies.
In this workflow, individual mutation information is used to determine the
relatedness between individuals and data from The Cancer Genome Atlas landmark
paper on most common `AML mutations`_ is used to reproduce some of the
figures in this publication.

.. graphviz::
   :caption: Diagram of mutation analysis workflow.

   digraph MUTATION {
       fontname="sans-serif";
       compound="true";
       penwidth="0.1";
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
           label="Read data";
           edge [comment="Wildcard node added automatic in EG."];
           node [comment="Wildcard node added automatic in EG."];
           c0_input1 [shape="invhouse",
                      label="AML genotype"];
           c0_input2 [shape="invhouse",
                      label="AML meta"];
           c0_read [label=<Read data<BR/><I>read.delim()</I>>];
           c0_check [shape="box",
                     label=<Validate data<BR/><I>length(), intersect(), order()</I>>];
           c0_out [label="list()"];
           c0_input1 -> c0_read;
           c0_read -> c0_out;
           c0_input2 -> c0_read;
           c0_out -> c0_check;
           {
               rank=same;
               edge [comment="Wildcard node added automatic in EG."];
               node [comment="Wildcard node added automatic in EG."];
               c0_input1;
               c0_input2;
           }

           {
               rank=same;
               edge [comment="Wildcard node added automatic in EG."];
               node [comment="Wildcard node added automatic in EG."];
               c0_out;
               c0_check;
           }

       }

       subgraph cluster1 {
           label="reproduce Figure 1a";
           edge [comment="Wildcard node added automatic in EG."];
           node [comment="Wildcard node added automatic in EG."];
           c1_in [shape="invhouse",
                  label=<aggragate and format data<BR/><I>rbind(), split(), table('tier')</I>>];
           c1_meta [shape="invhouse",
                    label="Meta data \n input$meta"];
           c1_mrg [label=<Merge with meta information<BR/><I>merge(x, y)</I>>];
           c1_subs [label="subset('tier1')"];
           c1_plot [label="plot()"];
           c0_out -> c1_in  [ltail=cluster0,
                             lhead=cluster1];
           c0_out -> c1_meta  [ltail=cluster0,
                               lhead=cluster1];
           c1_in -> c1_mrg  [label="x"];
           c1_meta -> c1_mrg  [label="y"];
           c1_mrg -> c1_subs;
           c1_subs -> c1_plot;
       }

       subgraph cluster2 {
           label="reproduce Figure 1b";
           edge [comment="Wildcard node added automatic in EG."];
           node [comment="Wildcard node added automatic in EG."];
           c2_in [shape="invhouse",
                  label=<aggragate and format data<BR/><I>rbind(), split(), table('gene_name')</I>>];
           c0_out -> c2_in  [ltail=cluster0,
                             lhead=cluster2];
           c2_ordr [label="order(sample_count)"];
           c2_subs [label="subset('tier1')"];
           c2_head [label="head(100)\ndata.frame(gene_name, tier, sample_count)"];
           c2_in -> c2_ordr;
           c2_ordr -> c2_subs;
           c2_subs -> c2_head;
       }

       subgraph cluster3 {
           label="Load family data";
           edge [comment="Wildcard node added automatic in EG."];
           node [comment="Wildcard node added automatic in EG."];
           c3_input [label="Read genotyping data\n read.delim(genotypes)"];
           c3_genot [label=<Summarize genotype freq's <BR/><I>rbind(), table('genotype')</I>>];
           c3_compu [shape="box",
                     label=<create mutation likelihood matrix <BR/><I>matrix(ncol(genotypes), nrow(genotypes))</I>>];
           c3_match [shape="box",
                     label=<compute likelihood scores<BR/><I>number of matching alleles * allele occurence</I>>];
           c3_ibdve [shape="box",
                     label=<Run sliding window over mutations and<BR/>compute 3rd Quantile of mean marker likelihoods<BR/>for different window sizes using<BR/><I>lapply(), summary(), base:::simplify2array(), seq()</I>>];
           c3_input -> c3_genot;
           c3_genot -> c3_compu;
           c3_compu -> c3_match;
           c3_match -> c3_ibdve;
       }

   }


.. _AML mutations: http://www.doi.org/10.1056/NEJMoa1301689

Packages and Dependencies
-------------------------

There are 2 core packages used in this workflow, which have no dependencies.

**Used packages:**

* *Core*: stats, utils

Data
------

The familial data was obtained from the nice people at `Genomes Unzipped`_, who
make `their own genomic data`_ publicly available.

The dataset we are using comes from the `23andme v2`_ sequencing service.

Though the individuals are not related, this data can still be used to perform
some typical tests carried out on pedigree studies, such as determining
"relatedness" between individuals.

+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| member            | dataset id | link                                                                                                                       |
+===================+============+============================================================================================================================+
| Daniel MacArthur  | DGM001     | `http://s3.amazonaws.com/gnz.genotypes/DGM001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/DGM001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Luke Jostins      | LXJ001     | `http://s3.amazonaws.com/gnz.genotypes/LXJ001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/LXJ001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Dan Vorhaus       | DBV001     | `http://s3.amazonaws.com/gnz.genotypes/DBV001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/DBV001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Caroline Wright   | CFW001     | `http://s3.amazonaws.com/gnz.genotypes/CFW001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/CFW001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Kate Morley       | KIM001     | `http://s3.amazonaws.com/gnz.genotypes/KIM001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/KIM001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Vincent Plagnol   | VXP001     | `http://s3.amazonaws.com/gnz.genotypes/VXP001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/VXP001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Jeff Barrett      | JCB001     | `http://s3.amazonaws.com/gnz.genotypes/JCB001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/JCB001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Jan Aerts         | JXA001     | `http://s3.amazonaws.com/gnz.genotypes/JXA001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/JXA001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Joe Pickrell      | JKP001     | `http://s3.amazonaws.com/gnz.genotypes/JKP001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/JKP001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Don Conrad        | DFC001     | `http://s3.amazonaws.com/gnz.genotypes/JKP001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/JKP001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Carl Anderson     | CAA001     | `http://s3.amazonaws.com/gnz.genotypes/CAA001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/CAA001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+
| Ilana Fisher      | IPF001     | `http://s3.amazonaws.com/gnz.genotypes/IPF001_genotypes.zip <http://s3.amazonaws.com/gnz.genotypes/IPF001_genotypes.zip>`_ |
+-------------------+------------+----------------------------------------------------------------------------------------------------------------------------+

The population study data is from the TCGA consortium publication `TCGA,
2013`_, publication `data archive`_, mutation and annotation (`maf`_), and
`patient meta data`_.

.. _Genomes Unzipped: http://genomesunzipped.org/members
.. _their own genomic data: http://genomesunzipped.org/data
.. _23andme v2: https://www.23andme.com/
.. _TCGA, 2013: http://www.doi.org/10.1056/NEJMoa1301689
.. _data archive: https://tcga-data.nci.nih.gov/docs/publications/laml_2012/
.. _maf: http://tcga-data.nci.nih.gov/docs/publications/laml_2012/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.12.0.tar.gz
.. _patient meta data: http://tcga-data.nci.nih.gov/docs/publications/laml_2012/clinical_patient_laml.tsv

License
-------

* Copyright (c) 2015 Ieuan Clay based on code from `genbench <https://github.com/biolion/genbench>`_
* Copyright (c) 2015-2016 BeDataDriven B.V.  License: `GPL version 2 or higher`_

.. _GPL version 2 or higher: http://www.gnu.org/licenses/gpl.html
