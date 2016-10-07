
Mutation
========

Two main types of genetic studies are population and familial/pedigree studies.
In this workflow, individual mutation information is used to determine the
relatedness between individuals and data from The Cancer Genome Atlas landmark
paper on most common `AML mutations`_ is used to reproduce some of the
figures in this publication.

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

