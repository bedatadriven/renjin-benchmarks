
RNAseq Preprocess
==================

All the functions that take place within a cell are performed through proteins.
These proteins are coded within the DNA (Deoxyribonucleic acid) of the cell.
A gene is a sequence of DNA that encodes for a particular protein. In order to
make the necessary proteins, the transcriptional machinary of a cell makes
special copies of the respective genes that can be translated to protein
sequences. These special copies are called messenger RNAs (Ribonucleic acids).

The amount of mRNA produced by a specific gene is used as a surrogate marker for
quantification of the gene activity. With Current highthroughput sequencing
technologies (such as in this case RNAseq), all (fragmented) RNA molecules
from a biological sample are sequenced. These sequences are then matched to
the annotated genome sequence to identify to which gene each sequenced fragment
belongs. More sequenced fragments mapping to a gene in one samples as compared
to another, means that the specific gene had a higher activity.

This workflow uses 'Rsubread' package to perform mapping of sequencing reads to
the genome and to quantify the number of reads mapped to each gene. This
preprocessing step is mostly done outside R using tools such as `STAR aligner`_, `TopHat`_,
`HISAT2`_, `kallisto`_, or `Salmon`_.


Packages and Dependencies
------------------------------
There are 1 packages used in this workflow, which depend
on ... additional packages (dependencies).

**Used packages:**

- **Bioconductor**:

- **CRAN**:

**Package dependencies:**

* *Bioconductor*:

* *CRAN*:

Data
------

RNAseq fastq-files from `Solaimani Kartalaei P, (2014) <http://www.doi.org/10.1084/jem.20140767>`_
samples from SRA id: SRR1653121 to SRR1653136


.. _STAR aligner: https://github.com/alexdobin/STAR
.. _HISAT2: http://ccb.jhu.edu/software/hisat2/index.shtml
.. _kallisto: https://pachterlab.github.io/kallisto/
.. _Salmon: https://github.com/COMBINE-lab/salmon
.. _TopHat: https://ccb.jhu.edu/software/tophat/index.shtml

********************
License
********************
| Copyright (c) 2016 BeDataDriven B.V.
| License: Artistic 2.0
