
Generate count (From RNAseq BAM files)
======================================

The amount of RNA produced by a specific gene is used as a surrogate marker for detection the activity,
of that gene. With Current highthroughput sequencing technologies (such as in this case RNAseq), all (fragmented) RNA molecules
from a biological sample are sequenced. These sequences are then matched to the annotated genome sequence
to identify to which gene each sequenced fragment belongs. More sequenced fragments mapping to a gene in
one samples as compared to another, means that the specific gene had a higher activity.

The first step in this workflow is reading the relevant parts of genome
annotation file and converting it the a list of ranges of feature of interest
such as exon (IRanges, GenomicRanges). Sequencing products previously aligned
to genome are imported using Rsamtools. 
Mapping of mapped
 sequences to genes as defined in the anotation file is the next step. Followed by aggragation of data by
 summing the number of mapped fragments per gene.


Packages and Dependencies
-------------------------

There are 5 packages used in this workflow, which depend
on 6 additional packages from Bioconductor and CRAN (dependencies)

Used packages:
^^^^^^^^^^^^^^

- **Bioconductor**: BiocGenerics, GenomicRanges, IRanges, Rsamtools

- **CRAN**: parallel

Package dependencies:
^^^^^^^^^^^^^^^^^^^^^

- **Bioconductor**: Biostrings, GenomeInfoDb, XVector, S4Vectors, zlibbioc

- **CRAN**: bitops

Data:
^^^^^

Bam files are single-end RNAseq data from Kartalaei PS et al, 2015.

License
-------

Copyright (c) 2015-2016 BeDataDriven B.V.
License: http://www.gnu.org/licenses/gpl.html
