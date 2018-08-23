
############
Introduction
############

************************
Benchmarks for Renjin
************************


The purpose of this project is to collect *benchmarks* for `Renjin`_, the R
interpreter in the Java Virtual Machine. This document includes a general
introduction to each of these benchmarks. For the purpose of this project, a
*benchmark* is a program (or a collection of scripts) which performs a
particular calculation in R. The reference for these benchmarks is the
performance of the reference interpreter for the R programming language, namely
`GNU R`_. We sometimes also compare with other interpreters such as `TERR`_ and
`pqR`_.

Benchmarks are used to evaluate a number of metrics:

* can Renjin execute the full benchmark? In other words: is Renjin's
  interpreter sufficiently compatible with GNU R to complete the full
  benchmark?

* how does Renjin's performance compare with the performance of the GNU R
  interpreter?

* where are the computational bottlenecks in the benchmark?

* what is the nature of these bottlenecks and how can this be addressed in
  Renjin?

For Renjin, we are mostly interested in longer running programs which you
typically find in daily situations and less in micro-benchmarks which are
commonly used in computer programming. The :ref:`Bioinformatics` chapter contains a
collection of bio-informatic workflows which represent a wide range of
calculations and some of which require large data sets as their input.

.. _Renjin: http://www.renjin.org
.. _GNU R: https://www.r-project.org
.. _TERR: https://docs.tibco.com/products/tibco-enterprise-runtime-for-r
.. _pqR: http://www.pqr-project.org/

************************
Updates to this document
************************

This is a live document which means that we add and update information as the
project progresses. A recent version of the contents of this document can be
browsed online at http://bedatadriven.github.io/renjin-benchmarks/.

The source code of this document is available at https://github.com/bedatadriven/renjin-benchmarks.

Benchmarks are executed by our automated build system using a variety of
interpreters and the results are reported at
http://packages.renjin.org/benchmarks.

