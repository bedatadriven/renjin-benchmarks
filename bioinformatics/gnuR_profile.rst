Profiling benchmarks with GNU R
===============================


Profiling environment
---------------------

We used the following environment to test the profile and analyse the benchmarks with GNU R.

+-----------------------+------------------------------------------------------+
| R version	        | R version 3.3.2 (2016-10-31); svn rev 71607          |
+-----------------------+------------------------------------------------------+
| "proftools" package	| 0.99-2                                               |
+-----------------------+------------------------------------------------------+
| perf_events version	| 4.4.35                                               |
+-----------------------+------------------------------------------------------+
| Machine	        | Intel® Core™ i5-3317U CPU @ 1.70GHz × 4; 16 GB (2x8) |
+-----------------------+------------------------------------------------------+
| Operating System	| Linux 4.4.0-57-generic; Ubuntu 14.04.1 x86_64        |
+-----------------------+------------------------------------------------------+



GNU R instrumentation results
-----------------------------

We instrumented GNU R by addition of timers to functions that are used for calling external
C, C++, and Fortran code (.C, .Call, and .External functions). This information provides us
with an overview of where the most time is spent in each workflow.

.. figure:: ../docs/_static/results_instrumentation.pdf
   :scale: 100 %
   :alt: Results instrumentation GNU R
   :figwidth: 75 %

   Running of benchmarks on instrumented GNU R. (Left) Distribution of time spent in 
   executing R code or native code through  .C/.Fortran, .Call, 
   .External interfaces. (Right) Running time of each of the benchmarks; red dashed
   line is 10% trimmed mean of all running times.

   +---+-------------------------+----+-------------------------+----+-------------------------+
   | # | Benchmark name          | #  | Benchmark name          | #  | Benchmark name          |
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 1 | affy                    | 6  | integration_livercohort | 11 | simulated_geo_matrix    | 
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 2 | clinical_eslII          | 7  | microarray              | 12 | single_cell_lun         | 
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 3 | clinical_livercohort    | 8  | mutation                | 13 | survival_simple         | 
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 4 | generate_count          | 9  | rnaseq_deseq2           | 14 | survival_tcga           |
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 5 | integration_igraph      | 10 | rppa                    | 15 | tcga_browser            | 
   +---+-------------------------+----+-------------------------+----+-------------------------+

The running time of the workflows ranges from 0.5 second (clinical_livercohort) to 37 minutes
(generate_count). In workflows with long running times, most of time is spent in native code such as C/Fortran
(rppa, survival_tcga) or C++ (generate_count, rnaseq_deseq2, singlecell_lun). Indeed, it has become a common
practice to rewrite slow parts of R code in native languages. However, having to rewrite an algorithm (partially)
in other languages is an extra burden for scientists in terms of learning new programming languages and maintaining 
such hybrid code base.

CPU Profiling with GNU R
------------------------

Analysis run fastest when the number of operations can be reduced and by making optimal use of the available (hardware) 
resources (CPU, memory, etc). Therefore there are two kind of optimizations that can be built into the interpreter, 
reducing the number of operations (overhead), and improving usage of for example CPUs. Therefore we have analyzed the 
CPU, and CPU cache memory usage while running each of the workflows. 

.. figure:: ../docs/_static/results_cpu_performance.pdf
   :scale: 80 %
   :alt: Results instrumentation GNU R

   Analysis of the CPU and CPU cache memory usage by benchmarks. (Top left) Instructions per CPU cycle (IPC) with 
   red dashed line indicating the acceptable efficiency threshold of 1 IPC. (Top right) Stalled cycles per instruction (SCPI)
   indicating how many operations waited for necessary data to be loaded in the cache. Red dashed line is the mean SCPI (trimmed 
   10%). (Bottom left) Percentage of L1-dcache-load-misses in each benchmark with red dashed line indicating 10% trimmed mean of
   all benchmarks. (Bottom right) The ratio of number of lookups in slowest cache with number of L1-dcache-load-misses, indicating
   what part of data that were not found in fastest cache memory were not found in intermediate caches and had to be looked up in 
   the slowest cache (LLC).

   +---+-------------------------+----+-------------------------+----+-------------------------+
   | # | Benchmark name          | #  | Benchmark name          | #  | Benchmark name          |
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 1 | affy                    | 6  | integration_livercohort | 11 | simulated_geo_matrix    | 
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 2 | clinical_eslII          | 7  | microarray              | 12 | single_cell_lun         | 
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 3 | clinical_livercohort    | 8  | mutation                | 13 | survival_simple         | 
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 4 | generate_count          | 9  | rnaseq_deseq2           | 14 | survival_tcga           |
   +---+-------------------------+----+-------------------------+----+-------------------------+
   | 5 | integration_igraph      | 10 | rppa                    | 15 | tcga_browser            | 
   +---+-------------------------+----+-------------------------+----+-------------------------+


For each operation the required data is loaded in a super fast memory on the CPU chip named L1-cache. If the data 
necessary for the operation is not found, the CPU searches on the slower (but larger) caches such as L2, L3, and LCC 
(each level is slower but larger than previous one). If the data is not found on L1 cache and need to be looked up 
in other caches, this is called L1-dcache-load-misses.

We observed between 2-81% L1-cache-load-misses while running the workflows with lowest cache-misses being observed 
in ‘tcga_browser’ and the highest in ‘rppa’. In addition, between 13-71% of L1-dcache-load-misses had to be looked 
up in the slowest LCC cache. All this suggest 
significant gains can be achieved by optimizing processor cache usage. The metrics instructions per CPU cycle (IPC) 
and stalled cycles per instruction (SCPI) indicates how efficient the processor cycles are being used. Ideally you 
want to reduce SCPI and increase the IPC. An IPC lower than 1 is considered inefficient code. This is the case for 
only one of the workflows (integration_igraph) with ‘rppa’ being on the border (1.01 instructions per cycle). Likewise 
‘integration_igraph’ and ‘rppa’ have the highest SCPI values.


Summary
----------
The collected benchmarks cover a variety of analysis commonly performed in the field of biomedical research. The 
benchmarks are highly variable in their running time, efficient use of resources, size of input data, and implementation 
of algorithms in native languages such as C, C++, and Fortran. These benchmarks could, therefore, assist in identification 
of computational bottlenecks, causes of inefficient uses of available hardware, incompatibility with the standard R 
R interpreter. In addition, the efficiently running benchmarks could help us to detect regressions during development 
of Renjin. Improvement of IPC/SCPI, running time, and reduction CPU cache misses are valuable metrics for this purpose.
