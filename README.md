[![Build Status](http://build.renjin.org/job/Benchmark/job/GNUR-GCE/badge/icon)](http://build.renjin.org/job/Benchmark/job/GNUR-GCE/)

# Renjin Benchmarks

This repository contains a collection of bioinformatics workflows/pipelines
which we use to benchmark Renjin's performance.

## Documentation

These benchmarks are documented using
[Sphinx](http://www.sphinx-doc.org/en/stable/index.html). Use the `Makefile` in
the root of the repository to build the documentation. For example, `make html`
(executed in the root of the project) will generate HTML documentation in
a `build/html` folder inside the `docs` folder.

The Sphinx configuration file is `conf.py` and the index of the documentation
is in the `index.rst` file.
