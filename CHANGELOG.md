# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 2.0.0

This is a complete re-write in Rust. In addition to changing the language, the
the following modifications to the algorithm have been made:

1. The heart of this algorithm is $O(n^2)$. In order to approximate the null distribution,
this $O(n^2)$ operation is performed n times. These permutations are now parallelized and
will scale linearly down to the time it takes to run a single search through the
threshold space. Run time on human data with ~16k genes on 30 CPU is less than an hour
now, and uses less than 3G space.
1. The cmd line input has been greatly simplified. This additionally represents a 
significant change to the protocol:
    - Previously, while in the documentation the operation was described as 'ranking',
    it would be more appropriate to call it sorting or ordering. Ties were not
    addressed. In version 2.0.0, we now expect that the input is a ranked list where
    the first column is the feature identifier, and the second column is the rank. It
    is up to the user to appropriately rank their data. Examples and recommendations
    are provided in the documentation.
1. There are messages printed to stderr when there is more than one set of thresholds
with the same minimum p-value and the same intersect size. This occurs due to the
nature of the hypergeometric p-value.

### Added

1. `profiling/` stores runtime and memory usage information from
[hyperfine](https://github.com/sharkdp/hyperfine) and
[heaptrack](https://github.com/KDE/heaptrack) respectively
1. github actions CI has been added to run tests on pushes to `dev` and `main`
1. Semantic versioning and github releases have been added
1. The package is distributed through crates.io and bioconda
1. Docstrings with examples and module level documentation
1. tests
1. an MPI implementation to parallelize across multiple machines

## 1.0.0 -- Initial release

This version was written by [Yiming Kang](https://github.com/yiming-kang) and is the
version which was used to produce the results in
[Dual threshold optimization and network inference reveal convergent evidence from TF binding locations and TF perturbation responses](https://doi.org/10.1101/gr.259655.119)