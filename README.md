# Dual Threshold Optimization


This library provides a comprehensive toolkit for performing
[Dual Threshold Optimization](https://doi.org/10.1101/gr.259655.119) (DTO)
originally proposed by Kang et al.  

DTO compares two ranked lists of features (e.g. genes) to determine the rank
threshold for each list that minimizes the hypergeometric p-value of the overlap
of features. This method was originally applied to comparing the results of a
paired set of binding assay results and perturbation assay results for a given
transcription factor, but could be used with any two ranked lists of features. 

DTO provides the following statistics on the optimal threshold pair and overlap:
- An empirical p-value of the optimal overlap which is derived from a
  permutation-based null distribution.
- An FDR estimation based on the derivations detailed in the paper linked above.

This crate offers both a library to incorporate DTO into other workflows, 
and a command-line binary for standalone use.

Note: Starting with version 2.0.0, the method was fully re-implemented in Rust. For
the original implementation, which is the version used in the paper linked above,
see version 1.0.0, implemented by Yiming Kang.

## Table of Contents
- [Getting Started](#user-installation)
    - [Using the cmd line](#cmdline-usage)
        - [Output](#output)
    - [Using the library](#library-usage)
    - [Development](#developer-installation-and-usage)
- [Algorithmic Details](#algorithmic-details)
- [Troubleshooting](#troubleshooting)

## Getting started

Binaries for Linux, MacOS and Windows are provided in in the `release` tab. There are
two flavors of releases for each OS:

1. **The standard release**: this is suitable for almost every user

2. **An MPI enabled version**: only if you want to parallelize across multiple machines.
**NOTE**: For this version, MPI must be installed on the host system.


### Using the cmd line

With the correct binary, you can print the help message like so:

```bash
dual_threshold_optimization --help
```

The following will be provided. Please note that you can find examples of the input
lists and background here:

- input list examples: [list1](test_data/ranklist1.csv), [list2](test_data/ranklist2.csv)
- background example: [background](test_data/background.txt)

```bash
Dual Threshold Optimization CLI

Usage: dual_threshold_optimization [OPTIONS] --ranked-list1 <FILE> --ranked-list2 <FILE>

Options:
  -1, --ranked-list1 <FILE>
          Path to the first ranked feature list (CSV format).
          
          This should have two columns: "feature" and "rank". There should be **NO HEADER**.
          
          Rank is expected to be an integer. It is recommended that ties are handed with the `min` or `max` method.

  -2, --ranked-list2 <FILE>
          Path to the second ranked feature list (CSV format)
          
          This should have two columns: "feature" and "rank". There should be **NO HEADER**.
          
          Rank is expected to be an integer. It is recommended that ties are handed with the `min` or `max` method.

  -b, --background <FILE>
          Path to the background feature list (one feature per line, optional)

  -p, --permutations <PERMUTATIONS>
          Number of permutations to perform
          
          [default: 1000]

  -t, --threads <THREADS>
          Number of threads to use per task. For single-node parallelization, this will be the number of threads available on the machine.
          
          For multi-node parallelization, this will be the number of threads per task.
          
          Example: If you submit via Slurm with `-ntasks 4 --cpus-per-task 10`, then this value should be set to 10. This configuration will run 40 permutations in parallel.
          
          [default: 1]

  -m, --multi-node
          Enable multi-node mode using MPI. This requires that the program has been built with the `mpi` feature enabled

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version

```

#### Output

The output from the cmd line is a json to stdout. To redirect this to a file, you 
would do the following:

```bash
dual_threshold_optimization -1 list1.csv -2 $list2.csv -p 1000 -t 24 > output.json
```

This is what the output will look like:

```json
{
    "empirical_pvalue": 0.0,
    "fdr": 0.3481243441762854,
    "population_size": 14295,
    "rank1": 1016,
    "rank2": 896,
    "set1_len": 1016,
    "set2_len": 896,
    "unpermuted_intersection_size": 127,
    "unpermuted_pvalue": 1.5719992077514072e-14
}
```

Where the fields are the following:

- **empirical_pvalue**: The quantile of the unpermuted minimum p-value in relation to
the series of permuted minimal p-values
- **fdr**: The false discovery rate where the sensitivity is set to 0.8. See the
[DTO paper](https://doi.org/10.1101/gr.259655.119) for more details
- **population_size**: The size of the background. If no background is explicity
provided, this is the length of the input lists (when no background is provided, 
the lists must contain the same set of features)
- **rank1**: The optimal rank for the unpermuted minimum p-value for list 1
- **rank2**: The optimal rank for the unpermuted minimum p-value for list 2
- **set1_len**: The number of features with rank less than or equal to **rank1** in
list1
- **set2_len**: The number of features with rank less than or equal to **rank2** in
list2
- **unpermuted_intersection_size**: The number of genes in the intersection of
list1 and list2 with rank less than or equal to their respective optimal ranks
- **unpermuted p-value** the optimal p-value of the unpermuted lists. For analysis
purposes, the empirical p-value should be used.

### Using the library

To use the library, you can `cargo add dual_threshold_optimization` in your rust
project. See the crates.io documentation for more information about what is provided
in each of the submodules.

### Developer installation and usage

It is assumed that you have the
[rust toolchain](https://www.rust-lang.org/tools/install) already installed.

1. git pull this repository
2. `cd` into the repo

For any of the commands below, you can add `--features mpi` to include the MPI
feature. But, remember that this requires that MPI exist in your environment
(e.g. [openMPI](https://www.open-mpi.org/))

At this point, you can run the tests with:

```bash
cargo test
```

or

```bash
cargo test
```

you can run the binary with

```bash
cargo run -- --help
```

and you can guild with

```bash
cargo build
```

Note that there is a build profile for profiling which will build a release version
with the debug flags on:

```bash
cargo build --profile release-debug
```

## Algorithmic details

The following provides details on the DTO algorithm, step by step.

1. Initialize two ranked feature lists

    Begin with two ranked lists of features, e.g. genes, where each feature has
    an id, e.g. a unique identifier for the gene, and a rank. The rank must be an
    integer and is expected to have ties handled with a method such as "min" or "max"
    where ties all are assigned the same rank.
    
1. Create a series of thresholds for each list based on the ranks

    For each list, generate a series of thresholds T1, T2, ... . These thresholds are
    used to generate sets of features from each list to compare the overlap.
    The thresholds are calculated by the recurrence relation

    $$
    T_1 = 1 \\
    Tn = Floor(T_{n-1} * 1.01 + 1)
    $$
    
    The stopping condition is when the threshold meets or exceeds the largest rank.
    The final threshold is always set to the max rank. This series provides finer
    spacing at higher ranks, allowing more granular selection among top-ranked genes.

1. Conduct a brute force search of the threshold pairs to find an optimal overlap

    For each possible pair of thresholds (one from each list’s threshold series),
    select the genes from each list that rank above the respective threshold. Calculate
    the hypergeometric p-value by intersecting the feature sets

1. Select optimal threshold pair

    Track the threshold pair that produces the minimum P-value across all tested pairs.
    This threshold combination is considered optimal for identifying significant
    overlap between the two lists.

    **CAVEAT**: We have discovered that the minimal p-value may not be unique. There
    are possibly multiple sets that yield the same p-value, including the minimal
    p-value. When this occurs on the minimal p-value, the threshold pair that yields
    the largest overlap is selected. When there are multiple threshold pairs that
    have the same p-value and the same intersect size, the first in the set is
    chosen arbitrarily.

1. Use permutations to generate a null distribution for the minimal p-value

    To assess the statistical significance of the identified overlap, run DTO multiple
    times (e.g., 1000 runs) on randomized versions of the ranked lists. This creates a
    null distribution of the minimal p-value. This null distribution allows for
    evaluating the observed minimum P-value relative to random chance.

1. Calculate false discovery rate (FDR)

    In the [DTO paper](https://doi.org/10.1101/gr.259655.119), an FDR is derived. This
    FDR is estimated for the optimal threshold pair.:


## Troubleshooting

If you are using the MPI binary, then you must have MPI in your environment. If you
do have MPI installed, but you get an error similar to the one below:

```bash
./dual_threshold_optimization: error while loading shared libraries: libmpi.so.40: cannot open shared object file: No such file or directory
```

Then you need to find where the `libmpi.so.40` file lives and add
it to your `LD_LIBRARY_PATH` manually. E.g.

```bash
export LD_LIBRARY_PATH=/ref/mblab/software/spack-0.22.2/opt/spack/linux-rocky9-x86_64/gcc-11.4.1/openmpi-5.0.3-vjscapwoywmullqs3lj2mmdf7vyge4rk/lib:$LD_LIBRARY_PATH
```