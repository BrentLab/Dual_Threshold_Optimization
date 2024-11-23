# Dual Threshold Optimization


1. Initialize two ranked gene lists

    Begin with two ranked lists of genes, where each gene has an id and a value. These
    lists are ordered by value in descending order, with the highest values at the top
    of each list.
    
    **NOTE**: How *should* ties be handled? Current implementation is
    *sorted* **not** ranked. Ties are not considered.
    **NOTE**: Actually rank -- rank method doesn't matter  (i don't think) -- treat all
    that are tied as a single class, meaning the threshold takes all of them

1. Generate threshold series for each list

    For each list, generate a series of thresholds T1, T2, ..., which determine where
    to “cut off” the list. The thresholds are calculated by the recurrence relation

    $$
    T_1 = 1 \\
    Tn = Floor(T_{n-1} * 1.01 + 1)
    $$
    
    This series provides finer spacing at higher ranks, allowing more granular
    selection among top-ranked genes.

    **NOTE**: Implement this as the cartesian product of the threshold vectors

    **NOTE**: diff lengths do not matter

    **NOTE**: offer method of placing length limit on series based on, eg, binding
    value (do not examine thresholds below 0.1 for instance). This should be possible
    for both lists, eg cut the binding list where binding is  > 0.1 and cut the
    perturbation list where log2fc < 0.5 (which also requires that the ranking be 
    configurable as either asc/desc)

1. Select genes above each threshold pair

    For each possible pair of thresholds (one from each list’s threshold series),
    select the genes from each list that rank above the respective threshold. These
    genes are marked as “positives” in their respective lists. This process creates a
    subset of genes from each list that are the “top-ranked” genes according to the
    current threshold pair.

1. Calculate overlap and hypergeometric p-value between positive gene sets

    For each threshold pair, compute the overlap between the “positive” gene sets from
    each list (genes that are above the thresholds in both lists). This overlap
    represents genes present in both sets at the given thresholds.

    **NOTE**: this is basically one sided -- so the null is that the size of the
    intersect is smaller than or equal to random

1. Compute hypergeometric p-value for overlap

    For each overlap calculated in the previous step, compute a hypergeometric P-value.
    This value quantifies the statistical significance of the overlap, indicating how
    likely it is that the observed overlap would occur by chance.

1. Select optimal threshold pair

    Track the threshold pair that produces the minimum P-value across all tested pairs.
    This threshold combination is considered optimal for identifying significant
    overlap between the two lists.

1. Generate null distribution of p-values (this can be done simultaneously)

    To assess the statistical significance of the identified overlap, run DTO multiple
    times (e.g., 1000 runs) on randomized versions of the ranked lists. This creates a
    null distribution of P-values. This null distribution allows for evaluating the
    observed minimum P-value relative to random chance.

1. Calculate false discovery rate (FDR)

    Using the null distribution, compute the False Discovery Rate (FDR) to estimate the
    likelihood that the observed overlap is a true signal rather than a false positive.
    The FDR is derived based on the expected fraction of genes that would overlap due to
    non-functional binding, as modeled by the formula provided in the derivation
    (see paper)

    **NOTE**: quantile, from the empirical null, of the optimal p-value

1. Return optimal threshold pair and significance:

    The final output is the threshold pair that minimizes the nominal P-value, along with
    the significance of this P-value as assessed against the null distribution. This result
    indicates the gene sets in each list that overlap most significantly, with an estimated
    FDR to quantify confidence in the result.

    **NOTE**: provide option to return output of the FDR equation by user provided
    Sn (default to 0.8)
    **NOTE**: Nominal pvalue does not matter. Can output, but the valuable number is
    the quantile from the empirical null

## Software Design

### Objects/functions

- Gene -- struct with one attr: name

- GeneList -- array of Genes

- RankList -- array of tuples (rank, ptr) wh

- ThresholdGrid -- array of 

- DTO (output is a pvalue and label of either input/empirical_null)

### Process

1. Read in/check input sets
    - Input: user
    - Output: 2 GeneList, 2 RankList
1. Optionally permute the 
1. Double for loop over RankList1 and RankList2
1. Generate n_random_permutations of both GeneLists
1. create array of (GeneList1, Genelist2, ThresholdGrid), ...
    - These can be run in parallel