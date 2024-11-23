use rand::seq::SliceRandom;
use rand::thread_rng;

use crate::{hypergeometric_pvalue, intersect_genes, Gene, GeneList, RankedGeneList};
#[derive(Debug, Clone)]
pub enum GeneSets {
    Both(Vec<Gene>, Vec<Gene>),
    None,
}

/// A struct to record the result of dual threshold optimization.
#[derive(Debug, Clone)]
pub struct OptimizationResultRecord {
    pub rank1: usize,
    pub rank2: usize,
    pub intersection_size: usize,
    pub pvalue: f64,
    pub permuted: bool,
    pub gene_sets: GeneSets,
}

pub enum OptimizationResult {
    Debug(Vec<OptimizationResultRecord>), // All results
    Best(OptimizationResultRecord),       // Single best result
}

/// Perform dual threshold optimization.
///
/// This function optimizes thresholds for overlap between two ranked gene lists
/// while computing p-values for the observed intersection. It optionally supports
/// debug mode to store all intermediate results, including gene sets.
///
/// # Arguments
/// - `ranked_gene_list1`: The first `RankedGeneList` to optimize.
/// - `ranked_gene_list2`: The second `RankedGeneList` to optimize.
/// - `permute`: Whether to shuffle the ranked gene lists (useful for generating null distributions).
/// - `background`: An optional `GeneList` defining the background population. If absent,
///   the two ranked lists must have identical genes.
/// - `debug`: Whether to store and return all intermediate results for analysis.
///
/// # Returns
/// - An `OptimizationResult` containing:
///   - If `debug` is `true`, all intermediate results are returned as a `Vec` of `OptimizationResultRecord`s
///     under the `OptimizationResult::Debug` variant. These include the gene sets for each threshold pair.
///   - If `debug` is `false`, the single best `OptimizationResultRecord` is returned under the
///     `OptimizationResult::Best` variant, determined by p-value and intersection size.
///
/// # Panics
/// - If `background` is not provided and the two ranked gene lists do not have identical genes.
/// - If `background` is provided and either ranked gene list contains genes not in the background.
///
/// # Example
/// ```rust
/// use dual_threshold_optimization::{dual_threshold_optimization, Gene, GeneList, OptimizationResult, RankedGeneList};
///
/// let genes1 = GeneList::from(vec![
///     Gene::from("gene1"),
///     Gene::from("gene2"),
///     Gene::from("gene3"),
/// ]);
/// let ranks1 = vec![1, 2, 3];
/// let ranked_gene_list1 = RankedGeneList::from(genes1, ranks1).unwrap();
///
/// let genes2 = GeneList::from(vec![
///     Gene::from("gene1"),
///     Gene::from("gene2"),
///     Gene::from("gene4"),
/// ]);
/// let background = GeneList::from(vec![
///    Gene::from("gene1"),
///    Gene::from("gene2"),
///    Gene::from("gene3"),
///    Gene::from("gene4"),
/// ]);
/// let ranks2 = vec![1, 2, 3];
/// let ranked_gene_list2 = RankedGeneList::from(genes2, ranks2).unwrap();
///
/// let result = dual_threshold_optimization(&ranked_gene_list1, &ranked_gene_list2, false, Some(&background), true);
///
/// match result {
///     OptimizationResult::Debug(results) => {
///         assert!(!results.is_empty());
///         for record in results {
///             println!("Rank1: {}, Rank2: {}, p-value: {}", record.rank1, record.rank2, record.pvalue);
///         }
///     }
///     OptimizationResult::Best(best) => {
///         println!("Best result: Rank1 {}, Rank2 {}, p-value {}", best.rank1, best.rank2, best.pvalue);
///         assert!(best.pvalue > 0.0);
///     }
/// }
/// ```
pub fn dual_threshold_optimization<'a>(
    ranked_gene_list1: &RankedGeneList,
    ranked_gene_list2: &RankedGeneList,
    permute: bool,
    background: Option<&GeneList>,
    debug: bool,
) -> OptimizationResult {
    // Set the population_size based on the background
    let population_size = match background {
        // If background is None, then verify that that ranked_gene_list1 and
        // ranked_gene_list2 have the same set of genes. If they do not, panic!
        // otherwise, use the length of the intersection (which is the length) of
        // either ranked_gene_list1 or ranked_gene_list2) as the population size.
        None => {
            // Calculate intersection size
            let intersection_size = ranked_gene_list1
                .genes()
                .intersect(ranked_gene_list2.genes())
                .len();

            if intersection_size != ranked_gene_list1.len()
                || intersection_size != ranked_gene_list2.len()
            {
                panic!("If no background is provided, then the gene lists must have the same genes. Filter the gene lists and try again.");
            }

            intersection_size
        }
        // If background is Some, then verify that the ranked_gene_list1
        // and ranked_gene_list2 do not contain genes that are not in the background.
        // If they do, panic! otherwise, use the length of the background as
        // the population size.
        Some(bg) => {
            let list1_diff = ranked_gene_list1.genes().difference(bg, false);
            if list1_diff.len() > 0 {
                panic!(
                    "The following genes in the first ranked gene list are not in the background set: {:?}",
                    list1_diff.iter().map(|g| g.id()).collect::<Vec<_>>()
                );
            }

            let list2_diff = ranked_gene_list2.genes().difference(bg, false);
            if list2_diff.len() > 0 {
                panic!(
                    "The following genes in the second ranked gene list are not in the background set: {:?}",
                    list2_diff.iter().map(|g| g.id()).collect::<Vec<_>>()
                );
            }

            bg.genes().len()
        }
    };

    // Optionally permute the indices. This serves to permute the genes in relation to
    // the ranks. This is used over x iterations, eg 1000, to calculate the null
    // distribution of p-values
    let mut rng = thread_rng();
    let mut index1: Vec<usize> = (0..ranked_gene_list1.genes().len()).collect();
    let mut index2: Vec<usize> = (0..ranked_gene_list2.genes().len()).collect();
    if permute {
        index1.shuffle(&mut rng);
        index2.shuffle(&mut rng);
    }

    // extract a reference to the
    let mut results = Vec::new();

    // Double loop over thresholds
    for &threshold1 in ranked_gene_list1.thresholds().iter() {
        // Select genes with ranks <= threshold1
        let genes1: Vec<Gene> = index1
            .iter()
            .filter(|&&idx| ranked_gene_list1.ranks()[idx] <= threshold1)
            .map(|&idx| ranked_gene_list1.get(idx).unwrap().gene().clone())
            .collect();

        for &threshold2 in ranked_gene_list2.thresholds().iter() {
            // Select genes with ranks <= threshold2
            let genes2: Vec<Gene> = index2
                .iter()
                .filter(|&&idx| ranked_gene_list2.ranks()[idx] <= threshold2)
                .map(|&idx| ranked_gene_list2.get(idx).unwrap().gene().clone())
                .collect();

            // Calculate intersection size
            let intersection_size = intersect_genes(&genes1, &genes2);

            // Calculate hypergeometric p-value. note that the
            // `successes_in_population` and `sample_size` are symmetric -- it doesn't
            // matter which is which.
            let pvalue = hypergeometric_pvalue(
                population_size as u64,
                genes1.len() as u64,
                genes2.len() as u64,
                intersection_size as u64,
            );

            let gene_sets = if debug {
                GeneSets::Both(
                    genes1
                        .iter()
                        .map(|gene| (*gene).clone())
                        .collect::<Vec<Gene>>(),
                    genes2
                        .iter()
                        .map(|gene| (*gene).clone())
                        .collect::<Vec<Gene>>(),
                )
            } else {
                GeneSets::None
            };

            results.push(OptimizationResultRecord {
                rank1: threshold1 as usize,
                rank2: threshold2 as usize,
                intersection_size,
                pvalue,
                permuted: permute,
                gene_sets,
            });
        }
    }

    if debug {
        // Return all results if debug is true
        return OptimizationResult::Debug(results);
    } else {
        // Filter for the best result based on current criteria
        let min_pvalue = results
            .iter()
            .map(|r| r.pvalue)
            .fold(f64::INFINITY, f64::min);
        let mut best_results: Vec<_> = results
            .into_iter()
            .filter(|res| res.pvalue == min_pvalue)
            .collect();

        // if the best pvalue is not unique, then select the thresholds that yield the
        // largest intersection.
        // If the best p-value is not unique, then select the thresholds that yield the
        // largest intersection.
        if best_results.len() > 1 {
            println!(
                "Multiple results with the same minimum p-value ({:.15}). \
                Choosing the result with the largest intersection size.",
                min_pvalue
            );

            let max_intersect = best_results
                .iter()
                .map(|res| res.intersection_size)
                .max()
                .unwrap();

            best_results = best_results
                .into_iter()
                .filter(|res| res.intersection_size == max_intersect)
                .collect();

            // Notify the user if ties remain after filtering by intersection size
            if best_results.len() > 1 {
                println!(
                    "Multiple results with the same maximum intersection size ({}). \
                    Choosing an arbitrary result based on the order of thresholds.",
                    max_intersect
                );
            }

            // If the intersection is also not unique, then choose the threshold that
            // is first (sorted by rank1 and rank2 lexicographically).
            best_results.sort_by_key(|res| (res.rank1, res.rank2));
        }

        // Return the single best result
        OptimizationResult::Best(best_results.into_iter().next().unwrap())
    }
}
