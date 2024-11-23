use rand::seq::SliceRandom;
use rand::thread_rng;

use crate::{Gene, GeneList, hypergeometric_pvalue, intersect_genes, RankedGeneList};

/// A struct to record the result of dual threshold optimization.
#[derive(Debug, Clone)]
pub struct OptimizationResult {
    pub i: usize,
    pub j: usize,
    pub intersection_size: usize,
    pub pvalue: f64,
    pub permuted: bool,
}


/// Perform dual threshold optimization.
///
/// # Arguments
/// - `ranked_gene_list1`: A reference to the first `RankedGeneList`.
/// - `ranked_gene_list2`: A reference to the second `RankedGeneList`.
/// - `permute`: Whether to permute the indices of the genes in the ranked lists.
///
/// # Returns
/// - The `OptimizationResult` with the minimum p-value.
///
/// # Examples
/// ```
/// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList, dual_threshold_optimization};
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
/// 
/// let background = GeneList::from(vec![
///    Gene::from("gene1"),
///    Gene::from("gene2"),
///    Gene::from("gene3"),
///    Gene::from("gene4"),
/// ]);
/// let ranks2 = vec![1, 2, 3];
/// let ranked_gene_list2 = RankedGeneList::from(genes2, ranks2).unwrap();
///
/// let result = dual_threshold_optimization(&ranked_gene_list1, &ranked_gene_list2, false, Some(&background));
/// assert_eq!(result.i, 1);
/// assert_eq!(result.j, 1);
/// assert_eq!(result.intersection_size, 2);
/// assert_eq!(result.pvalue, 0.16666666666666666);
/// ```
pub fn dual_threshold_optimization(
    ranked_gene_list1: &RankedGeneList,
    ranked_gene_list2: &RankedGeneList,
    permute: bool,
    background: Option<&GeneList>,
) -> OptimizationResult {
    // If the background set is Some(), check that all genes in the RankedGeneLists
    // are in the background set. If not, raise an error. In the error, print the
    // genes that are not in the background set. Tell the user that they can omit
    // the background set to use the intersection of the RankedGeneLists as the
    // the population.
    let population_size = match background {
        None => {

            // Calculate intersection size
            let intersection_size = ranked_gene_list1.genes().intersect(ranked_gene_list2.genes()).len();

            if intersection_size != ranked_gene_list1.len() || intersection_size != ranked_gene_list2.len() {
                panic!("If no background is provided, then the gene lists must have the same genes. Filter the gene lists and try again.");
            }

            intersection_size

        },
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

    // Optionally permute the indices
    let mut rng = thread_rng();
    let mut index1: Vec<usize> = (0..ranked_gene_list1.genes().len()).collect();
    let mut index2: Vec<usize> = (0..ranked_gene_list2.genes().len()).collect();
    if permute {
        index1.shuffle(&mut rng);
        index2.shuffle(&mut rng);
    }

    let thresholds1 = ranked_gene_list1.thresholds();
    let thresholds2 = ranked_gene_list2.thresholds();
    let mut results = Vec::new();

    // Double loop over thresholds
    for (i, &threshold1) in thresholds1.iter().enumerate() {
        // Select genes with ranks <= threshold1
        let genes1: Vec<Gene> = index1
            .iter()
            .filter(|&&idx| ranked_gene_list1.ranks()[idx] <= threshold1)
            .map(|&idx| ranked_gene_list1.get(idx).unwrap().gene().clone())
            .collect();

        for (j, &threshold2) in thresholds2.iter().enumerate() {
            // Select genes with ranks <= threshold2
            let genes2: Vec<Gene> = index2
                .iter()
                .filter(|&&idx| ranked_gene_list2.ranks()[idx] <= threshold2)
                .map(|&idx| ranked_gene_list2.get(idx).unwrap().gene().clone())
                .collect();

            // Calculate intersection size
            let intersection_size = intersect_genes(&genes1, &genes2);

            // Calculate hypergeometric p-value
            let pvalue = hypergeometric_pvalue(
                population_size as u64,                       // Total population size
                genes1.len() as u64, // Successes in population (this and sample_size are symmetric)
                genes2.len() as u64,             // Sample size
                intersection_size as u64,   // Observed overlap
            );

            println!(
                "i: {} j: {} population_size: {}, successes_in_pop: {}, sample_size: {}, intersect: {}, pvalue: {}",
                i, j, ranked_gene_list1.genes().len(), genes1.len(), genes2.len(), intersection_size, pvalue
            );

            // Record result
            results.push(OptimizationResult {
                i,
                j,
                intersection_size,
                pvalue,
                permuted: permute,
            });
        }
    }

    // Find the result with the minimum p-value
    let min_pvalue = results.iter().map(|r| r.pvalue).fold(f64::INFINITY, f64::min);
    let mut best_results: Vec<_> = results.into_iter().filter(|res| res.pvalue == min_pvalue).collect();

    // If there are multiple results with the same p-value, choose the one with the greatest intersection size
    if best_results.len() > 1 {
        let max_intersect = best_results
            .iter()
            .map(|res| res.intersection_size)
            .max()
            .unwrap();
        best_results = best_results
            .into_iter()
            .filter(|res| res.intersection_size == max_intersect)
            .collect();

        // If there's still more than one result, raise an error
        if best_results.len() > 1 {
            panic!("Multiple results with the same p-value and intersection size. Unable to determine the best result.");
        }
    }

    // If there are multiple results with the same p-value, choose the one with the greatest intersection size
    if best_results.len() > 1 {
        let max_intersect = best_results
            .iter()
            .map(|res| res.intersection_size)
            .max()
            .unwrap();
        best_results = best_results
            .into_iter()
            .filter(|res| res.intersection_size == max_intersect)
            .collect();

        // If there's still more than one result, raise an error
        if best_results.len() > 1 {
            best_results.sort_by_key(|res| (res.i, res.j));
            println!("Multiple results with the same p-value and intersection size. Returning the first one.");
            return best_results.into_iter().next().unwrap();
        }
    }

    best_results.into_iter().next().unwrap()
}
