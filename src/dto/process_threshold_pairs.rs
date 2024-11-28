//! # Compare ranked feature lists over all threshold pairs.
//! 
//! This is the core logic of the `optimize` function. It performs a double for loop
//! over all threshold pairs. For each threshold pair, it extracts the feature sets
//! and computes the pvalue of the intersection. The asymptotic complexity of this
//! function is $O(n^2)$, where $n$ is the number of thresholds
//! in the ranked feature lists.
use std::collections::HashMap;

use crate::collections::{FeatureSetProvider, PermutedRankedFeatureList, RankedFeatureList};

use crate::stat_operations::{hypergeometric_pvalue, intersect_genes};
use crate::dto::{FeatureSets, OptimizationResultRecord};

/// Performs threshold pair optimization.
///
/// Computes intersections and p-values for all threshold pairs between two ranked feature lists.
///
/// # Arguments
/// - `ranked_feature_list1`: The first ranked feature list.
/// - `ranked_feature_list2`: The second ranked feature list.
/// - `use_permutation`: Whether to use permutation.
/// - `population_size`: The size of the population for p-value calculation.
/// - `debug`: Whether to include feature sets in the result records.
///
/// # Returns
/// A vector of `OptimizationResultRecord`s.
/// 
/// # Example
/// 
/// ```
/// use dual_threshold_optimization::collections::{Feature, FeatureList, PermutedRankedFeatureList, RankedFeatureList};
/// use dual_threshold_optimization::dto::{optimize, OptimizationResult, process_threshold_pairs};
///
/// let genes1 = FeatureList::from(vec![
///     Feature::from("gene1"),
///     Feature::from("gene2"),
///     Feature::from("gene3"),
/// ]);
/// let ranks1 = vec![1, 2, 3];
/// let ranked_feature_list1 = RankedFeatureList::from(genes1, ranks1).unwrap();
///
/// let genes2 = FeatureList::from(vec![
///     Feature::from("gene1"),
///     Feature::from("gene2"),
///     Feature::from("gene4"),
/// ]);
/// let background = FeatureList::from(vec![
///    Feature::from("gene1"),
///    Feature::from("gene2"),
///    Feature::from("gene3"),
///    Feature::from("gene4"),
/// ]);
/// let ranks2 = vec![1, 2, 3];
/// let ranked_feature_list2 = RankedFeatureList::from(genes2, ranks2).unwrap();
/// 
/// let permuted_list1 = PermutedRankedFeatureList::new(&ranked_feature_list1);
/// let permuted_list2 = PermutedRankedFeatureList::new(&ranked_feature_list2);
/// 
/// let results = process_threshold_pairs(
///    &ranked_feature_list1,
///    &ranked_feature_list2,
///    true,
///    4,
///    false,
/// );
/// 
/// assert_eq!(results.len(), 9);
/// 
/// ```
pub fn process_threshold_pairs(
    ranked_feature_list1: &RankedFeatureList,
    ranked_feature_list2: &RankedFeatureList,
    use_permutation: bool,
    population_size: usize,
    debug: bool,
) -> Vec<OptimizationResultRecord> {
    let mut feature_sets_cache = HashMap::new();
    let mut results = Vec::new();

    let permuted_list1 = PermutedRankedFeatureList::new(ranked_feature_list1);
    let permuted_list2 = PermutedRankedFeatureList::new(ranked_feature_list2);

    for &threshold1 in ranked_feature_list1.thresholds().iter() {
        let genes1 = if use_permutation {
            permuted_list1.get_feature_set_by_threshold(threshold1)
        } else {
            ranked_feature_list1.get_feature_set_by_threshold(threshold1)
        };

        for &threshold2 in ranked_feature_list2.thresholds().iter() {
            let genes2 = feature_sets_cache
                .entry(threshold2)
                .or_insert_with(|| {
                    if use_permutation {
                        permuted_list2.get_feature_set_by_threshold(threshold2)
                    } else {
                        ranked_feature_list2.get_feature_set_by_threshold(threshold2)
                    }
            });

            let intersection_size = intersect_genes(&genes1, &genes2);

            let population_size = population_size as u64;
            
            let pvalue = hypergeometric_pvalue(
                population_size,
                genes1.len() as u64,
                genes2.len() as u64,
                intersection_size as u64,
            );

            let feature_sets = if debug {
                FeatureSets::Both(genes1.clone(), genes2.clone())
            } else {
                FeatureSets::None
            };

            results.push(OptimizationResultRecord {
                rank1: threshold1 as usize,
                rank2: threshold2 as usize,
                set1_len: genes1.len(),
                set2_len: genes2.len(),
                population_size,
                intersection_size,
                pvalue,
                permuted: use_permutation,
                feature_sets,
            });
        }
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collections::{Feature, FeatureList, RankedFeatureList};

    #[test]
    fn test_process_threshold_pairs() {
        let genes1 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks1 = vec![1, 2, 3];
        let ranked_feature_list1 = RankedFeatureList::from(genes1, ranks1).unwrap();

        let genes2 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene4"),
        ]);
        let ranks2 = vec![1, 2, 3];
        let ranked_feature_list2 = RankedFeatureList::from(genes2, ranks2).unwrap();

        let results = process_threshold_pairs(
            &ranked_feature_list1,
            &ranked_feature_list2,
            true,
            4,
            false,
        );

        assert!(!results.is_empty());
        for record in results {
            assert!(record.pvalue >= 0.0);
        }
    }

    #[test]
    fn test_process_threshold_pairs_with_and_without_permutation() {
        let genes1 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
            Feature::from("gene4"),
            Feature::from("gene5"),
        ]);
        let ranks1 = vec![1, 2, 3, 4, 5];
        let ranked_feature_list1 = RankedFeatureList::from(genes1, ranks1).unwrap();

        let genes2 = FeatureList::from(vec![
            Feature::from("gene3"),
            Feature::from("gene5"),
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene4"),
        ]);
        let ranks2 = vec![5, 4, 3, 2, 1];
        let ranked_feature_list2 = RankedFeatureList::from(genes2, ranks2).unwrap();

        // Run without permutation
        let results_without_permutation = process_threshold_pairs(
            &ranked_feature_list1,
            &ranked_feature_list2,
            false, // No permutation
            5,
            false,
        );

        // Run with permutation
        let results_with_permutation = process_threshold_pairs(
            &ranked_feature_list1,
            &ranked_feature_list2,
            true, // Use permutation
            5,
            false,
        );

        // Ensure results are different
        assert_ne!(
            results_without_permutation, results_with_permutation,
            "Results with and without permutation should be different."
        );

        // Check basic properties of both results
        assert!(!results_without_permutation.is_empty());
        assert!(!results_with_permutation.is_empty());

        for record in results_without_permutation.iter() {
            assert!(record.pvalue >= 0.0);
        }

        for record in results_with_permutation.iter() {
            assert!(record.pvalue >= 0.0);
        }
    }

}
