//! Calculate an empirical p-value for an unpermuted optimization result.
//! 
//! Given a set of `OptimizationResult`s, where exactly one is expected to be
//! unpermuted, calculate a null distribution of p-values from the permuted set
//! and compare the unpermuted p-value to the null distribution to calculate
//! an empirical p-value. If only an unpermuted result is provided, the empirical
//! is set to 1.0 arbitrarily.
use serde_json::json;
use crate::dto::OptimizationResult;
use crate::stat_operations::fdr;

/// Calculates the empirical p-value based on an unpermuted optimization result
/// and a collection of permuted results.
///
/// The empirical p-value is computed as the proportion of permuted p-values that are
/// less than or equal to the p-value of the unpermuted result. The function also
/// returns additional fields, including thresholds, dataset sizes, and population data.
///
/// # Arguments
/// - `unpermuted_result`: The optimization result from the unpermuted data.
/// - `permuted_results`: A vector of optimization results from permuted data.
///
/// # Returns
/// A `serde_json::Value` containing:
/// - `rank1`: The first threshold from the unpermuted result.
/// - `rank2`: The second threshold from the unpermuted result.
/// - `set1_len`: Size of the first dataset.
/// - `set2_len`: Size of the second dataset.
/// - `population_size`: The total size of the population.
/// - `unpermuted_intersection_size`: Size of the intersection for the unpermuted result.
/// - `unpermuted_pvalue`: The p-value of the unpermuted result.
/// - `empirical_pvalue`: The calculated empirical p-value.
///
/// # Panics
/// - Panics if the `unpermuted_result` is not of type `OptimizationResult::Best`.
/// - Panics if any `OptimizationResult` in `permuted_results` is not of type `Best`.
///
/// # Examples
/// ```rust
/// use serde_json::json;
/// use dual_threshold_optimization::dto::{OptimizationResult, OptimizationResultRecord, FeatureSets};
/// use dual_threshold_optimization::stat_operations::empirical_pvalue;
///
/// let unpermuted_result = OptimizationResult::Best(OptimizationResultRecord {
///     rank1: 1,
///     rank2: 2,
///     pvalue: 0.05,
///     set1_len: 50,
///     set2_len: 40,
///     population_size: 100,
///     intersection_size: 10,
///     feature_sets: FeatureSets::None,
///     permuted: false,
/// });
///
/// let input_vec = vec![
///     unpermuted_result.clone(),
///     OptimizationResult::Best(OptimizationResultRecord {
///         rank1: 1,
///         rank2: 2,
///         pvalue: 0.10,
///         set1_len: 50,
///         set2_len: 40,
///         population_size: 100,
///         intersection_size: 8,
///         feature_sets: FeatureSets::None,
///         permuted: true,
///     }),
///     OptimizationResult::Best(OptimizationResultRecord {
///         rank1: 1,
///         rank2: 2,
///         pvalue: 0.03,
///         set1_len: 50,
///         set2_len: 40,
///         population_size: 100,
///         intersection_size: 9,
///         feature_sets: FeatureSets::None,
///         permuted: true,
///     }),
///     OptimizationResult::Best(OptimizationResultRecord {
///         rank1: 1,
///         rank2: 2,
///         pvalue: 0.07,
///         set1_len: 50,
///         set2_len: 40,
///         population_size: 100,
///         intersection_size: 7,
///         feature_sets: FeatureSets::None,
///         permuted: true,
///     }),
/// ];
///
/// let result = empirical_pvalue(input_vec);
/// assert_eq!(
///     result,
///     json!({
///         "rank1": 1,
///         "rank2": 2,
///         "set1_len": 50,
///         "set2_len": 40,
///         "population_size": 100,
///         "unpermuted_intersection_size": 10,
///         "unpermuted_pvalue": 0.05,
///         "empirical_pvalue": 0.3333333333333333, // 1/3 of permuted p-values <= 0.05
///         "fdr": 1.03125,
///     })
/// );
/// ```
pub fn empirical_pvalue(results: Vec<OptimizationResult>) -> serde_json::Value {
    // Separate the unpermuted result and the permuted results
    let (unpermuted_results, permuted_results): (Vec<_>, Vec<_>) = results
        .into_iter()
        .partition(|result| match result {
            OptimizationResult::Best(record) => !record.permuted,
            _ => false,
        });

    // Check the number of unpermuted results
    if unpermuted_results.len() != 1 {
        eprintln!(
            "Warning: Expected exactly one unpermuted result, but found {}.",
            unpermuted_results.len()
        );
    }

    // If there's no unpermuted result, return an empty JSON or error
    if unpermuted_results.is_empty() {
        panic!("No unpermuted result found in the provided results.");
    }

    // Extract the unpermuted result
    let unpermuted_result = match unpermuted_results.into_iter().next().unwrap() {
        OptimizationResult::Best(record) => record,
        _ => panic!("Unexpected result type for unpermuted optimization"),
    };


    let estimated_fdr: f64 = fdr(
        unpermuted_result.set1_len,
        unpermuted_result.set2_len,
        unpermuted_result.intersection_size,
        unpermuted_result.population_size,
        0.8);

    // If there's only one result, set empirical p-value to 1.0 and return
    if permuted_results.is_empty() {
        return json!({
            "rank1": unpermuted_result.rank1,
            "rank2": unpermuted_result.rank2,
            "set1_len": unpermuted_result.set1_len,
            "set2_len": unpermuted_result.set2_len,
            "population_size": unpermuted_result.population_size,
            "unpermuted_intersection_size": unpermuted_result.intersection_size,
            "unpermuted_pvalue": unpermuted_result.pvalue,
            "empirical_pvalue": 1.0,
            "fdr": estimated_fdr,
        });
    }

    // Calculate empirical p-value using only permuted results
    let unpermuted_pvalue = unpermuted_result.pvalue;
    let permuted_pvalues: Vec<f64> = permuted_results
        .into_iter()
        .filter_map(|result| match result {
            OptimizationResult::Best(record) => Some(record.pvalue),
            _ => None,
        })
        .collect();

    let empirical_pvalue = permuted_pvalues
        .iter()
        .filter(|&&p| p <= unpermuted_pvalue)
        .count() as f64
        / permuted_pvalues.len() as f64;
    
    // Return results as JSON
    json!({
        "rank1": unpermuted_result.rank1,
        "rank2": unpermuted_result.rank2,
        "set1_len": unpermuted_result.set1_len,
        "set2_len": unpermuted_result.set2_len,
        "population_size": unpermuted_result.population_size,
        "unpermuted_intersection_size": unpermuted_result.intersection_size,
        "unpermuted_pvalue": unpermuted_pvalue,
        "empirical_pvalue": empirical_pvalue,
        "fdr": estimated_fdr,
    })
}



// pub fn empirical_pvalue(
//     unpermuted_result: OptimizationResult,
//     permuted_results: Vec<OptimizationResult>,
// ) -> serde_json::Value {
//     let (unpermuted_pvalue, unpermuted_rank1, unpermuted_rank2,
//         set1_len, set2_len, population_size, unpermuted_intersect_size) = match unpermuted_result {
//             OptimizationResult::Best(record) => (
//                 record.pvalue,
//                 record.rank1,
//                 record.rank2,
//                 record.set1_len,
//                 record.set2_len,
//                 record.population_size,
//                 record.intersection_size),
//         _ => panic!("Unexpected result type for unpermuted optimization"),
//     };

//     let permuted_pvalues: Vec<f64> = permuted_results
//         .into_iter()
//         .filter_map(|result| match result {
//             OptimizationResult::Best(record) => Some(record.pvalue),
//             _ => None,
//         })
//         .collect();

//     let empirical_pvalue = permuted_pvalues
//         .iter()
//         .filter(|&&p| p <= unpermuted_pvalue)
//         .count() as f64
//         / permuted_pvalues.len() as f64;

//     json!({
//         "rank1": unpermuted_rank1,
//         "rank2": unpermuted_rank2,
//         "set1_len": set1_len,
//         "set2_len": set2_len,
//         "population_size": population_size,
//         "unpermuted_intersection_size": unpermuted_intersect_size,
//         "unpermuted_pvalue": unpermuted_pvalue,
//         "empirical_pvalue": empirical_pvalue,
//     })
// }