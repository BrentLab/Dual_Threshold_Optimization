use crate::collections::RankedFeatureList;
use crate::dto::{process_threshold_pairs, OptimizationResult};

pub fn optimize<'a>(
    ranked_feature_list1: &RankedFeatureList,
    ranked_feature_list2: &RankedFeatureList,
    permute: bool,
    population_size: usize,
    debug: bool,
) -> OptimizationResult {
    let results = process_threshold_pairs(
        ranked_feature_list1,
        ranked_feature_list2,
        permute,
        population_size,
        debug,
    );

    if debug {
        // Return all results if debug is true
        OptimizationResult::Debug(results)
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
            eprintln!(
                "Multiple results with the same minimum p-value ({:.15}). \
                Choosing the result with the largest intersection size.",
                min_pvalue
            );

            let max_intersect = best_results
                .iter()
                .map(|res| res.intersection_size)
                .max()
                .unwrap();

            best_results.retain(|res| res.intersection_size == max_intersect);

            // Notify the user if ties remain after filtering by intersection size
            if best_results.len() > 1 {
                eprintln!(
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

#[cfg(test)]
mod tests {
    use crate::collections::{Feature, FeatureList, RankedFeatureList};
    use crate::dto::{optimize, OptimizationResult};

    #[test]
    fn test_full_optimization_with_permutation() {
        // First ranked feature list with 10 genes
        let genes1 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
            Feature::from("gene4"),
            Feature::from("gene5"),
            Feature::from("gene6"),
            Feature::from("gene7"),
            Feature::from("gene8"),
            Feature::from("gene9"),
            Feature::from("gene10"),
        ]);
        let ranks1 = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let ranked1 = RankedFeatureList::from(genes1, ranks1).unwrap();

        // Second ranked feature list with the same genes but different order
        let genes2 = FeatureList::from(vec![
            Feature::from("gene7"),
            Feature::from("gene3"),
            Feature::from("gene9"),
            Feature::from("gene1"),
            Feature::from("gene5"),
            Feature::from("gene10"),
            Feature::from("gene4"),
            Feature::from("gene8"),
            Feature::from("gene6"),
            Feature::from("gene2"),
        ]);
        let ranks2 = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let ranked2 = RankedFeatureList::from(genes2, ranks2).unwrap();

        let population_size: usize = 10;

        // Run the optimization function with permutation enabled and no background
        let result = optimize(&ranked1, &ranked2, false, population_size, false);

        // Ensure the optimization returned a valid best result
        assert!(matches!(result, OptimizationResult::Best(_)));

        if let OptimizationResult::Best(best) = result {
            assert_eq!(best.rank1, 3);
            assert_eq!(best.rank2, 4);
            assert_eq!(best.pvalue, 0.33333333333333337);
            println!(
                "Best result: Rank1 {}, Rank2 {}, P-value {}",
                best.rank1, best.rank2, best.pvalue
            );
        }

        // Run the optimization function with permutation enabled and no background
        let result_permuted = optimize(&ranked1, &ranked2, true, population_size, false);

        // Ensure the optimization returned a valid best result
        assert!(matches!(result_permuted, OptimizationResult::Best(_)));

        if let OptimizationResult::Best(best) = result_permuted {
            // assert_ne!(best.rank1, 2);
            // assert_ne!(best.rank2, 1);
            // assert_ne!(best.pvalue, 0.2);
            println!(
                "Best result: Rank1 {}, Rank2 {}, P-value {}",
                best.rank1, best.rank2, best.pvalue
            );
        }
    }
}
