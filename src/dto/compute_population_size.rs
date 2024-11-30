//! # Compute the population size of two sets
//!
//! `compute_population_size` validates the compatibility of the two
//! `RankedFeatureList`s and, optionally, a background `FeatureList`. If the background
//! is not provided, then the two ranked feature lists must have identical feature sets.
//! If the background is provided, then the two lists must be subsets of the background.
use crate::collections::{FeatureList, RankedFeatureList};

/// Computes the population size and validates the compatibility of the background.
///
/// If a background is provided, it ensures that all features in the ranked feature lists
/// are present in the background. Otherwise, it checks that the two ranked feature lists
/// have identical feature sets.
///
/// # Arguments
/// - `ranked_feature_list1`: The first ranked feature list.
/// - `ranked_feature_list2`: The second ranked feature list.
/// - `background`: Optional background feature list.
///
/// # Returns
/// The population size as `usize`.
///
/// # Panics
/// - If the background is provided and features in the ranked lists are missing from it.
/// - If no background is provided and the ranked feature lists have differing feature sets.
///
/// # Example
/// ```
/// use dual_threshold_optimization::collections::{RankedFeatureList, FeatureList, Feature};
/// use dual_threshold_optimization::dto::compute_population_size;
///
/// // RankedFeatureList 1
/// let features1 = FeatureList::from(vec![
///     Feature::from("gene1"),
///     Feature::from("gene2"),
///     Feature::from("gene3"),
/// ]);
/// let ranks1 = vec![1, 2, 3];
/// let ranked_feature_list1 = RankedFeatureList::from(features1, ranks1).unwrap();
///
/// // RankedFeatureList 2
/// let features2 = FeatureList::from(vec![
///     Feature::from("gene1"),
///     Feature::from("gene2"),
///     Feature::from("gene3"),
/// ]);
/// let ranks2 = vec![1, 2, 3];
/// let ranked_feature_list2 = RankedFeatureList::from(features2, ranks2).unwrap();
///
/// // Background feature list
/// let background = FeatureList::from(vec![
///     Feature::from("gene1"),
///     Feature::from("gene2"),
///     Feature::from("gene3"),
///     Feature::from("gene4"),
/// ]);
///
/// // Compute population size with a background
/// let pop_size = compute_population_size(&ranked_feature_list1, &ranked_feature_list2, Some(&background));
/// assert_eq!(pop_size, 4);
///
/// // Compute population size without a background (identical feature sets required)
/// let pop_size = compute_population_size(&ranked_feature_list1, &ranked_feature_list2, None);
/// assert_eq!(pop_size, 3);
/// ```
pub fn compute_population_size(
    ranked_feature_list1: &RankedFeatureList,
    ranked_feature_list2: &RankedFeatureList,
    background: Option<&FeatureList>,
) -> usize {
    match background {
        None => {
            let intersection_size = ranked_feature_list1
                .genes()
                .intersect(ranked_feature_list2.genes())
                .len();

            if intersection_size != ranked_feature_list1.len()
                || intersection_size != ranked_feature_list2.len()
            {
                panic!(
                    "If no background is provided, the feature lists must have identical genes."
                );
            }
            intersection_size
        }
        Some(bg) => {
            for (ranked_list, list_name) in [
                (ranked_feature_list1, "first"),
                (ranked_feature_list2, "second"),
            ] {
                let diff = ranked_list.genes().difference(bg, false);
                if !diff.is_empty() {
                    panic!(
                        "The following genes in the {} ranked feature list are not in the background: {:?}",
                        list_name,
                        diff.iter().map(|g| g.id()).collect::<Vec<_>>()
                    );
                }
            }
            bg.genes().len()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collections::{Feature, FeatureList, RankedFeatureList};

    #[test]
    #[should_panic(
        expected = "If no background is provided, the feature lists must have identical genes."
    )]
    fn test_compute_population_size_no_background_mismatched_lists() {
        let features1 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks1 = vec![1, 2, 3];
        let ranked_feature_list1 = RankedFeatureList::from(features1, ranks1).unwrap();

        let features2 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene4"), // Different feature
        ]);
        let ranks2 = vec![1, 2, 3];
        let ranked_feature_list2 = RankedFeatureList::from(features2, ranks2).unwrap();

        // No background and mismatched lists should panic
        compute_population_size(&ranked_feature_list1, &ranked_feature_list2, None);
    }

    #[test]
    #[should_panic(
        expected = "The following genes in the first ranked feature list are not in the background:"
    )]
    fn test_compute_population_size_missing_genes_in_background_list1() {
        let features1 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks1 = vec![1, 2, 3];
        let ranked_feature_list1 = RankedFeatureList::from(features1, ranks1).unwrap();

        let features2 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks2 = vec![1, 2, 3];
        let ranked_feature_list2 = RankedFeatureList::from(features2, ranks2).unwrap();

        let background = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            // Missing gene3
        ]);

        compute_population_size(
            &ranked_feature_list1,
            &ranked_feature_list2,
            Some(&background),
        );
    }

    #[test]
    #[should_panic(
        expected = "The following genes in the second ranked feature list are not in the background:"
    )]
    fn test_compute_population_size_missing_genes_in_background_list2() {
        let features1 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks1 = vec![1, 2, 3];
        let ranked_feature_list1 = RankedFeatureList::from(features1, ranks1).unwrap();

        let features2 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene4"), // Extra gene
        ]);
        let ranks2 = vec![1, 2, 3];
        let ranked_feature_list2 = RankedFeatureList::from(features2, ranks2).unwrap();

        let background = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);

        compute_population_size(
            &ranked_feature_list1,
            &ranked_feature_list2,
            Some(&background),
        );
    }

    #[test]
    fn test_compute_population_size_valid_background() {
        let features1 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks1 = vec![1, 2, 3];
        let ranked_feature_list1 = RankedFeatureList::from(features1, ranks1).unwrap();

        let features2 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks2 = vec![1, 2, 3];
        let ranked_feature_list2 = RankedFeatureList::from(features2, ranks2).unwrap();

        let background = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
            Feature::from("gene4"), // Extra gene in the background
        ]);

        let population_size = compute_population_size(
            &ranked_feature_list1,
            &ranked_feature_list2,
            Some(&background),
        );
        assert_eq!(population_size, 4); // Background size
    }

    #[test]
    fn test_compute_population_size_no_background_identical_lists() {
        let features1 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks1 = vec![1, 2, 3];
        let ranked_feature_list1 = RankedFeatureList::from(features1, ranks1).unwrap();

        let features2 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks2 = vec![1, 2, 3];
        let ranked_feature_list2 = RankedFeatureList::from(features2, ranks2).unwrap();

        let population_size =
            compute_population_size(&ranked_feature_list1, &ranked_feature_list2, None);
        assert_eq!(population_size, 3); // Size of the feature lists
    }
}
