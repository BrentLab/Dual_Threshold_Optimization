//! A module for providing a permuted view of a `RankedFeatureList`.
use rand::seq::SliceRandom;
use rand::thread_rng;

use crate::collections::{Feature, FeatureSetProvider, RankedFeatureList};

/// A structure that provides a permuted view of a `RankedFeatureList`.
///
/// This struct uses shuffled indices to permute the genes while keeping the original
/// `RankedFeatureList` unchanged. It supports efficient operations like fetching genes
/// below a given rank threshold.
///
/// # Example
/// ```rust
/// use dual_threshold_optimization::collections::{Feature, FeatureList, FeatureSetProvider, RankedFeatureList, PermutedRankedFeatureList};
///
/// // Create a RankedFeatureList
/// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2"), Feature::from("gene3")]);
/// let ranks = vec![1, 2, 3];
/// let ranked_feature_list = RankedFeatureList::from(genes, ranks).unwrap();
///
/// // Create a permuted view of the RankedFeatureList
/// let permuted_list = PermutedRankedFeatureList::new(&ranked_feature_list);
///
/// // Get genes below a threshold
/// let genes_below_threshold = permuted_list.get_feature_set_by_threshold(2);
///
/// assert_eq!(genes_below_threshold.len(), 2); // Two genes have ranks <= 2
/// ```
pub struct PermutedRankedFeatureList<'a> {
    /// The original `RankedFeatureList`.
    original: &'a RankedFeatureList,
    /// A vector of shuffled indices for the genes in the original list.
    indices: Vec<usize>,
}

impl<'a> PermutedRankedFeatureList<'a> {
    /// Creates a new `PermutedRankedFeatureList` from an existing `RankedFeatureList`.
    ///
    /// This constructor shuffles the indices of the original feature list to provide
    /// a random permutation. The original `RankedFeatureList` remains unchanged.
    ///
    /// # Arguments
    /// - `original`: A reference to the `RankedFeatureList` to be permuted.
    ///
    /// # Example
    /// ```rust
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList, PermutedRankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2"), Feature::from("gene3")]);
    /// let ranks = vec![1, 2, 3];
    /// let ranked_feature_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// let permuted_list = PermutedRankedFeatureList::new(&ranked_feature_list);
    /// ```
    pub fn new(original: &'a RankedFeatureList) -> Self {
        let mut indices: Vec<usize> = (0..original.genes().len()).collect();
        indices.shuffle(&mut thread_rng());
        Self { original, indices }
    }

    pub fn thresholds(&self) -> &[u32] {
        self.original.thresholds()
    }

    pub fn ranks(&self) -> &[u32] {
        &self.original.ranks()
    }
}

impl<'a> FeatureSetProvider for PermutedRankedFeatureList<'a> {
    /// Retrieves genes from the `PermutedRankedFeatureList` with ranks <= the given threshold.
    ///
    /// This implementation ensures the genes are accessed in a randomized order, based on
    /// the permuted indices.
    ///
    /// # Example
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList, PermutedRankedFeatureList, FeatureSetProvider};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2"), Feature::from("gene3")]);
    /// let ranks = vec![1, 2, 3];
    /// let ranked_feature_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// let permuted_list = PermutedRankedFeatureList::new(&ranked_feature_list);
    /// let feature_set = permuted_list.get_feature_set_by_threshold(2);
    ///
    /// assert_eq!(feature_set.len(), 2);
    /// ```
    fn get_feature_set_by_threshold(&self, threshold: u32) -> Vec<Feature> {
        let mut result = Vec::new();

        // iterate along ranks. While teh rank is less than the threshold, add a
        // randomized gene tot he result
        for (rank, index) in self.ranks().iter().zip(self.indices.iter()) {
            if *rank <= threshold {
                result.push(self.original.genes()[*index].clone());
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::collections::{Feature, FeatureList, RankedFeatureList};

    #[test]
    fn test_permuted_ranked_feature_list_with_thresholds() {
        // Create an example RankedFeatureList. Note that it is possible, but highly
        // unlikely that the shuffled indices will result in the same order as the
        // input.
        let genes = FeatureList::from(vec![
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
        let ranks = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let ranked_feature_list = RankedFeatureList::from(genes.clone(), ranks).unwrap();

        // Create a permuted version of the RankedFeatureList
        let permuted = PermutedRankedFeatureList::new(&ranked_feature_list);

        // Fetch genes with ranks <= 2
        let genes_below_threshold = permuted.get_feature_set_by_threshold(10);

        // Assert that the number of genes is correct
        assert_eq!(genes_below_threshold.len(), 10);

        // Assert that the order is not the same as the input order
        let expected_genes: Vec<Feature> = genes.iter().take(10).cloned().collect();
        assert_ne!(
            genes_below_threshold, expected_genes,
            "Features should not be in the same order after permutation"
        );
    }
}
