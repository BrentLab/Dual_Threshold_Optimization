//! Traits for working with ranked feature lists.
use crate::collections::Feature;

/// A trait to retrieve a set of genes with ranks below or equal to a given threshold.
///
/// This trait provides a unified interface for accessing genes based on a rank threshold,
/// allowing the same logic to apply to different representations of ranked feature lists,
/// such as `RankedFeatureList` or `PermutedRankedFeatureList`.
pub trait FeatureSetProvider {
    fn get_feature_set_by_threshold(&self, threshold: u32) -> Vec<Feature>;
}
