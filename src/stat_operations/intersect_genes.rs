//! Find the intersection between two slices of `Feature` objects.
use crate::collections::Feature;
use std::collections::HashSet;

/// Computes the size of the intersection of two slices of `&Feature` based on their IDs.
///
/// This function takes two slices of references to `Feature` objects and calculates
/// the number of `Feature` objects that are present in both input slices. The comparison
/// is performed using the unique `id` of each `Feature`.
///
/// # Arguments
///
/// - `genes1`: A slice of references to `Feature` objects.
/// - `genes2`: Another slice of references to `Feature` objects.
///
/// # Returns
///
/// The size of the intersection between the two input slices, as a `usize`.
///
/// # Examples
///
/// ```rust
/// use dual_threshold_optimization::collections::Feature;
/// use dual_threshold_optimization::stat_operations::intersect_genes;
///
/// let gene1 = Feature::from("gene1");
/// let gene2 = Feature::from("gene2");
/// let gene3 = Feature::from("gene3");
/// let gene4 = Feature::from("gene4");
///
/// let genes1 = vec![gene1.clone(), gene2.clone(), gene3.clone()];
/// let genes2 = vec![gene2.clone(), gene3.clone(), gene4.clone()];
///
/// let intersection_size = intersect_genes(&genes1, &genes2);
///
/// assert_eq!(intersection_size, 2);
/// ```
pub fn intersect_genes<'a, T1, T2>(genes1: T1, genes2: T2) -> usize
where
    T1: AsRef<[Feature]>,
    T2: AsRef<[Feature]>,
{
    // Collect IDs from the first list into a HashSet
    let ids1: HashSet<_> = genes1
        .as_ref()
        .iter()
        .map(|feature: &Feature| feature.id())
        .collect();

    // Count the intersection with the second list
    genes2
        .as_ref()
        .iter()
        .filter(|feature| ids1.contains(feature.id()))
        .count()
}
