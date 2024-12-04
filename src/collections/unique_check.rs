//! A utility for checking if a list of features is unique
use std::collections::HashMap;

/// Enum to represent the result of the uniqueness check for feature IDs in a `FeatureList`.
///
/// - `Unique`: Indicates all feature IDs are unique.
/// - `Duplicates`: Contains a map where the key is a duplicated feature ID, and the value is a vector of indices where the duplicates occur.
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::collections::{Feature, FeatureList, UniqueCheck};
///
/// let genes = vec![
///     Feature::from("gene1"),
///     Feature::from("gene2"),
///     Feature::from("gene1"),
/// ];
/// let feature_list = FeatureList::from(genes);
///
/// match feature_list.check_unique() {
///     UniqueCheck::Unique => println!("All feature IDs are unique."),
///     UniqueCheck::Duplicates(duplicates) => {
///         for (id, indices) in duplicates {
///             println!("Duplicate ID: {}, Indices: {:?}", id, indices);
///         }
///     }
/// }
/// ```
#[derive(Debug, PartialEq)]
pub enum UniqueCheck {
    Unique,
    Duplicates(HashMap<String, Vec<usize>>),
}
