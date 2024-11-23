use std::collections::HashMap;

/// Enum to represent the result of the uniqueness check for gene IDs in a `GeneList`.
///
/// - `Unique`: Indicates all gene IDs are unique.
/// - `Duplicates`: Contains a map where the key is a duplicated gene ID, and the value is a vector of indices where the duplicates occur.
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::{Gene, GeneList, UniqueCheck};
///
/// let genes = vec![
///     Gene::from("gene1"),
///     Gene::from("gene2"),
///     Gene::from("gene1"),
/// ];
/// let gene_list = GeneList::from(genes);
///
/// match gene_list.check_unique() {
///     UniqueCheck::Unique => println!("All gene IDs are unique."),
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
