//! # A structure for storing and manipulating a list of `Feature` instances.
//!
//! The `FeatureList` struct provides methods for managing and validating a
//! collection of `Feature` structs which may represent, for example, genes or
//! transcripts.
use std::collections::HashMap;
use std::collections::HashSet;
use std::ops::Index;

use serde::{Deserialize, Serialize};

use crate::collections::{Feature, UniqueCheck};

/// A collection of `Feature` objects, with methods for manipulation and validation.
///
/// The `FeatureList` struct provides methods for managing and validating a collection of genes.
/// It ensures that operations like adding, removing, or checking for unique identifiers are simple and efficient.
///
/// # Fields
///
/// - `genes`: A vector of `Feature` instances.
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::collections::{Feature, FeatureList};
///
/// let genes = vec![Feature::from("gene1"), Feature::from("gene2")];
/// let feature_list = FeatureList::from(genes);
/// assert_eq!(feature_list.len(), 2);
/// ```
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct FeatureList {
    genes: Vec<Feature>,
}

impl Default for FeatureList {
    fn default() -> Self {
        Self::new()
    }
}

impl FeatureList {
    /// Create an empty `FeatureList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::FeatureList;
    ///
    /// let feature_list = FeatureList::new();
    /// assert_eq!(feature_list.len(), 0);
    /// ```
    pub fn new() -> Self {
        Self { genes: Vec::new() }
    }

    /// Create a `FeatureList` from an initial vector of genes.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let genes = vec![Feature::from("gene1"), Feature::from("gene2")];
    /// let feature_list = FeatureList::from(genes);
    /// assert_eq!(feature_list.len(), 2);
    /// ```
    pub fn from(genes: Vec<Feature>) -> Self {
        Self { genes }
    }

    /// Add a `Feature` to the list.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let mut feature_list = FeatureList::new();
    /// feature_list.push(Feature::from("gene1"));
    /// assert_eq!(feature_list.len(), 1);
    /// ```
    pub fn push(&mut self, feature: Feature) {
        self.genes.push(feature);
    }

    /// Remove the last `Feature` from the list and return it.
    /// Returns `None` if the list is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let mut feature_list = FeatureList::new();
    /// feature_list.push(Feature::from("gene1"));
    /// let popped_feature = feature_list.pop();
    /// assert_eq!(popped_feature.unwrap().id(), "gene1");
    /// assert!(feature_list.is_empty());
    /// ```
    pub fn pop(&mut self) -> Option<Feature> {
        self.genes.pop()
    }

    /// Check if all `Feature` IDs in the list are unique.
    ///
    /// - Returns `UniqueCheck::Unique` if all IDs are unique.
    /// - Returns `UniqueCheck::Duplicates` with a map of duplicate IDs and their indices if duplicates are found.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, UniqueCheck};
    ///
    /// let genes = vec![Feature::from("gene1"), Feature::from("gene2"), Feature::from("gene1")];
    /// let feature_list = FeatureList::from(genes);
    ///
    /// match feature_list.check_unique() {
    ///     UniqueCheck::Unique => println!("All IDs are unique."),
    ///     UniqueCheck::Duplicates(duplicates) => {
    ///         for (id, indices) in duplicates {
    ///             println!("Duplicate ID: {}, Indices: {:?}", id, indices);
    ///         }
    ///     }
    /// }
    /// ```
    pub fn check_unique(&self) -> UniqueCheck {
        let mut seen: HashMap<&str, Vec<usize>> = HashMap::new();
        for (index, feature) in self.genes.iter().enumerate() {
            seen.entry(feature.id()).or_default().push(index);
        }

        let duplicates: HashMap<String, Vec<usize>> = seen
            .into_iter()
            .filter(|(_, indices)| indices.len() > 1)
            .map(|(id, indices)| (id.to_string(), indices))
            .collect();

        if duplicates.is_empty() {
            UniqueCheck::Unique
        } else {
            UniqueCheck::Duplicates(duplicates)
        }
    }

    /// Get a reference to the internal list of genes.
    ///
    /// # Returns
    ///
    /// A reference to the vector of `Feature` instances stored in the `FeatureList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let genes = vec![Feature::from("gene1"), Feature::from("gene2")];
    /// let feature_list = FeatureList::from(genes);
    ///
    /// let internal_genes = feature_list.genes();
    /// assert_eq!(internal_genes.len(), 2);
    /// assert_eq!(internal_genes[0].id(), "gene1");
    /// ```
    pub fn genes(&self) -> &[Feature] {
        &self.genes
    }

    /// Get the number of genes in the list.
    ///
    /// # Returns
    ///
    /// The total number of `Feature` instances in the `FeatureList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let genes = vec![Feature::from("gene1"), Feature::from("gene2")];
    /// let feature_list = FeatureList::from(genes);
    ///
    /// assert_eq!(feature_list.len(), 2);
    /// ```
    pub fn len(&self) -> usize {
        self.genes.len()
    }

    /// Check if the feature list is empty.
    ///
    /// # Returns
    ///
    /// `true` if the `FeatureList` contains no `Feature` instances, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let empty_feature_list = FeatureList::new();
    /// assert!(empty_feature_list.is_empty());
    ///
    /// let genes = vec![Feature::from("gene1")];
    /// let non_empty_feature_list = FeatureList::from(genes);
    /// assert!(!non_empty_feature_list.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.genes.is_empty()
    }

    /// Remove a `Feature` from the list by index.
    /// Returns the removed `Feature` if successful,
    /// or an error message if the index is out of bounds.
    ///
    /// # Arguments
    ///
    /// * `index` - The index of the `Feature` to remove.
    ///
    /// # Returns
    ///
    /// - `Ok(Feature)` if the `Feature` was successfully removed.
    /// - `Err(String)` if the index is out of bounds.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let mut feature_list = FeatureList::from(vec![
    ///    Feature::from("gene1"),
    ///   Feature::from("gene2"),
    ///   Feature::from("gene3"),
    /// ]);
    ///
    /// let removed_feature = feature_list.remove(1);
    /// match removed_feature {
    ///     Ok(feature) => println!("Removed feature: {}", feature.id()),
    ///     Err(msg) => println!("Error: {}", msg.to_string()),
    /// }
    /// assert_eq!(feature_list.len(), 2);
    /// assert_eq!(feature_list.genes()[0].id(), "gene1");
    /// assert_eq!(feature_list.genes()[1].id(), "gene3");
    /// ```
    pub fn remove(&mut self, index: usize) -> Result<Feature, String> {
        if index >= self.genes.len() {
            return Err(format!(
                "Index out of bounds: {} (length: {})",
                index,
                self.genes.len()
            ));
        }
        Ok(self.genes.remove(index))
    }

    /// Computes the intersection of two `FeatureList`s.
    ///
    /// This method identifies genes that are present in both `FeatureList`s
    /// based on their IDs and returns a new `FeatureList` containing the common genes.
    ///
    /// # Arguments
    /// - `other`: The other `FeatureList` to intersect with.
    ///
    /// # Returns
    /// A new `FeatureList` containing the intersection of the two `FeatureList`s.
    ///
    /// # Examples
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let list1 = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let list2 = FeatureList::from(vec![Feature::from("gene2"), Feature::from("gene3")]);
    ///
    /// let intersection = list1.intersect(&list2);
    ///
    /// assert_eq!(intersection.len(), 1);
    /// assert_eq!(intersection[0].id(), "gene2");
    /// ```
    pub fn intersect<'a>(&'a self, other: &FeatureList) -> Vec<&'a Feature> {
        let ids1: HashSet<_> = self.genes().iter().map(|g| g.id()).collect();
        let ids2: HashSet<_> = other.genes().iter().map(|g| g.id()).collect();

        self.genes()
            .iter()
            .filter(|feature| ids1.contains(feature.id()) && ids2.contains(feature.id()))
            .collect()
    }

    /// Computes the genes that are either exclusive to one of the two `FeatureList`s
    /// (symmetric difference) or only in `self` but not in `other` (set difference).
    ///
    /// # Arguments
    /// - `other`: The other `FeatureList` to compare against.
    /// - `symmetric`: If `true`, computes the symmetric difference (genes in one list but not both).
    ///                If `false`, computes the set difference (genes in `self` but not in `other`).
    ///
    /// # Returns
    /// A new `FeatureList` containing the resulting genes.
    ///
    /// # Examples
    /// ## Symmetric Difference:
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let list1 = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let list2 = FeatureList::from(vec![Feature::from("gene2"), Feature::from("gene3")]);
    ///
    /// let exclusive = list1.difference(&list2, true);
    ///
    /// assert_eq!(exclusive.len(), 2);
    /// assert!(exclusive.iter().any(|g| g.id() == "gene1"));
    /// assert!(exclusive.iter().any(|g| g.id() == "gene3"));
    /// ```
    /// ## Set Difference:
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let list1 = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let list2 = FeatureList::from(vec![Feature::from("gene2"), Feature::from("gene3")]);
    ///
    /// let exclusive = list1.difference(&list2, false);
    ///
    /// assert_eq!(exclusive.len(), 1);
    /// assert_eq!(exclusive[0].id(), "gene1");
    /// ```
    pub fn difference<'a>(&'a self, other: &'a FeatureList, symmetric: bool) -> Vec<&'a Feature> {
        let ids1: HashSet<_> = self.genes().iter().map(|g| g.id()).collect();
        let ids2: HashSet<_> = other.genes().iter().map(|g| g.id()).collect();

        if symmetric {
            self.genes()
                .iter()
                .chain(other.genes().iter())
                .filter(|feature| {
                    ids1.symmetric_difference(&ids2)
                        .any(|&id| id == feature.id())
                })
                .collect()
        } else {
            self.genes()
                .iter()
                .filter(|feature| ids1.difference(&ids2).any(|&id| id == feature.id()))
                .collect()
        }
    }

    /// Returns an iterator over the features in the `FeatureList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList};
    ///
    /// let feature_list = FeatureList::from(vec![
    ///     Feature::from("gene1"),
    ///     Feature::from("gene2"),
    /// ]);
    ///
    /// for feature in feature_list.iter() {
    ///     println!("Feature: {}", feature.id());
    /// }
    /// ```
    pub fn iter(&self) -> FeatureListIterator {
        FeatureListIterator {
            feature_list: self,
            index: 0,
        }
    }
}

pub struct FeatureListIterator<'a> {
    feature_list: &'a FeatureList,
    index: usize,
}

impl<'a> Iterator for FeatureListIterator<'a> {
    type Item = &'a Feature;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.feature_list.len() {
            let item = &self.feature_list.genes[self.index];
            self.index += 1;
            Some(item)
        } else {
            None
        }
    }
}

impl<'a> IntoIterator for &'a FeatureList {
    type Item = &'a Feature;
    type IntoIter = FeatureListIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl From<Vec<String>> for FeatureList {
    fn from(ids: Vec<String>) -> Self {
        let genes: Vec<Feature> = ids
            .into_iter()
            .map(|id| Feature::from(id.as_str()))
            .collect();
        FeatureList { genes }
    }
}

impl From<Vec<&str>> for FeatureList {
    fn from(ids: Vec<&str>) -> Self {
        let genes: Vec<Feature> = ids.into_iter().map(Feature::from).collect();
        FeatureList { genes }
    }
}

impl Index<usize> for FeatureList {
    type Output = Feature;

    fn index(&self, index: usize) -> &Self::Output {
        &self.genes[index]
    }
}
