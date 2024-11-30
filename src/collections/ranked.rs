//! A struct for storing a `FeatureList` along with their corresponding ranks.
//!
//! This module supports `RankedFeatureList` to store ranked genes and offers methods
//! to manage, sort, and analyze ranked feature lists.
use serde::{Deserialize, Serialize};

use crate::collections::{Feature, FeatureList, FeatureSetProvider};

/// A struct that represents a single feature in a `RankedFeatureList`.
/// This struct is used to provide access to the gene, its rank,
/// and its index in the list.
///
/// # Fields
///
/// - `gene`: A reference to the feature.
/// - `rank`: The rank of the feature.
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList, RankedFeatureListItem};
///
/// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
/// let ranks = vec![1, 2];
/// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
///
/// let item = ranked_list.get(1).unwrap();
///
/// assert_eq!(item.feature().id(), "gene2");
/// assert_eq!(item.rank(), 2);
/// ```
#[derive(Debug, PartialEq)]
pub struct RankedFeatureListItem<'a> {
    feature: &'a Feature,
    rank: u32,
    index: usize,
}

impl RankedFeatureListItem<'_> {
    /// Returns a reference to the feature.
    pub fn feature(&self) -> &Feature {
        self.feature
    }
    /// Returns the rank
    pub fn rank(&self) -> u32 {
        self.rank
    }
    /// Returns the index
    pub fn index(&self) -> usize {
        self.index
    }
}

/// ThresholdState is an enum used in RankedFeatureList to track whether there are
/// valid thresholds for the ranks in the list. This is used to determine whether
/// the thresholds need to be recalculated
///
/// # Fields
///
/// - `Unthresholded`: The list has no thresholds
/// - `Thresholded`: The list has thresholds
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::collections::ThresholdState;
///
/// let state = ThresholdState::Unthresholded;
///
/// assert_eq!(state, ThresholdState::Unthresholded);
///
/// let state = ThresholdState::Thresholded;
///
/// assert_eq!(state, ThresholdState::Thresholded);
/// ```
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub enum ThresholdState {
    Unthresholded,
    Thresholded,
}

/// An enum to represent whether a single index or multiple indices should be removed.
/// This is used in the `remove` method of `RankedFeatureList`.
///
/// # Fields
///
/// - `Single(usize)`: Remove a single index
/// - `Multiple(Vec<usize>)`: Remove multiple indices
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::collections::RemoveIndices;
///
/// let single_index = RemoveIndices::Single(5);
/// let multiple_indices = RemoveIndices::Multiple(vec![1, 2, 3]);
/// ```
pub enum RemoveIndices {
    Single(usize),
    Multiple(Vec<usize>),
}

/// A struct that represents a collection of `Feature` instances (`FeatureList`) and
/// their corresponding ranks.
///
/// This struct is designed to store ranked genes and offers utility functions
/// to manage, sort, and analyze ranked feature lists.
///
/// # Fields
///
/// - `genes`: A `FeatureList` containing the `Feature` instances.
/// - `ranks`: A vector of ranks corresponding to each `Feature` in the list.
/// - `thresholds`: Threshold values derived from the ranks for downstream analysis.
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
///
/// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
/// let ranks = vec![1, 2];
/// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
/// assert_eq!(ranked_list.genes().len(), 2);
/// ```
#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub struct RankedFeatureList {
    /// A list of `Feature` instances.
    genes: FeatureList,
    /// A vector of ranks corresponding to the genes.
    ranks: Vec<u32>,
    /// Precomputed thresholds derived from the ranks.
    thresholds: Vec<u32>,
    threshold_state: ThresholdState,
}

impl Default for RankedFeatureList {
    fn default() -> Self {
        Self::new()
    }
}

impl RankedFeatureList {
    /// Create an empty `RankedFeatureList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::RankedFeatureList;
    ///
    /// let ranked_list = RankedFeatureList::new();
    /// assert_eq!(ranked_list.genes().len(), 0);
    /// ```
    pub fn new() -> Self {
        Self {
            genes: FeatureList::new(),
            ranks: Vec::new(),
            thresholds: Vec::new(),
            threshold_state: ThresholdState::Unthresholded,
        }
    }

    /// Create a `RankedFeatureList` from a `FeatureList` and a vector of ranks.
    ///
    /// # Errors
    /// - Returns `Err` if the lengths of `genes` and `ranks` do not match.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    /// ```
    pub fn from(genes: FeatureList, ranks: Vec<u32>) -> Result<Self, String> {
        Self::check_lengths(&genes, &ranks)?;

        let mut ranked_list = Self {
            genes,
            ranks,
            thresholds: Vec::new(),
            threshold_state: ThresholdState::Unthresholded,
        };

        ranked_list.sort_genes_and_ranks();
        ranked_list.generate_thresholds();
        ranked_list.threshold_state = ThresholdState::Thresholded;

        Ok(ranked_list)
    }

    /// Get the genes in this `RankedFeatureList`.
    ///
    /// # Returns
    ///
    /// A slice of `Feature` instances.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// let genes = ranked_list.genes();
    /// assert_eq!(genes.len(), 2);
    /// assert_eq!(genes[0].id(), "gene1");
    /// ```
    pub fn genes(&self) -> &FeatureList {
        &self.genes
    }

    /// Get the ranks in this `RankedFeatureList`.
    ///
    /// # Returns
    ///
    /// A slice of rank values (`u32`).
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// let ranks = ranked_list.ranks();
    /// assert_eq!(ranks.len(), 2);
    /// assert_eq!(ranks[0], 1);
    /// ```
    pub fn ranks(&self) -> &[u32] {
        &self.ranks
    }

    /// Get the thresholds derived from the ranks.
    ///
    /// # Returns
    ///
    /// A slice of threshold values (`u32`).
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    ///
    /// let ranks = vec![1, 2, 3];
    /// let genes = FeatureList::from(vec![
    ///     Feature::from("gene1"),
    ///     Feature::from("gene2"),
    ///     Feature::from("gene3"),
    /// ]);
    /// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// let thresholds = ranked_list.thresholds();
    /// assert_eq!(thresholds[0], 1);
    /// ```
    pub fn thresholds(&self) -> &[u32] {
        &self.thresholds
    }

    /// Retrieve the details of a ranked feature at a specific index.
    ///
    /// This method provides access to a `RankedFeatureListItem` that includes a reference
    /// to the gene, its rank, and its threshold (if the list is thresholded). If the
    /// index is out of bounds, the method returns `None`.
    ///
    /// # Arguments
    /// - `index`: The index of the feature to retrieve.
    ///
    /// # Returns
    /// - `Some(RankedFeatureListItem)` if the index is valid.
    /// - `None` if the index is out of bounds.
    ///
    /// # Notes
    /// - If the `RankedFeatureList` is in an unthresholded state, the `threshold` field
    ///   in the returned `RankedFeatureListItem` will be `None`.
    ///
    /// # Examples
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList, ThresholdState};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// let item = ranked_list.get(1).unwrap();
    /// assert_eq!(item.feature().id(), "gene2");
    /// assert_eq!(item.rank(), 2);
    ///
    /// // Simulate an unthresholded state
    /// ranked_list.remove(0);
    /// let item = ranked_list.get(0).unwrap();
    ///
    /// assert!(ranked_list.get(10).is_none()); // Out-of-bounds index
    /// ```
    pub fn get(&self, index: usize) -> Option<RankedFeatureListItem> {
        if index >= self.genes.len() {
            return None;
        }

        Some(RankedFeatureListItem {
            feature: &self.genes[index],
            rank: self.ranks[index],
            index,
        })
    }

    /// Get the length of the `RankedFeatureList`.
    ///
    /// # Returns
    ///
    /// The number of genes in the list.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    ///
    /// let ranks = vec![1, 2];
    ///
    /// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// assert_eq!(ranked_list.len(), 2);
    pub fn len(&self) -> usize {
        self.genes.len()
    }

    /// Check if the `RankedFeatureList` is empty.
    /// 
    /// # Returns
    /// 
    /// `true` if the list is empty, `false` otherwise.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    /// 
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    /// 
    /// assert!(!ranked_list.is_empty());
    /// 
    /// let empty_list = RankedFeatureList::new();
    /// 
    /// assert!(empty_list.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.genes.is_empty()
    }

    pub fn generate_thresholds(&mut self) {
        let mut thresholds = Vec::new();
        let mut current = 1;
        let max_rank = *self.ranks.last().unwrap_or(&0);

        while current <= max_rank {
            thresholds.push(current);
            current = ((current as f64) * 1.01 + 1.0).floor() as u32;
        }

        // set the final threshold to the maximum rank
        if let Some(last) = self.thresholds.last_mut() {
            *last = max_rank;
        }

        self.thresholds = thresholds;
    }

    /// Removes genes based on a single index or multiple indices.
    ///
    /// This method provides flexibility to remove either a single feature by its index
    /// or multiple genes by a list of indices. Indices are processed in reverse order
    /// to prevent index shifting during removal of multiple genes.
    ///
    /// # Arguments
    ///
    /// - `indices`: Can be a single `usize` or a slice of `usize` (`&[usize]`) indicating
    ///   the indices of the genes to be removed.
    ///
    /// # Returns
    ///
    /// - `Ok(())` if all specified genes are successfully removed.
    /// - `Err(String)` if any of the indices are out of bounds or removal fails.
    ///
    /// # Errors
    ///
    /// - Returns an error if any of the provided indices are out of bounds or removal
    ///   fails for any reason.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![
    ///     Feature::from("gene1"),
    ///     Feature::from("gene2"),
    ///     Feature::from("gene3")
    /// ]);
    /// let ranks = vec![3, 2, 1];
    /// let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// // Remove a single gene
    /// assert!(ranked_list.remove(1).is_ok());
    /// assert_eq!(ranked_list.genes().len(), 2);
    /// assert_eq!(ranked_list.genes()[0].id(), "gene3");
    /// assert_eq!(ranked_list.genes()[1].id(), "gene1");
    ///
    /// // Remove multiple genes
    /// assert!(ranked_list.remove(vec![0, 1]).is_ok());
    /// assert_eq!(ranked_list.genes().len(), 0);
    ///
    /// // Attempt to remove an out-of-bounds index
    /// assert!(ranked_list.remove(10).is_err());
    /// ```
    pub fn remove<T>(&mut self, indices: T) -> Result<(), String>
    where
        T: Into<RemoveIndices>,
    {
        match indices.into() {
            RemoveIndices::Single(index) => self.remove_index(index),
            RemoveIndices::Multiple(mut indices) => {
                indices.sort_unstable_by(|a, b| b.cmp(a)); // Reverse order for safe removal
                for index in indices {
                    self.remove_index(index)?; // Propagate error if any
                }
                Ok(())
            }
        }
    }

    /// Filters the `RankedFeatureList` based on a user-provided closure and removes the matching elements.
    ///
    /// This method applies a filter function to each `RankedFeatureListItem` in the list, collects
    /// the indices of the items for which the filter function returns `true`, and removes those items.
    ///
    /// # Arguments
    /// - `filter_fn`: A closure that takes a reference to a `RankedFeatureListItem` and returns a `bool`.
    ///
    /// # Examples
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![
    ///     Feature::from("gene1"),
    ///     Feature::from("gene2"),
    ///     Feature::from("gene3"),
    /// ]);
    /// let ranks = vec![3, 2, 1];
    /// let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// // Remove all items with a rank less than 3
    /// ranked_list.filter_and_remove(|item| item.rank() < 3);
    ///
    /// assert_eq!(ranked_list.len(), 1);
    /// assert_eq!(ranked_list.get(0).unwrap().feature().id(), "gene1");
    /// ```
    pub fn filter_and_remove<F>(&mut self, filter_fn: F)
    where
        F: Fn(&RankedFeatureListItem) -> bool,
    {
        let indices_to_remove: Vec<usize> = self
            .iter()
            .enumerate()
            .filter_map(|(idx, item)| if filter_fn(&item) { Some(idx) } else { None })
            .collect();

        for idx in indices_to_remove.into_iter().rev() {
            self.remove(idx).unwrap();
        }
    }

    // Private Helper Methods

    /// Removes the feature and its corresponding rank at the specified index.
    ///
    /// This method removes a feature from the `RankedFeatureList` at the given index, ensuring
    /// that the associated rank is also removed. If the feature list is in a thresholded state,
    /// it transitions to an unthresholded state upon successful removal.
    ///
    /// # Arguments
    ///
    /// - `index`: The index of the feature to remove.
    ///
    /// # Returns
    ///
    /// - `Ok(())` if the feature and rank are successfully removed.
    /// - `Err(String)` if the `index` is out of bounds or if the feature removal fails.
    ///
    /// # Errors
    ///
    /// - Returns an error if the provided `index` is greater than or equal to the length of the feature list.
    /// - Returns an error if the removal operation on the `FeatureList` fails for any reason.
    fn remove_index(&mut self, index: usize) -> Result<(), String> {
        if index >= self.genes.len() {
            return Err(format!("Index {} is out of bounds", index));
        }
        self.genes
            .remove(index)
            .map_err(|_| "Feature removal failed".to_string())?;
        self.ranks.remove(index);
        if self.threshold_state == ThresholdState::Thresholded {
            self.threshold_state = ThresholdState::Unthresholded;
        }
        Ok(())
    }

    fn check_lengths(genes: &FeatureList, ranks: &[u32]) -> Result<(), String> {
        if genes.len() != ranks.len() {
            return Err(format!(
                "Feature list length ({}) does not match ranks length ({})",
                genes.len(),
                ranks.len()
            ));
        }
        Ok(())
    }

    fn sort_genes_and_ranks(&mut self) {
        let mut combined: Vec<(u32, &Feature)> = self
            .ranks
            .iter()
            .zip(self.genes.genes())
            .map(|(rank, feature)| (*rank, feature))
            .collect();
        combined.sort_by_key(|(rank, _)| *rank);
        self.ranks = combined.iter().map(|(rank, _)| *rank).collect();
        self.genes = FeatureList::from(
            combined
                .into_iter()
                .map(|(_, feature)| feature.clone())
                .collect::<Vec<_>>(),
        );
    }
}

impl FeatureSetProvider for RankedFeatureList {
    /// Retrieves genes from the `RankedFeatureList` with ranks <= the given threshold.
    ///
    /// # Example
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, FeatureSetProvider, RankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2"), Feature::from("gene3")]);
    /// let ranks = vec![1, 2, 3];
    /// let ranked_feature_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// let feature_set = ranked_feature_list.get_feature_set_by_threshold(2);
    /// assert_eq!(feature_set.len(), 2);
    /// assert!(feature_set.iter().any(|g| g.id() == "gene1"));
    /// assert!(feature_set.iter().any(|g| g.id() == "gene2"));
    /// ```
    fn get_feature_set_by_threshold(&self, threshold: u32) -> Vec<Feature> {
        self.ranks
            .iter()
            .enumerate()
            .filter(|(_, &rank)| rank <= threshold)
            .map(|(idx, _)| self.genes[idx].clone())
            .collect()
    }
}

// Implement `Into<RemoveIndices>` for `usize` and `Vec<usize>` for flexible input.
impl From<usize> for RemoveIndices {
    fn from(index: usize) -> Self {
        RemoveIndices::Single(index)
    }
}

impl From<&[usize]> for RemoveIndices {
    fn from(indices: &[usize]) -> Self {
        RemoveIndices::Multiple(indices.to_vec())
    }
}

impl From<Vec<usize>> for RemoveIndices {
    fn from(indices: Vec<usize>) -> Self {
        RemoveIndices::Multiple(indices)
    }
}

// Iterators for RankedFeatureList
pub struct RankedFeatureListIterator<'a> {
    ranked_feature_list: &'a RankedFeatureList,
    index: usize,
}

impl RankedFeatureList {
    /// Create an iterator for the `RankedFeatureList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// for item in ranked_list.iter() {
    ///     println!("Feature: {}, Rank: {}", item.feature().id(), item.rank());
    /// }
    /// ```
    pub fn iter(&self) -> RankedFeatureListIterator {
        RankedFeatureListIterator {
            ranked_feature_list: self,
            index: 0,
        }
    }
}

impl<'a> Iterator for RankedFeatureListIterator<'a> {
    type Item = RankedFeatureListItem<'a>;

    /// Iterator for the `RankedFeatureList` that yields references to genes and their ranks.
    ///
    /// This iterator allows simultaneous traversal of the genes and their associated ranks
    /// in the `RankedFeatureList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::{Feature, FeatureList, RankedFeatureList};
    ///
    /// let genes = FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();
    ///
    /// for item in ranked_list.iter() {
    ///     println!("Feature: {}, Rank: {}", item.feature().id(), item.rank());
    /// }
    /// ```
    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.ranked_feature_list.genes().len() {
            let item = RankedFeatureListItem {
                feature: &self.ranked_feature_list.genes()[self.index],
                rank: self.ranked_feature_list.ranks()[self.index],
                index: self.index,
            };
            self.index += 1;
            Some(item)
        } else {
            None
        }
    }
}

impl<'a> IntoIterator for &'a RankedFeatureList {
    type Item = RankedFeatureListItem<'a>;
    type IntoIter = RankedFeatureListIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        RankedFeatureListIterator {
            ranked_feature_list: self,
            index: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test the iterator over a `RankedFeatureList` with multiple elements.
    #[test]
    fn test_ranked_feature_list_iterator() {
        let genes = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks = vec![3, 2, 1];
        let ranked_list = RankedFeatureList::from(genes, ranks).unwrap();

        let mut iter = ranked_list.into_iter();
        let item = iter.next().unwrap();
        assert_eq!(item.feature().id(), "gene3");
        assert_eq!(item.rank(), 1);

        let item = iter.next().unwrap();
        assert_eq!(item.feature().id(), "gene2");
        assert_eq!(item.rank(), 2);

        let item = iter.next().unwrap();
        assert_eq!(item.feature().id(), "gene1");
        assert_eq!(item.rank(), 3);

        assert_eq!(iter.next(), None); // This will now work without error
    }

    /// Test that an empty `RankedFeatureList` produces an empty iterator.
    #[test]
    fn test_ranked_feature_list_iterator_empty() {
        let ranked_list = RankedFeatureList::new();
        let mut iter = ranked_list.iter();

        assert_eq!(iter.next(), None);
    }

    /// Test removal of a single valid index.
    #[test]
    fn test_remove_single_index() {
        let genes = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks = vec![3, 2, 1];
        let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();

        assert!(ranked_list.remove(1).is_ok()); // Remove "gene2"
        assert_eq!(ranked_list.genes().len(), 2);

        // Verify remaining genes
        assert_eq!(ranked_list.genes()[0].id(), "gene3");
        assert_eq!(ranked_list.genes()[1].id(), "gene1");
        assert_eq!(ranked_list.ranks(), &[1, 3]);
    }

    /// Test removal of multiple valid indices in reverse order.
    #[test]
    fn test_remove_multiple_indices() {
        let genes = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks = vec![3, 2, 1];
        let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();

        let _ = ranked_list.remove(vec![0, 2]); // Remove "gene1" and "gene3"
        assert_eq!(ranked_list.genes().len(), 1);

        // Verify the remaining gene
        assert_eq!(ranked_list.genes()[0].id(), "gene2");
        assert_eq!(ranked_list.ranks(), &[2]);
    }

    /// Test removing indices out of bounds, which should panic.
    #[test]
    fn test_remove_invalid_index() {
        let genes = FeatureList::from(vec![Feature::from("gene1")]);
        let ranks = vec![1];
        let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();

        let res = ranked_list.remove(5);

        match res {
            Err(e) => assert_eq!(e, "Index 5 is out of bounds"),
            _ => panic!("Expected an error"),
        }
    }

    /// Test that removal of all indices leaves an empty `RankedFeatureList`.
    #[test]
    fn test_remove_all_indices() {
        let genes = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks = vec![3, 2, 1];
        let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();

        let _ = ranked_list.remove(vec![0, 1, 2]); // Remove all indices
        assert!(ranked_list.genes().is_empty());
        assert!(ranked_list.ranks().is_empty());
    }

    /// Test removal of indices with repeated ranks to ensure stability.
    #[test]
    fn test_remove_indices_with_duplicate_ranks() {
        let genes = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks = vec![2, 2, 1];
        let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();

        let _ = ranked_list.remove(vec![0, 2]); // Remove "gene1" and "gene3"
        assert_eq!(ranked_list.genes().len(), 1);

        // Verify the remaining gene
        assert_eq!(ranked_list.genes()[0].id(), "gene1");
        assert_eq!(ranked_list.ranks(), &[2]);
    }

    /// Test the iterator after removals to ensure consistency.
    #[test]
    fn test_iterator_after_removals() {
        let genes = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
        ]);
        let ranks = vec![3, 2, 1];
        let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();

        let _ = ranked_list.remove(vec![1]); // Remove "gene2"

        let mut iter = ranked_list.into_iter();

        // Check the first item
        let item = iter.next().unwrap();
        assert_eq!(item.feature().id(), "gene3");
        assert_eq!(item.rank(), 1);

        // Check the second item
        let item = iter.next().unwrap();
        assert_eq!(item.feature().id(), "gene1");
        assert_eq!(item.rank(), 3);

        // Ensure no more items
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_filter_and_remove() {
        let genes = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
            Feature::from("gene4"),
        ]);
        let ranks = vec![4, 3, 2, 1];
        let mut ranked_list = RankedFeatureList::from(genes, ranks).unwrap();

        // Filter out all items with a rank greater than 2
        ranked_list.filter_and_remove(|item| item.rank() > 2);

        // Check the resulting list
        assert_eq!(ranked_list.len(), 2);
        let remaining_items: Vec<_> = ranked_list.iter().collect();
        assert_eq!(remaining_items[0].feature().id(), "gene4");
        assert_eq!(remaining_items[0].rank(), 1);
        assert_eq!(remaining_items[1].feature().id(), "gene3");
        assert_eq!(remaining_items[1].rank(), 2);
    }

    #[test]
    fn test_filter_by_another_ranked_feature_list() {
        let genes1 = FeatureList::from(vec![
            Feature::from("gene1"),
            Feature::from("gene2"),
            Feature::from("gene3"),
            Feature::from("gene4"),
        ]);
        let ranks1 = vec![4, 3, 2, 1];
        let mut ranked_list1 = RankedFeatureList::from(genes1, ranks1).unwrap();

        let genes2 = FeatureList::from(vec![Feature::from("gene2"), Feature::from("gene4")]);
        let ranks2 = vec![1, 2];
        let ranked_list2 = RankedFeatureList::from(genes2, ranks2).unwrap();

        // Use the second list to filter the first
        ranked_list1.filter_and_remove(|item| {
            !ranked_list2
                .iter()
                .any(|other_item| other_item.feature().id() == item.feature().id())
        });

        // Check that only the genes from the second list remain in the first
        assert_eq!(ranked_list1.len(), 2);
        let remaining_items: Vec<_> = ranked_list1.iter().collect();
        assert_eq!(remaining_items[0].feature().id(), "gene4");
        assert_eq!(remaining_items[1].feature().id(), "gene2");

        // Verify ranks are intact
        assert_eq!(remaining_items[0].rank(), 1);
        assert_eq!(remaining_items[1].rank(), 3);
    }
}
