use crate::{Gene, GeneList};
#[derive(Debug, PartialEq)]
pub struct RankedGeneListItem<'a> {
    gene: &'a Gene,
    rank: u32,
    threshold: Option<u32>,
    index: usize,
}

impl<'a> RankedGeneListItem<'a> {
    /// Returns a reference to the gene.
    pub fn gene(&self) -> &Gene {
        self.gene
    }
    /// Returns the rank
    pub fn rank(&self) -> u32 {
        self.rank
    }
    /// Returns the threshold
    pub fn threshold(&self) -> Option<u32> {
        self.threshold
    }
    /// Returns the index
    pub fn index(&self) -> usize {
        self.index
    }
}
/// ThresholdState is an enum used in RankedGeneList to track whether there are
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
/// use dual_threshold_optimization::ThresholdState;
/// 
/// let state = ThresholdState::Unthresholded;
/// 
/// assert_eq!(state, ThresholdState::Unthresholded);
/// 
/// let state = ThresholdState::Thresholded;
/// 
/// assert_eq!(state, ThresholdState::Thresholded);
/// ```
#[derive(Debug, PartialEq)]
pub enum ThresholdState {
    Unthresholded,
    Thresholded,
}

/// An enum to represent whether a single index or multiple indices should be removed.
/// This is used in the `remove` method of `RankedGeneList`.
/// 
/// # Fields
/// 
/// - `Single(usize)`: Remove a single index
/// - `Multiple(Vec<usize>)`: Remove multiple indices
/// 
/// # Examples
/// 
/// ```
/// use dual_threshold_optimization::RemoveIndices;
/// 
/// let single_index = RemoveIndices::Single(5);
/// let multiple_indices = RemoveIndices::Multiple(vec![1, 2, 3]);
/// ```
pub enum RemoveIndices {
    Single(usize),
    Multiple(Vec<usize>),
}

/// Iterator for the `RankedGeneList` that yields references to genes and their ranks.
///
/// This iterator allows simultaneous traversal of the genes and their associated ranks
/// in the `RankedGeneList`.
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
///
/// let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
/// let ranks = vec![1, 2];
/// let ranked_list = RankedGeneList::from(genes, ranks).unwrap();
///
/// for item in ranked_list.iter() {
///     println!("Gene: {}, Rank: {}", item.gene().id(), item.rank());
/// }
/// ```
pub struct RankedGeneListIterator<'a> {
    ranked_gene_list: &'a RankedGeneList,
    index: usize,
}

impl<'a> Iterator for RankedGeneListIterator<'a> {
    type Item = RankedGeneListItem<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.ranked_gene_list.genes().len() {
            let item = RankedGeneListItem {
                gene: &self.ranked_gene_list.genes()[self.index],
                rank: self.ranked_gene_list.ranks()[self.index],
                threshold: if self.ranked_gene_list.threshold_state == ThresholdState::Thresholded {
                    self.ranked_gene_list.thresholds().get(self.index).copied()
                } else {
                    None
                },
                index: self.index,
            };
            self.index += 1;
            Some(item)
        } else {
            None
        }
    }
}

impl<'a> IntoIterator for &'a RankedGeneList {
    type Item = RankedGeneListItem<'a>;
    type IntoIter = RankedGeneListIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        RankedGeneListIterator {
            ranked_gene_list: self,
            index: 0,
        }
    }
}

/// A struct that represents a collection of `Gene` instances (`GeneList`) and
/// their corresponding ranks.
///
/// This struct is designed to store ranked genes and offers utility functions
/// to manage, sort, and analyze ranked gene lists.
///
/// # Fields
///
/// - `genes`: A `GeneList` containing the `Gene` instances.
/// - `ranks`: A vector of ranks corresponding to each `Gene` in the list.
/// - `thresholds`: Threshold values derived from the ranks for downstream analysis.
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
///
/// let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
/// let ranks = vec![1, 2];
/// let ranked_list = RankedGeneList::from(genes, ranks).unwrap();
/// assert_eq!(ranked_list.genes().len(), 2);
/// ```
#[derive(Debug, PartialEq)]
pub struct RankedGeneList {
    /// A list of `Gene` instances.
    genes: GeneList,
    /// A vector of ranks corresponding to the genes.
    ranks: Vec<u32>,
    /// Precomputed thresholds derived from the ranks.
    thresholds: Vec<u32>,
    threshold_state: ThresholdState,
}

impl RankedGeneList {
    /// Create an empty `RankedGeneList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::RankedGeneList;
    ///
    /// let ranked_list = RankedGeneList::new();
    /// assert_eq!(ranked_list.genes().len(), 0);
    /// ```
    pub fn new() -> Self {
        Self {
            genes: GeneList::new(),
            ranks: Vec::new(),
            thresholds: Vec::new(),
            threshold_state: ThresholdState::Unthresholded,
        }
    }

    /// Create a `RankedGeneList` from a `GeneList` and a vector of ranks.
    ///
    /// # Errors
    /// - Returns `Err` if the lengths of `genes` and `ranks` do not match.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
    ///
    /// let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedGeneList::from(genes, ranks).unwrap();
    /// ```
    pub fn from(genes: GeneList, ranks: Vec<u32>) -> Result<Self, String> {
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

    /// Get the genes in this `RankedGeneList`.
    ///
    /// # Returns
    ///
    /// A slice of `Gene` instances.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
    ///
    /// let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedGeneList::from(genes, ranks).unwrap();
    ///
    /// let genes = ranked_list.genes();
    /// assert_eq!(genes.len(), 2);
    /// assert_eq!(genes[0].id(), "gene1");
    /// ```
    pub fn genes(&self) -> &GeneList {
        &self.genes
    }

    /// Get the ranks in this `RankedGeneList`.
    ///
    /// # Returns
    ///
    /// A slice of rank values (`u32`).
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
    ///
    /// let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedGeneList::from(genes, ranks).unwrap();
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
    /// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
    ///
    /// let ranks = vec![1, 2, 3];
    /// let genes = GeneList::from(vec![
    ///     Gene::from("gene1"),
    ///     Gene::from("gene2"),
    ///     Gene::from("gene3"),
    /// ]);
    /// let ranked_list = RankedGeneList::from(genes, ranks).unwrap();
    ///
    /// let thresholds = ranked_list.thresholds();
    /// assert_eq!(thresholds[0], 1);
    /// ```
    pub fn thresholds(&self) -> &[u32] {
        &self.thresholds
    }

    /// Retrieve the details of a ranked gene at a specific index.
    ///
    /// This method provides access to a `RankedGeneListItem` that includes a reference
    /// to the gene, its rank, and its threshold (if the list is thresholded). If the
    /// index is out of bounds, the method returns `None`.
    ///
    /// # Arguments
    /// - `index`: The index of the gene to retrieve.
    ///
    /// # Returns
    /// - `Some(RankedGeneListItem)` if the index is valid.
    /// - `None` if the index is out of bounds.
    ///
    /// # Notes
    /// - If the `RankedGeneList` is in an unthresholded state, the `threshold` field
    ///   in the returned `RankedGeneListItem` will be `None`.
    ///
    /// # Examples
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList, ThresholdState};
    ///
    /// let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();
    ///
    /// let item = ranked_list.get(1).unwrap();
    /// assert_eq!(item.gene().id(), "gene2");
    /// assert_eq!(item.rank(), 2);
    /// assert_eq!(item.threshold(), Some(2));
    ///
    /// // Simulate an unthresholded state
    /// ranked_list.remove(0);
    /// let item = ranked_list.get(0).unwrap();
    /// assert_eq!(item.threshold(), None);
    ///
    /// assert!(ranked_list.get(10).is_none()); // Out-of-bounds index
    /// ```
    pub fn get(&self, index: usize) -> Option<RankedGeneListItem> {
        if index >= self.genes.len() {
            return None;
        }

        Some(RankedGeneListItem {
            gene: &self.genes[index],
            rank: self.ranks[index],
            threshold: if self.threshold_state == ThresholdState::Thresholded {
                Some(self.thresholds[index])
            } else {
                None
            },
            index,
        })
    }

    /// Get the length of the `RankedGeneList`.
    /// 
    /// # Returns
    /// 
    /// The number of genes in the list.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
    /// 
    /// let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
    /// 
    /// let ranks = vec![1, 2];
    /// 
    /// let ranked_list = RankedGeneList::from(genes, ranks).unwrap();
    /// 
    /// assert_eq!(ranked_list.len(), 2);
    pub fn len(&self) -> usize {
        self.genes.len()
    }

    pub fn generate_thresholds(&mut self) {
        let mut thresholds = Vec::new();
        let mut current = 1;
        let max_rank = *self.ranks.last().unwrap_or(&0);

        while current <= max_rank {
            thresholds.push(current);
            current = ((current as f64) * 1.01 + 1.0).floor() as u32;
        }

        self.thresholds = thresholds;
    }

    /// Create an iterator for the `RankedGeneList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
    ///
    /// let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
    /// let ranks = vec![1, 2];
    /// let ranked_list = RankedGeneList::from(genes, ranks).unwrap();
    ///
    /// for item in ranked_list.iter() {
    ///     println!("Gene: {}, Rank: {}", item.gene().id(), item.rank());
    /// }
    /// ```
    pub fn iter(&self) -> RankedGeneListIterator {
        RankedGeneListIterator {
            ranked_gene_list: self,
            index: 0,
        }
    }

    /// Removes genes based on a single index or multiple indices.
    ///
    /// This method provides flexibility to remove either a single gene by its index
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
    /// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
    ///
    /// let genes = GeneList::from(vec![
    ///     Gene::from("gene1"),
    ///     Gene::from("gene2"),
    ///     Gene::from("gene3")
    /// ]);
    /// let ranks = vec![3, 2, 1];
    /// let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();
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

    /// Filters the `RankedGeneList` based on a user-provided closure and removes the matching elements.
    ///
    /// This method applies a filter function to each `RankedGeneListItem` in the list, collects
    /// the indices of the items for which the filter function returns `true`, and removes those items.
    ///
    /// # Arguments
    /// - `filter_fn`: A closure that takes a reference to a `RankedGeneListItem` and returns a `bool`.
    ///
    /// # Examples
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList, RankedGeneList};
    ///
    /// let genes = GeneList::from(vec![
    ///     Gene::from("gene1"),
    ///     Gene::from("gene2"),
    ///     Gene::from("gene3"),
    /// ]);
    /// let ranks = vec![3, 2, 1];
    /// let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();
    ///
    /// // Remove all items with a rank less than 3
    /// ranked_list.filter_and_remove(|item| item.rank() < 3);
    ///
    /// assert_eq!(ranked_list.len(), 1);
    /// assert_eq!(ranked_list.get(0).unwrap().gene().id(), "gene1");
    /// ```
    pub fn filter_and_remove<F>(&mut self, filter_fn: F)
    where
        F: Fn(&RankedGeneListItem) -> bool,
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

    /// Removes the gene and its corresponding rank at the specified index.
    ///
    /// This method removes a gene from the `RankedGeneList` at the given index, ensuring
    /// that the associated rank is also removed. If the gene list is in a thresholded state,
    /// it transitions to an unthresholded state upon successful removal.
    ///
    /// # Arguments
    ///
    /// - `index`: The index of the gene to remove.
    ///
    /// # Returns
    ///
    /// - `Ok(())` if the gene and rank are successfully removed.
    /// - `Err(String)` if the `index` is out of bounds or if the gene removal fails.
    ///
    /// # Errors
    ///
    /// - Returns an error if the provided `index` is greater than or equal to the length of the gene list.
    /// - Returns an error if the removal operation on the `GeneList` fails for any reason.
    fn remove_index(&mut self, index: usize) -> Result<(), String> {
        if index >= self.genes.len() {
            return Err(format!("Index {} is out of bounds", index));
        }
        self.genes.remove(index).map_err(|_| "Gene removal failed".to_string())?;
        self.ranks.remove(index);
        if self.threshold_state == ThresholdState::Thresholded {
            self.threshold_state = ThresholdState::Unthresholded;
        }
        Ok(())
    }

    fn check_lengths(genes: &GeneList, ranks: &[u32]) -> Result<(), String> {
        if genes.len() != ranks.len() {
            return Err(format!(
                "Gene list length ({}) does not match ranks length ({})",
                genes.len(),
                ranks.len()
            ));
        }
        Ok(())
    }

    fn sort_genes_and_ranks(&mut self) {
        let mut combined: Vec<(u32, &Gene)> = self
            .ranks
            .iter()
            .zip(self.genes.genes())
            .map(|(rank, gene)| (*rank, gene))
            .collect();
        combined.sort_by_key(|(rank, _)| *rank);
        self.ranks = combined.iter().map(|(rank, _)| *rank).collect();
        self.genes = GeneList::from(
            combined
                .into_iter()
                .map(|(_, gene)| gene.clone())
                .collect::<Vec<_>>(),
        );
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Test the iterator over a `RankedGeneList` with multiple elements.
    #[test]
    fn test_ranked_gene_list_iterator() {
        let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2"), Gene::from("gene3")]);
        let ranks = vec![3, 2, 1];
        let ranked_list = RankedGeneList::from(genes, ranks).unwrap();

        let mut iter = ranked_list.into_iter();
        let item = iter.next().unwrap();
        assert_eq!(item.gene().id(), "gene3");
        assert_eq!(item.rank(), 1);

        let item = iter.next().unwrap();
        assert_eq!(item.gene().id(), "gene2");
        assert_eq!(item.rank(), 2);

        let item = iter.next().unwrap();
        assert_eq!(item.gene().id(), "gene1");
        assert_eq!(item.rank(), 3);

        assert_eq!(iter.next(), None); // This will now work without error
    }


    /// Test that an empty `RankedGeneList` produces an empty iterator.
    #[test]
    fn test_ranked_gene_list_iterator_empty() {
        let ranked_list = RankedGeneList::new();
        let mut iter = ranked_list.iter();

        assert_eq!(iter.next(), None);
    }

    /// Test removal of a single valid index.
    #[test]
    fn test_remove_single_index() {
        let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2"), Gene::from("gene3")]);
        let ranks = vec![3, 2, 1];
        let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();

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
        let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2"), Gene::from("gene3")]);
        let ranks = vec![3, 2, 1];
        let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();

        let _ = ranked_list.remove(vec![0, 2]); // Remove "gene1" and "gene3"
        assert_eq!(ranked_list.genes().len(), 1);

        // Verify the remaining gene
        assert_eq!(ranked_list.genes()[0].id(), "gene2");
        assert_eq!(ranked_list.ranks(), &[2]);
    }

    /// Test removing indices out of bounds, which should panic.
    #[test]
    fn test_remove_invalid_index() {
        let genes = GeneList::from(vec![Gene::from("gene1")]);
        let ranks = vec![1];
        let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();

        let res = ranked_list.remove(5);

        match res {
            Err(e) => assert_eq!(e, "Index 5 is out of bounds"),
            _ => panic!("Expected an error"),
        }
    }

    /// Test that removal of all indices leaves an empty `RankedGeneList`.
    #[test]
    fn test_remove_all_indices() {
        let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2"), Gene::from("gene3")]);
        let ranks = vec![3, 2, 1];
        let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();

        let _ = ranked_list.remove(vec![0, 1, 2]); // Remove all indices
        assert!(ranked_list.genes().is_empty());
        assert!(ranked_list.ranks().is_empty());
    }

    /// Test removal of indices with repeated ranks to ensure stability.
    #[test]
    fn test_remove_indices_with_duplicate_ranks() {
        let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2"), Gene::from("gene3")]);
        let ranks = vec![2, 2, 1];
        let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();

        let _ = ranked_list.remove(vec![0, 2]); // Remove "gene1" and "gene3"
        assert_eq!(ranked_list.genes().len(), 1);

        // Verify the remaining gene
        assert_eq!(ranked_list.genes()[0].id(), "gene1");
        assert_eq!(ranked_list.ranks(), &[2]);
    }

    /// Test the iterator after removals to ensure consistency.
    #[test]
    fn test_iterator_after_removals() {
        let genes = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2"), Gene::from("gene3")]);
        let ranks = vec![3, 2, 1];
        let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();

        let _ = ranked_list.remove(vec![1]); // Remove "gene2"

        let mut iter = ranked_list.into_iter();

        // Check the first item
        let item = iter.next().unwrap();
        assert_eq!(item.gene().id(), "gene3");
        assert_eq!(item.rank(), 1);

        // Check the second item
        let item = iter.next().unwrap();
        assert_eq!(item.gene().id(), "gene1");
        assert_eq!(item.rank(), 3);

        // Ensure no more items
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn test_filter_and_remove() {
        let genes = GeneList::from(vec![
            Gene::from("gene1"),
            Gene::from("gene2"),
            Gene::from("gene3"),
            Gene::from("gene4"),
        ]);
        let ranks = vec![4, 3, 2, 1];
        let mut ranked_list = RankedGeneList::from(genes, ranks).unwrap();

        // Filter out all items with a rank greater than 2
        ranked_list.filter_and_remove(|item| item.rank() > 2);

        // Check the resulting list
        assert_eq!(ranked_list.len(), 2);
        let remaining_items: Vec<_> = ranked_list.iter().collect();
        assert_eq!(remaining_items[0].gene().id(), "gene4");
        assert_eq!(remaining_items[0].rank(), 1);
        assert_eq!(remaining_items[1].gene().id(), "gene3");
        assert_eq!(remaining_items[1].rank(), 2);
    }

    #[test]
    fn test_filter_by_another_ranked_gene_list() {
        let genes1 = GeneList::from(vec![
            Gene::from("gene1"),
            Gene::from("gene2"),
            Gene::from("gene3"),
            Gene::from("gene4"),
        ]);
        let ranks1 = vec![4, 3, 2, 1];
        let mut ranked_list1 = RankedGeneList::from(genes1, ranks1).unwrap();

        let genes2 = GeneList::from(vec![
            Gene::from("gene2"),
            Gene::from("gene4"),
        ]);
        let ranks2 = vec![1, 2];
        let ranked_list2 = RankedGeneList::from(genes2, ranks2).unwrap();

        // Use the second list to filter the first
        ranked_list1.filter_and_remove(|item| {
            !ranked_list2.iter().any(|other_item| other_item.gene().id() == item.gene().id())
        });

        // Check that only the genes from the second list remain in the first
        assert_eq!(ranked_list1.len(), 2);
        let remaining_items: Vec<_> = ranked_list1.iter().collect();
        assert_eq!(remaining_items[0].gene().id(), "gene4");
        assert_eq!(remaining_items[1].gene().id(), "gene2");

        // Verify ranks are intact
        assert_eq!(remaining_items[0].rank(), 1);
        assert_eq!(remaining_items[1].rank(), 3);
    }

}
