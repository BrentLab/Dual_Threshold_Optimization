use std::collections::HashMap;
use std::collections::HashSet;
use std::ops::Index;

use crate::{Gene, UniqueCheck};

/// A collection of `Gene` objects, with methods for manipulation and validation.
///
/// The `GeneList` struct provides methods for managing and validating a collection of genes.
/// It ensures that operations like adding, removing, or checking for unique identifiers are simple and efficient.
///
/// # Fields
///
/// - `genes`: A vector of `Gene` instances.
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::{Gene, GeneList};
///
/// let genes = vec![Gene::from("gene1"), Gene::from("gene2")];
/// let gene_list = GeneList::from(genes);
/// assert_eq!(gene_list.len(), 2);
/// ```
#[derive(Debug, PartialEq)]
pub struct GeneList {
    genes: Vec<Gene>,
}

impl GeneList {
    /// Create an empty `GeneList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::GeneList;
    ///
    /// let gene_list = GeneList::new();
    /// assert_eq!(gene_list.len(), 0);
    /// ```
    pub fn new() -> Self {
        Self { genes: Vec::new() }
    }

    /// Create a `GeneList` from an initial vector of genes.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    ///
    /// let genes = vec![Gene::from("gene1"), Gene::from("gene2")];
    /// let gene_list = GeneList::from(genes);
    /// assert_eq!(gene_list.len(), 2);
    /// ```
    pub fn from(genes: Vec<Gene>) -> Self {
        Self { genes }
    }

    /// Add a `Gene` to the list.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    ///
    /// let mut gene_list = GeneList::new();
    /// gene_list.push(Gene::from("gene1"));
    /// assert_eq!(gene_list.len(), 1);
    /// ```
    pub fn push(&mut self, gene: Gene) {
        self.genes.push(gene);
    }

    /// Remove the last `Gene` from the list and return it.
    /// Returns `None` if the list is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    ///
    /// let mut gene_list = GeneList::new();
    /// gene_list.push(Gene::from("gene1"));
    /// let popped_gene = gene_list.pop();
    /// assert_eq!(popped_gene.unwrap().id(), "gene1");
    /// assert!(gene_list.is_empty());
    /// ```
    pub fn pop(&mut self) -> Option<Gene> {
        self.genes.pop()
    }

    /// Check if all `Gene` IDs in the list are unique.
    ///
    /// - Returns `UniqueCheck::Unique` if all IDs are unique.
    /// - Returns `UniqueCheck::Duplicates` with a map of duplicate IDs and their indices if duplicates are found.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList, UniqueCheck};
    ///
    /// let genes = vec![Gene::from("gene1"), Gene::from("gene2"), Gene::from("gene1")];
    /// let gene_list = GeneList::from(genes);
    ///
    /// match gene_list.check_unique() {
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
        for (index, gene) in self.genes.iter().enumerate() {
            seen.entry(gene.id()).or_default().push(index);
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
    /// A reference to the vector of `Gene` instances stored in the `GeneList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    ///
    /// let genes = vec![Gene::from("gene1"), Gene::from("gene2")];
    /// let gene_list = GeneList::from(genes);
    ///
    /// let internal_genes = gene_list.genes();
    /// assert_eq!(internal_genes.len(), 2);
    /// assert_eq!(internal_genes[0].id(), "gene1");
    /// ```
    pub fn genes(&self) -> &[Gene] {
        &self.genes
    }

    /// Get the number of genes in the list.
    ///
    /// # Returns
    ///
    /// The total number of `Gene` instances in the `GeneList`.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    ///
    /// let genes = vec![Gene::from("gene1"), Gene::from("gene2")];
    /// let gene_list = GeneList::from(genes);
    ///
    /// assert_eq!(gene_list.len(), 2);
    /// ```
    pub fn len(&self) -> usize {
        self.genes.len()
    }

    /// Check if the gene list is empty.
    ///
    /// # Returns
    ///
    /// `true` if the `GeneList` contains no `Gene` instances, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    ///
    /// let empty_gene_list = GeneList::new();
    /// assert!(empty_gene_list.is_empty());
    ///
    /// let genes = vec![Gene::from("gene1")];
    /// let non_empty_gene_list = GeneList::from(genes);
    /// assert!(!non_empty_gene_list.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.genes.is_empty()
    }

    /// Remove a `Gene` from the list by index.
    /// Returns the removed `Gene` if successful,
    /// or an error message if the index is out of bounds.
    /// 
    /// # Arguments
    /// 
    /// * `index` - The index of the `Gene` to remove.
    /// 
    /// # Returns
    /// 
    /// - `Ok(Gene)` if the `Gene` was successfully removed.
    /// - `Err(String)` if the index is out of bounds.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    /// 
    /// let mut gene_list = GeneList::from(vec![
    ///    Gene::from("gene1"),
    ///   Gene::from("gene2"),
    ///   Gene::from("gene3"),
    /// ]);
    /// 
    /// let removed_gene = gene_list.remove(1);
    /// match removed_gene {
    ///     Ok(gene) => println!("Removed gene: {}", gene.id()),
    ///     Err(msg) => println!("Error: {}", msg.to_string()),
    /// }
    /// assert_eq!(gene_list.len(), 2);
    /// assert_eq!(gene_list.genes()[0].id(), "gene1");
    /// assert_eq!(gene_list.genes()[1].id(), "gene3");
    /// ```
    pub fn remove(&mut self, index: usize) -> Result<Gene, String> {
    if index >= self.genes.len() {
        return Err(format!(
            "Index out of bounds: {} (length: {})",
            index,
            self.genes.len()
        ));
    }
    Ok(self.genes.remove(index))
    }

    /// Computes the intersection of two `GeneList`s.
    ///
    /// This method identifies genes that are present in both `GeneList`s
    /// based on their IDs and returns a new `GeneList` containing the common genes.
    ///
    /// # Arguments
    /// - `other`: The other `GeneList` to intersect with.
    ///
    /// # Returns
    /// A new `GeneList` containing the intersection of the two `GeneList`s.
    ///
    /// # Examples
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    ///
    /// let list1 = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
    /// let list2 = GeneList::from(vec![Gene::from("gene2"), Gene::from("gene3")]);
    ///
    /// let intersection = list1.intersect(&list2);
    ///
    /// assert_eq!(intersection.len(), 1);
    /// assert_eq!(intersection[0].id(), "gene2");
    /// ```
    pub fn intersect<'a>(&'a self, other: &GeneList) -> Vec<&'a Gene> {
        let ids1: HashSet<_> = self.genes().iter().map(|g| g.id()).collect();
        let ids2: HashSet<_> = other.genes().iter().map(|g| g.id()).collect();

        self
            .genes()
            .iter()
            .filter(|gene| ids1.contains(gene.id()) && ids2.contains(gene.id()))
            .collect()
    }


    /// Computes the genes that are either exclusive to one of the two `GeneList`s
    /// (symmetric difference) or only in `self` but not in `other` (set difference).
    ///
    /// # Arguments
    /// - `other`: The other `GeneList` to compare against.
    /// - `symmetric`: If `true`, computes the symmetric difference (genes in one list but not both).
    ///                If `false`, computes the set difference (genes in `self` but not in `other`).
    ///
    /// # Returns
    /// A new `GeneList` containing the resulting genes.
    ///
    /// # Examples
    /// ## Symmetric Difference:
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    ///
    /// let list1 = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
    /// let list2 = GeneList::from(vec![Gene::from("gene2"), Gene::from("gene3")]);
    ///
    /// let exclusive = list1.difference(&list2, true);
    ///
    /// assert_eq!(exclusive.len(), 2);
    /// assert!(exclusive.iter().any(|g| g.id() == "gene1"));
    /// assert!(exclusive.iter().any(|g| g.id() == "gene3"));
    /// ```
    /// ## Set Difference:
    /// ```
    /// use dual_threshold_optimization::{Gene, GeneList};
    ///
    /// let list1 = GeneList::from(vec![Gene::from("gene1"), Gene::from("gene2")]);
    /// let list2 = GeneList::from(vec![Gene::from("gene2"), Gene::from("gene3")]);
    ///
    /// let exclusive = list1.difference(&list2, false);
    ///
    /// assert_eq!(exclusive.len(), 1);
    /// assert_eq!(exclusive[0].id(), "gene1");
    /// ```
    pub fn difference<'a>(
        &'a self,
        other: &'a GeneList,
        symmetric: bool,
    ) -> Vec<&'a Gene> {
        let ids1: HashSet<_> = self.genes().iter().map(|g| g.id()).collect();
        let ids2: HashSet<_> = other.genes().iter().map(|g| g.id()).collect();

        if symmetric {
            self.genes()
                .iter()
                .chain(other.genes().iter())
                .filter(|gene| ids1.symmetric_difference(&ids2).any(|&id| id == gene.id()))
                .collect()
        } else {
            self.genes()
                .iter()
                .filter(|gene| ids1.difference(&ids2).any(|&id| id == gene.id()))
                .collect()
        }
    }
}

impl From<Vec<String>> for GeneList {
    fn from(ids: Vec<String>) -> Self {
        let genes: Vec<Gene> = ids.into_iter().map(|id| Gene::from(id.as_str())).collect();
        GeneList { genes }
    }
}

impl From<Vec<&str>> for GeneList {
    fn from(ids: Vec<&str>) -> Self {
        let genes: Vec<Gene> = ids.into_iter().map(Gene::from).collect();
        GeneList { genes }
    }
}

impl Index<usize> for GeneList {
    type Output = Gene;

    fn index(&self, index: usize) -> &Self::Output {
        &self.genes[index]
    }
}