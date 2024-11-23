/// A struct representing a single gene with a unique identifier.
///
/// The `Gene` struct is designed to represent individual genes using a unique identifier,
/// such as a gene name or accession number. It ensures easy creation and retrieval of the identifier
/// while maintaining flexibility for efficient data handling.
///
/// # Fields
///
/// - `id`: A unique identifier for the gene (e.g., a gene name or accession number).
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::Gene;
///
/// let gene = Gene::from("gene1");
/// assert_eq!(gene.id(), "gene1");
/// ```
#[derive(Debug, PartialEq, Clone)]
pub struct Gene {
    /// The unique identifier for the gene.
    id: String,
}

impl Gene {
    /// Create a new `Gene` from a string slice.
    ///
    /// This method takes a string slice (`&str`) as input and creates a new `Gene`
    /// instance with the given identifier.
    ///
    /// # Arguments
    ///
    /// - `id`: A string slice representing the unique identifier for the gene.
    ///
    /// # Returns
    ///
    /// A new `Gene` instance.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::Gene;
    ///
    /// let gene = Gene::from("gene1");
    /// assert_eq!(gene.id(), "gene1");
    /// ```
    pub fn from(id: &str) -> Self {
        Self { id: id.to_string() }
    }

    /// Retrieve the identifier of the `Gene` as a string slice.
    ///
    /// This method provides an efficient way to access the gene's identifier without
    /// creating a new string, making it ideal for situations where the identifier
    /// needs to be read but not modified.
    ///
    /// # Returns
    ///
    /// A string slice representing the gene's identifier.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::Gene;
    ///
    /// let gene = Gene::from("gene1");
    /// assert_eq!(gene.id(), "gene1");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }
}
