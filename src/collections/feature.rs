//! # Represents a single feature (e.g. gene)
//! 
//! This struct is used as the type which is stored in a vector in `FeatureList`
use serde::{Deserialize, Serialize};

/// A struct representing a single feature with a unique identifier.
///
/// The `Feature` struct is designed to represent individual genes using a unique identifier,
/// such as a feature name or accession number. It ensures easy creation and retrieval of the identifier
/// while maintaining flexibility for efficient data handling.
///
/// # Fields
///
/// - `id`: A unique identifier for the feature (e.g., a feature name or accession number).
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::collections::Feature;
///
/// let feature = Feature::from("gene1");
/// assert_eq!(feature.id(), "gene1");
/// ```
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Feature {
    /// The unique identifier for the feature.
    id: String,
}

impl Feature {
    /// Create a new `Feature` from a string slice.
    ///
    /// This method takes a string slice (`&str`) as input and creates a new `Feature`
    /// instance with the given identifier.
    ///
    /// # Arguments
    ///
    /// - `id`: A string slice representing the unique identifier for the feature.
    ///
    /// # Returns
    ///
    /// A new `Feature` instance.
    ///
    /// # Examples
    ///
    /// ```
    /// use dual_threshold_optimization::collections::Feature;
    ///
    /// let feature = Feature::from("gene1");
    /// assert_eq!(feature.id(), "gene1");
    /// ```
    pub fn from(id: &str) -> Self {
        Self { id: id.to_string() }
    }

    /// Retrieve the identifier of the `Feature` as a string slice.
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
    /// use dual_threshold_optimization::collections::Feature;
    ///
    /// let feature = Feature::from("gene1");
    /// assert_eq!(feature.id(), "gene1");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }
}
