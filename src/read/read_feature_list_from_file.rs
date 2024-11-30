//! A utility to read in a list of features (e.g. genes) from a file
//!
//! The file is expected to be a single column with no header where each row
//! contains a single feature identifier. This is used to read in the optional
//! background set.
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::vec::Vec;

use crate::collections::{Feature, FeatureList};

/// Reads a file containing a list of genes (one per line) to create a `FeatureList`.
///
/// # Arguments
/// - `filepath`: Path to the file. Each line must contain a single feature identifier.
///
/// # Returns
/// A `FeatureList` containing all the genes from the file.
///
/// # Panics
/// Panics if the file cannot be opened or if there are errors reading the lines.
///
/// # Examples
///
/// Writing to a temporary file and reading it back:
/// ```
/// use std::fs::File;
/// use std::io::Write;
/// use std::env::temp_dir;
/// use dual_threshold_optimization::collections::FeatureList;
/// use dual_threshold_optimization::read::read_feature_list_from_file;
///
/// let temp_path = temp_dir().join("temp_genes.csv");
/// let temp_file = temp_path.to_str().unwrap();
///
/// // Write test data to the temporary file
/// let mut file = File::create(temp_file).unwrap();
/// writeln!(file, "gene1").unwrap();
/// writeln!(file, "gene2").unwrap();
/// writeln!(file, "gene3").unwrap();
///
/// let feature_list = read_feature_list_from_file(temp_file);
/// assert_eq!(feature_list.len(), 3);
/// ```
pub fn read_feature_list_from_file(filepath: &str) -> FeatureList {
    let file = File::open(filepath).expect("Could not open file");
    let reader = BufReader::new(file);
    let genes: Vec<Feature> = reader
        .lines()
        .map(|line| Feature::from(line.expect("Could not read line").trim()))
        .collect();
    FeatureList::from(genes)
}
