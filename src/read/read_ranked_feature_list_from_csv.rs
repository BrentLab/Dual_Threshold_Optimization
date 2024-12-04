//! Read in a CSV file containing features (e.g. genes) and their ranks.
//!
//! This module provides a function to read in a CSV file containing features
//! and their ranks, returning a `RankedFeatureList` object. The file is expected to
//! be a CSV file with two columns: `gene` and `rank` **without a header**. The ranks
//! are expected to be integers where ties have been handled with either a `min` or
//! `max` method where all features with with identical underlying scores are assigned
//! the same integer rank.
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::vec::Vec;

use crate::collections::{Feature, FeatureList, RankedFeatureList};

/// Reads a CSV file containing genes and their ranks, returning a `RankedFeatureList`.
///
/// # Arguments
/// - `filepath`: Path to the CSV file. The file must have two columns: `gene` and `rank`.
///
/// # Returns
/// A `RankedFeatureList` containing the genes and their ranks.
///
/// # Panics
/// Panics if the file cannot be opened, has invalid formatting, or if creating the
/// `RankedFeatureList` fails.
///
/// # Examples
///
/// Writing to a temporary file and reading it back:
/// ```rust
/// use std::fs::File;
/// use std::io::Write;
/// use std::env::temp_dir;
/// use dual_threshold_optimization::read::read_ranked_feature_list_from_csv;
/// use dual_threshold_optimization::collections::RankedFeatureList;
///
/// let temp_path = temp_dir().join("temp_genes.csv");
/// let temp_file = temp_path.to_str().unwrap();
///
/// // Write test data to the temporary file
/// let mut file = File::create(temp_file).unwrap();
/// writeln!(file, "gene1,1").unwrap();
/// writeln!(file, "gene2,2").unwrap();
/// writeln!(file, "gene3,3").unwrap();
///
/// // Read and test the RankedFeatureList
/// let ranked_feature_list = read_ranked_feature_list_from_csv(temp_file);
/// assert_eq!(ranked_feature_list.len(), 3);
/// ```
pub fn read_ranked_feature_list_from_csv(filepath: &str) -> RankedFeatureList {
    let file = File::open(filepath).expect("Could not open file");
    let reader = BufReader::new(file);
    let mut genes = Vec::new();
    let mut ranks = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.expect("Could not read line");
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() != 2 {
            panic!("Invalid format in file {} at line {}", filepath, i + 1);
        }
        let feature = Feature::from(fields[0].trim());
        let rank: usize = fields[1].trim().parse().expect("Invalid rank value");
        genes.push(feature);
        ranks.push(rank as u32);
    }

    RankedFeatureList::from(FeatureList::from(genes), ranks)
        .expect("Failed to create RankedFeatureList")
}
