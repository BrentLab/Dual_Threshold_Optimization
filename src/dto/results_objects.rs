//! # Objects to represent the results of the optimization process.
use serde::{Deserialize, Serialize};

use crate::collections::Feature;

/// An enum to represent the gene1, gene2 sets output by optimize().
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub enum FeatureSets {
    Both(Vec<Feature>, Vec<Feature>),
    None,
}

/// An enum to represent the result of optimize().
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum OptimizationResult {
    Debug(Vec<OptimizationResultRecord>), // All results
    Best(OptimizationResultRecord),       // Single best result
}

/// A struct to record the result of optimize().
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct OptimizationResultRecord {
    pub rank1: usize,
    pub rank2: usize,
    pub set1_len: usize,
    pub set2_len: usize,
    pub population_size: u64,
    pub intersection_size: usize,
    pub pvalue: f64,
    pub permuted: bool,
    pub feature_sets: FeatureSets,
}