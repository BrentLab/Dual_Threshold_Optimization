//! # Structs and utilities for storing Features (e.g. genes) and manipulating ranked gene sets.
//!
//! This module provides core data structures and utilities for representing and
//! manipulating ranked objects, ie a ranked gene set, for the dual threshold
//! optimization algorithm.
pub mod feature;
pub mod feature_list;
pub mod permuted;
pub mod ranked;
pub mod unique_check;
pub mod traits;

pub use feature::Feature;
pub use feature_list::FeatureList;
pub use ranked::{RankedFeatureList, RankedFeatureListItem, RemoveIndices, ThresholdState};
pub use permuted::PermutedRankedFeatureList;
pub use unique_check::UniqueCheck;
pub use traits::FeatureSetProvider;