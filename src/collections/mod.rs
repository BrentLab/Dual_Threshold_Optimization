//! # Structs and utilities for storing Features (e.g. genes) and manipulating ranked gene sets.
//!
//! This module provides core data structures and utilities for representing and
//! manipulating ranked objects, ie a ranked gene set, for the dual threshold
//! optimization algorithm.
pub mod feature;
pub mod feature_list;
pub mod permuted;
pub mod ranked;
pub mod traits;
pub mod unique_check;

pub use feature::Feature;
pub use feature_list::FeatureList;
pub use permuted::PermutedRankedFeatureList;
pub use ranked::{RankedFeatureList, RankedFeatureListItem, RemoveIndices, ThresholdState};
pub use traits::FeatureSetProvider;
pub use unique_check::UniqueCheck;
