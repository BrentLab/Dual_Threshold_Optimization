#![doc = include_str!("../README.md")]
pub mod collections;
pub mod dto;
pub mod read;
pub mod run;
pub mod stat_operations;

// Re-export for easier access from the core module
pub use collections::{PermutedRankedFeatureList, RankedFeatureList};
pub use dto::optimize;
