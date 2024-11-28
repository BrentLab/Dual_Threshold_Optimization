//! # Core logic for the dual threshold optimization algorithm on two ranked lists.
//! 
//! The main function is `optimize()`. The asymptotic complexity of the algorithm
//! is $O(n^2)$, where $n$ is the length of the threshold lists. See 
//! `RankedFeatureList` for more information on the threshold lists.
pub mod compute_population_size;
pub mod optimize_main;
pub mod process_threshold_pairs;
pub mod results_objects;

pub use compute_population_size::compute_population_size;
pub use optimize_main::optimize;
pub use process_threshold_pairs::process_threshold_pairs;
pub use results_objects::{FeatureSets, OptimizationResult, OptimizationResultRecord};