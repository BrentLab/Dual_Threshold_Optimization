//! # Functions for execution over multiple permutations and threads.
//!
//! This offers both a method to parallelize on a single machine or across multiple
//! machines with MPI.
pub mod multi_node;
pub mod single_node;
pub mod task;

pub use multi_node::run as run_multi_node;
pub use single_node::run as run_single_node;
pub use task::Task;
