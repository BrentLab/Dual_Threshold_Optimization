//! # Statistical operations
//!
//! This module contains functions to calculate p-values and other statistical
//! measures.
pub mod empirical_pvalue;
pub mod fdr;
pub mod hypergeometric_pvalue;
pub mod intersect_genes;

pub use empirical_pvalue::empirical_pvalue;
pub use fdr::fdr;
pub use hypergeometric_pvalue::hypergeometric_pvalue;
pub use intersect_genes::intersect_genes;
