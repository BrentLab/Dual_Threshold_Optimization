pub mod gene;
pub mod gene_list;
pub mod hypergeometric_pvalue;
pub mod intersect_genes;
pub mod optimize;
pub mod ranked_gene_list;
pub mod unique_check;

// Re-export for easier access from the core module
pub use gene::Gene;
pub use gene_list::GeneList;
pub use hypergeometric_pvalue::hypergeometric_pvalue;
pub use intersect_genes::intersect_genes;
pub use optimize::*;
pub use ranked_gene_list::{RankedGeneList, RankedGeneListItem, RemoveIndices, ThresholdState};
pub use unique_check::UniqueCheck;
