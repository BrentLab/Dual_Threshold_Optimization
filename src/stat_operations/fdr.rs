/// Calculate the expected false discovery rate (E\[FDR\])
///
/// # Arguments
///
/// * `ranked_list_1_len` - The length of the first ranked list (genes bound by the TF in the assay)
/// * `ranked_list_2_len` - The length of the second ranked list (genes responsive when the TF is perturbed)
/// * `overlap_len` - The size of the overlap between the two ranked lists (R ∩ B)
/// * `population_size` - The total population size (all genes assayed in both binding and response experiments)
/// * `sensitivity` - Sensitivity of the intersection algorithm (default is 0.8)
///
/// # Returns
///
/// The expected false discovery rate (E\[FDR\]) as a `f64`.
///
/// # Examples
///
/// ```
/// use dual_threshold_optimization::stat_operations::fdr;
///
/// let ranked_list_1_len = 200; // Size of set B
/// let ranked_list_2_len = 150; // Size of set R
/// let overlap_len = 50; // Size of R ∩ B
/// let population_size = 1000; // Size of G
/// let sensitivity = 0.8; // Default sensitivity
///
/// let est_fdr = fdr(ranked_list_1_len, ranked_list_2_len, overlap_len, population_size, sensitivity);
/// assert!((est_fdr - 0.240625).abs() < 1e-4); // Expected value is approximately 0.240625
/// ```
pub fn fdr(
    ranked_list_1_len: usize,
    ranked_list_2_len: usize,
    overlap_len: usize,
    population_size: u64,
    sensitivity: f64,
) -> f64 {
    // Sensitivity cannot be zero to avoid division by zero
    if sensitivity <= 0.0 {
        panic!("Sensitivity must be greater than 0.");
    }

    // Calculate |DF| estimate using sensitivity (Sn)
    let df = (overlap_len as f64 / sensitivity).max(0.0);

    // Calculate max(0, |B| - |DF|) and max(0, |R| - |DF|)
    let b_minus_df = (ranked_list_1_len as f64 - df).max(0.0);
    let r_minus_df = (ranked_list_2_len as f64 - df).max(0.0);

    // Calculate the numerator: max(0, |B| - |R ∩ B| / Sn) * max(0, |R| - |R ∩ B| / Sn)
    let numerator = b_minus_df * r_minus_df;

    // Calculate the denominator: |G| * |B ∩ R|
    let denominator = (population_size as f64) * (overlap_len as f64);

    // Return the expected false discovery rate
    if denominator > 0.0 {
        numerator / denominator
    } else {
        0.0 // Avoid division by zero
    }
}
