//! Calculate the hypergeometric p-value for two sets
use statrs::distribution::{Hypergeometric, DiscreteCDF};

/// Calculate the hypergeometric p-value.
///
/// This function computes the probability of observing an intersection size at least
/// as large as the observed size (`k`) between two sets, assuming random selection.
///
/// # Arguments
///
/// * `population_size` (`usize`) - The total number of items in the population (`N`).
/// * `successes_in_population` (`usize`) - The size of one set (`K`).
/// * `sample_size` (`usize`) - The size of the other set (`n`).
/// * `observed_overlap` (`usize`) - The size of the intersection between the
///    two sets (`k`).
///
/// # Returns
///
/// The p-value representing the probability of observing an overlap size at least as large
/// as `observed_overlap` under the hypergeometric distribution.
///
/// # Example
///
/// ```
/// use dual_threshold_optimization::stat_operations::hypergeometric_pvalue;
///
/// let p_value = hypergeometric_pvalue(1000, 50, 60, 10);
/// assert_eq!(p_value, 0.00044068070222441115);
/// 
/// let p_value = hypergeometric_pvalue(6060, 5808, 154, 153);
/// assert_eq!(p_value, 0.010413637619010246);
/// ```
pub fn hypergeometric_pvalue(
    population_size: u64,
    successes_in_population: u64,
    sample_size: u64,
    observed_overlap: u64,
) -> f64 {
    // Create the hypergeometric distribution
    let hypergeom = Hypergeometric::new(population_size, successes_in_population, sample_size)
        .expect("Failed to create hypergeometric distribution");

        // Handle the edge case where observed_overlap is 0
    if observed_overlap == 0 {
        return 1.0; // P(X >= 0) is always 1
    }

    // Use the survival function (SF) for the upper tail probability
    hypergeom.sf(observed_overlap - 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hypergeometric_pvalue_scenario_1() {
        // Test case 1: i: 0, j: 0, population_size: 3, successes_in_pop: 1, sample_size: 1, intersect: 1
        let p_value = hypergeometric_pvalue(3, 1, 1, 1);
        assert!((p_value - 0.3333333333333333).abs() < 1e-10);

        // Test case 2: i: 0, j: 1, population_size: 3, successes_in_pop: 1, sample_size: 2, intersect: 1
        let p_value = hypergeometric_pvalue(3, 1, 2, 1);
        assert!((p_value - 0.6666666666666666).abs() < 1e-10);

        // Test case 3: i: 0, j: 2, population_size: 3, successes_in_pop: 1, sample_size: 3, intersect: 1
        let p_value = hypergeometric_pvalue(3, 1, 3, 1);
        assert!((p_value - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hypergeometric_pvalue_scenario_2() {
        // Test case 4: i: 1, j: 0, population_size: 3, successes_in_pop: 2, sample_size: 1, intersect: 1
        let p_value = hypergeometric_pvalue(3, 2, 1, 1);
        assert!((p_value - 0.6666666666666666).abs() < 1e-10);

        // Test case 5: i: 1, j: 1, population_size: 3, successes_in_pop: 2, sample_size: 2, intersect: 2
        let p_value = hypergeometric_pvalue(3, 2, 2, 2);
        assert!((p_value - 0.3333333333333333).abs() < 1e-10);
    }

    #[test]
    fn test_hypergeometric_pvalue_scenario_3() {
        // Test case 4: i: 1, j: 0, population_size: 3, successes_in_pop: 2, sample_size: 1, intersect: 1
        let p_value = hypergeometric_pvalue(10, 1, 1, 0);
        assert!((p_value - 1.0).abs() < 1e-10);
    }
}
