//! # Execute a set of optimization tasks in parallel on a single shared memory device.
//!
//! This function is used to execute a set of optimization tasks in parallel on a single
//! shared memory device. A `Task` is defined as a call to `optimize` and a set of
//! tasks is expected to be a list of at least length 1 where at least one task is
//! the unpermuted optimization and the rest are permutations which will be used to
//! calculate the empirical p-value of the unpermuted optimization.
use std::clone::Clone;
use std::sync::{Arc, Mutex};

use crate::collections::RankedFeatureList;
use crate::dto::OptimizationResult;
use crate::optimize;
use crate::run::Task;

/// Executes multiple tasks for dual threshold optimization in parallel.
///
/// This function distributes tasks across a specified number of threads, leveraging
/// multithreading for efficient parallel execution of the `optimize`
/// function. Each thread processes a chunk of tasks, and the results are aggregated
/// and returned as a single vector.
///
/// # Arguments
/// - `tasks`: A vector of `Task` objects defining the parameters for each optimization.
/// - `ranked_feature_list1`: The first `RankedFeatureList` for the optimization.
/// - `ranked_feature_list2`: The second `RankedFeatureList` for the optimization.
/// - `background`: An optional `FeatureList` defining the background population. If absent,
///   the two ranked lists must have identical genes.
/// - `num_threads`: The number of threads to use for parallel processing.
///
/// # Returns
/// A `Vec<OptimizationResult>` containing the results for all tasks.
///
/// # Panics
/// - If the number of tasks is smaller than the number of threads.
/// - If a thread fails to complete execution.
///
/// # Implementation Details
/// - The `tasks` vector is divided into chunks, and each thread processes one chunk.
/// - The `ranked_feature_list1`, `ranked_feature_list2`, and `background` are cloned for each thread.
/// - Results are stored in a thread-safe `Arc<Mutex<Vec<OptimizationResult>>>` to avoid race conditions.
///
/// # Examples
/// ```rust
/// use dual_threshold_optimization::run::{run_single_node, Task};
/// use dual_threshold_optimization::collections::{FeatureList, RankedFeatureList, Feature};
/// use dual_threshold_optimization::dto::OptimizationResult;
///
/// // Define RankedFeatureLists and background
/// let ranked_feature_list1 = RankedFeatureList::from(
///     FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2"), Feature::from("gene3")]),
///     vec![1, 2, 3],
/// ).unwrap();
///
/// let ranked_feature_list2 = RankedFeatureList::from(
///     FeatureList::from(vec![Feature::from("gene1"), Feature::from("gene2"), Feature::from("gene4")]),
///     vec![1, 2, 3],
/// ).unwrap();
///
/// let population_size: usize = 4;
///
/// // Define tasks
/// let tasks: Vec<Task> = vec![
///     Task { id: 1, permute: true }, // Replace with actual task definitions
///     Task { id: 2, permute: true },
/// ];
///
/// // Run with 2 threads
/// let results = run_single_node(tasks, ranked_feature_list1, ranked_feature_list2, population_size, 2);
///
/// // Check results
/// assert_eq!(results.len(), 2);
/// for result in results {
///     match result {
///         OptimizationResult::Debug(_) => println!("Debug results returned."),
///         OptimizationResult::Best(best) => println!(
///             "Best result: Rank1 {}, Rank2 {}, p-value {}",
///             best.rank1, best.rank2, best.pvalue
///         ),
///     }
/// }
/// ```
pub fn run(
    tasks: Vec<Task>,
    ranked_feature_list1: RankedFeatureList,
    ranked_feature_list2: RankedFeatureList,
    population_size: usize,
    num_threads: usize,
) -> Vec<OptimizationResult> {
    // Wrap ranked feature lists and background in Arc for thread-safe shared access
    let ranked_feature_list1 = Arc::new(ranked_feature_list1);
    let ranked_feature_list2 = Arc::new(ranked_feature_list2);

    // Split the tasks into chunks
    let task_chunks: Vec<Vec<Task>> = tasks
        .chunks(tasks.len().div_ceil(num_threads)) // Avoid zero chunks if num_threads > tasks.len()
        .map(|chunk| chunk.to_vec())
        .collect();

    // Shared results vector protected by a mutex
    let results = Arc::new(Mutex::new(Vec::new()));

    // Spawn threads
    let handles: Vec<_> = task_chunks
        .into_iter()
        .map(|chunk| {
            // Clone the Arc-wrapped data for each thread
            let ranked_feature_list1 = Arc::clone(&ranked_feature_list1);
            let ranked_feature_list2 = Arc::clone(&ranked_feature_list2);
            let results = Arc::clone(&results);

            std::thread::spawn(move || {
                let mut local_results = Vec::new();
                for task in chunk {
                    let result = optimize(
                        &ranked_feature_list1,
                        &ranked_feature_list2,
                        task.permute,
                        population_size,
                        false,
                    );
                    local_results.push(result);
                }
                // Lock and extend the shared results vector
                results.lock().unwrap().extend(local_results);
            })
        })
        .collect();

    // Wait for all threads to complete
    for handle in handles {
        handle.join().unwrap();
    }

    // Unwrap the Arc and Mutex to get the final results
    Arc::try_unwrap(results).unwrap().into_inner().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    };
    use std::thread;

    use crate::collections::{Feature, FeatureList};

    #[test]
    fn test_parallel_execution_with_multiple_threads() {
        // Create dummy RankedFeatureLists
        let ranked_feature_list1 = RankedFeatureList::from(
            FeatureList::from(vec![
                Feature::from("gene1"),
                Feature::from("gene2"),
                Feature::from("gene3"),
            ]),
            vec![1, 2, 3],
        )
        .unwrap();

        let ranked_feature_list2 = RankedFeatureList::from(
            FeatureList::from(vec![
                Feature::from("gene1"),
                Feature::from("gene2"),
                Feature::from("gene4"),
            ]),
            vec![1, 2, 3],
        )
        .unwrap();

        let population_size: usize = 4;

        // Create dummy tasks
        let tasks: Vec<Task> = (1..=10).map(|id| Task { id, permute: true }).collect();

        // Shared counter to track active threads
        let active_threads = Arc::new(AtomicUsize::new(0));
        let max_active_threads = Arc::new(Mutex::new(0));

        let task_chunks: Vec<Vec<Task>> = tasks
            .chunks(tasks.len() / 5) // 5 threads
            .map(|chunk| chunk.to_vec())
            .collect();

        let handles: Vec<_> = task_chunks
            .into_iter()
            .map(|chunk| {
                let ranked_feature_list1 = ranked_feature_list1.clone();
                let ranked_feature_list2 = ranked_feature_list2.clone();
                let active_threads = Arc::clone(&active_threads);
                let max_active_threads = Arc::clone(&max_active_threads);

                thread::spawn(move || {
                    active_threads.fetch_add(1, Ordering::SeqCst);
                    let current_active = active_threads.load(Ordering::SeqCst);

                    // Update max active threads
                    let mut max_active = max_active_threads.lock().unwrap();
                    if current_active > *max_active {
                        *max_active = current_active;
                    }

                    // Simulate work
                    for _ in chunk {
                        optimize(
                            &ranked_feature_list1,
                            &ranked_feature_list2,
                            true,
                            population_size,
                            false,
                        );
                    }

                    active_threads.fetch_sub(1, Ordering::SeqCst);
                })
            })
            .collect();

        for handle in handles {
            handle.join().unwrap();
        }

        // Assert that at some point we had multiple active threads
        let max_threads = *max_active_threads.lock().unwrap();
        assert!(
            max_threads > 1,
            "Expected multiple threads to run simultaneously, but found only one active thread."
        );
    }
}

// The following was a previous implementation of the run_single_node function
// which did not use the Arc<Mutex<Vec<OptimizationResult>>> for thread-safe
// access to the rankedlists and background. kept for comparison purposes
// pub fn run(
//     tasks: Vec<Task>,
//     ranked_feature_list1: RankedFeatureList,
//     ranked_feature_list2: RankedFeatureList,
//     background: Option<FeatureList>,
//     num_threads: usize,
// ) -> Vec<OptimizationResult> {
//     // Split the tasks into chunks and ensure ownership is transferred
//     let task_chunks: Vec<Vec<Task>> = tasks
//         .chunks(tasks.len() / num_threads)
//         .map(|chunk| chunk.to_vec())
//         .collect();

//     let results = Arc::new(Mutex::new(Vec::new()));

//     let handles: Vec<_> = task_chunks
//         .into_iter()
//         .map(|chunk| {
//             let ranked_feature_list1 = ranked_feature_list1.clone();
//             let ranked_feature_list2 = ranked_feature_list2.clone();
//             let background = background.clone();
//             let results = Arc::clone(&results);

//             std::thread::spawn(move || {
//                 let mut local_results = Vec::new();
//                 for _ in chunk {
//                     let result = optimize(
//                         &ranked_feature_list1,
//                         &ranked_feature_list2,
//                         true,
//                         background.as_ref(),
//                         false,
//                     );
//                     local_results.push(result);
//                 }
//                 results.lock().unwrap().extend(local_results);
//             })
//         })
//         .collect();

//     for handle in handles {
//         handle.join().unwrap();
//     }

//     Arc::try_unwrap(results).unwrap().into_inner().unwrap()
// }
