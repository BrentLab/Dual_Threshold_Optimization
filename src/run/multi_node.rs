//! # Execute a set of optimization tasks using MPI for distributed processing.
//! 
//! This function is used to execute a set of optimization tasks in parallel by first
//! dividing the tasks into roughly equal parts and broadcasting them to the MPI 
//! processes. On those MPI processes, if multiple threads are available, those tasks
//! can be further parallelized. A `Task` is defined as a call to `optimize` and a set of
//! tasks is expected to be a list of at least length 1 where at least one task is
//! the unpermuted optimization and the rest are permutations which will be used to
//! calculate the empirical p-value of the unpermuted optimization.
#[cfg(feature = "mpi")]
use bincode;
#[cfg(feature = "mpi")]
use mpi::traits::*;
#[cfg(feature = "mpi")]
use mpi::topology::Communicator;
#[cfg(feature = "mpi")]
use crate::read::read_ranked_feature_list_from_csv;

#[cfg(feature = "mpi")]
use crate::run::run_single_node;

use crate::dto::OptimizationResult;
use crate::run::Task;

/// Executes tasks for dual threshold optimization using MPI for distributed processing.
///
/// This function uses MPI to distribute tasks across multiple processes. The root process
/// divides the tasks into chunks and sends them to other processes. Each process computes
/// its assigned tasks and returns the results to the root. The root process aggregates
/// the results and returns them.
///
/// # Arguments
/// - `tasks`: A vector of `Task` objects defining the parameters for each optimization.
/// - `ranked_feature_list1`: The first `RankedFeatureList` for the optimization.
/// - `ranked_feature_list2`: The second `RankedFeatureList` for the optimization.
/// - `background`: An optional `FeatureList` defining the background population.
///
/// # Returns
/// A `Vec<OptimizationResult>` containing the results for all tasks. Only the root process
/// will have a non-empty vector; worker processes return an empty vector.
///
/// # Notes
/// - This function is designed to work in a multi-node MPI environment but can also run
///   on a single node for testing purposes.
/// - The root process handles task distribution and result aggregation.
#[cfg(feature = "mpi")]
pub fn run(
    tasks: Vec<Task>,
    ranked_list1_path: &str,
    ranked_list2_path: &str,
    population_size: usize,
    num_threads: usize,
) -> Vec<OptimizationResult> {
    // initialize MPI
    let universe = mpi::initialize().expect("Failed to initialize MPI");
    let world = universe.world();

    let rank = world.rank();
    eprint!("MPI Rank: {}\n", rank);

    let n_processes = world.size();
    eprint!("MPI Size: {}\n", n_processes);

    // Step 1: Prepare filenames and population size
    let file1 = if rank == 0 { ranked_list1_path.to_string() } else { String::new() };
    let file2 = if rank == 0 { ranked_list2_path.to_string() } else { String::new() };
    let mut population_size = if rank == 0 { population_size } else { 0 };

    // Step 2: Broadcast filename sizes
    let mut size1 = if rank == 0 { file1.len() } else { 0 };
    let mut size2 = if rank == 0 { file2.len() } else { 0 };

    world.process_at_rank(0).broadcast_into(&mut size1);
    world.process_at_rank(0).broadcast_into(&mut size2);

    // Step 3: Allocate buffers for filenames
    let mut buffer1 = vec![0; size1];
    let mut buffer2 = vec![0; size2];

    if rank == 0 {
        buffer1.copy_from_slice(file1.as_bytes());
        buffer2.copy_from_slice(file2.as_bytes());
    }

    // Step 4: Broadcast filenames and population size
    world.process_at_rank(0).broadcast_into(&mut buffer1);
    world.process_at_rank(0).broadcast_into(&mut buffer2);
    world.process_at_rank(0).broadcast_into(&mut population_size);

    // Step 5: Convert buffers back to strings
    let file1 = String::from_utf8(buffer1).expect("Failed to convert buffer1 to string");
    let file2 = String::from_utf8(buffer2).expect("Failed to convert buffer2 to string");

    eprintln!(
        "Rank {}: Received files: {}, {} with population size {}",
        rank, file1, file2, population_size
    );

    // Step 6: Read ranked feature lists locally
    let ranked_feature_list1 = read_ranked_feature_list_from_csv(&file1);
    let ranked_feature_list2 = read_ranked_feature_list_from_csv(&file2);

    // Step 7: Distribute tasks
    let task_chunks: Vec<Vec<Task>> = if rank == 0 {
        tasks
            .chunks((tasks.len() + n_processes as usize - 1) / n_processes as usize)
            .map(|chunk| chunk.to_vec())
            .collect()
    } else {
        Vec::new()
    };

    // Step 8: Send task chunks to worker processes
    let chunk = if rank == 0 {
        for (i, task_chunk) in task_chunks.iter().enumerate().skip(1) {
            let serialized_chunk = bincode::serialize(task_chunk).expect("Serialization failed");
            world.process_at_rank(i as i32).send(&serialized_chunk);
        }
        task_chunks.get(0).cloned().unwrap_or_default()
    } else {
        let (received_data, _status) = world.process_at_rank(0).receive_vec::<u8>();
        bincode::deserialize(&received_data).expect("Deserialization failed")
    };

    // print the rank and number of items in the chunk
    eprintln!("Rank {}: Chunk size: {}", rank, chunk.len());

    // Step 9: Run tasks in parallel using `run_single_node`
    let results = run_single_node(
        chunk,
        ranked_feature_list1.clone(),
        ranked_feature_list2.clone(),
        population_size,
        num_threads,
    );

    // Step 10: Gather results at root
    if rank == 0 {
        let mut all_results = results;
        for i in 1..n_processes {
            let (received_data, _status) = world.process_at_rank(i).receive_vec::<u8>();
            let worker_results: Vec<OptimizationResult> =
                bincode::deserialize(&received_data).expect("Deserialization failed");
            all_results.extend(worker_results);
        }
        all_results
    } else {
        let serialized_results = bincode::serialize(&results).expect("Serialization failed");
        world.process_at_rank(0).send(&serialized_results);
        Vec::new()
    }
}

#[cfg(not(feature = "mpi"))]
pub fn run(
    _tasks: Vec<Task>,
    _ranked_list1_path: &str,
    _ranked_list2_path: &str,
    _population_size: usize,
    _num_threads: usize,
) -> Vec<OptimizationResult> {
    panic!("Multi-node mode requires MPI. Rebuild with the 'mpi' feature enabled.");
}