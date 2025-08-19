use std::path::PathBuf;

use clap::Parser;

#[cfg(feature = "mpi")]
use dual_threshold_optimization::run::run_multi_node;

use dual_threshold_optimization::dto::compute_population_size;
use dual_threshold_optimization::read::{
    read_feature_list_from_file, read_ranked_feature_list_from_csv,
};
use dual_threshold_optimization::run::run_single_node;
use dual_threshold_optimization::run::Task;
use dual_threshold_optimization::stat_operations::empirical_pvalue;

/// Dual Threshold Optimization CLI
///
/// This program allows users to optimize thresholds for ranked feature lists,
/// optionally using permutations and parallelization.
#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Path to the first ranked feature list (CSV format).Rank 1 is the most important feature.
    ///
    /// This should have two columns: "feature" and "rank". There should be
    /// **NO HEADER**.  
    ///
    /// Rank is expected to be an integer. It is recommended that
    /// ties are handed with the `min` or `max` method.
    #[arg(short = '1', long, value_name = "FILE")]
    ranked_list1: PathBuf,

    /// Path to the second ranked feature list (CSV format). Rank 1 is the most important feature.
    ///
    /// This should have two columns: "feature" and "rank". There should be
    /// **NO HEADER**.  
    ///
    /// Rank is expected to be an integer. It is recommended that
    /// ties are handed with the `min` or `max` method.
    #[arg(short = '2', long, value_name = "FILE")]
    ranked_list2: PathBuf,

    /// Path to the background feature list (one feature per line, optional)
    #[arg(short, long, value_name = "FILE")]
    background: Option<PathBuf>,

    /// Number of permutations to perform
    #[arg(short, long, value_parser, default_value_t = 1000)]
    permutations: usize,

    /// Number of threads to use per task. For single-node parallelization, this will
    /// be the number of threads available on the machine.
    ///
    /// For multi-node parallelization, this will be the number of threads per task.
    ///
    /// Example: If you submit via Slurm with `-ntasks 4 --cpus-per-task 10`,
    /// then this value should be set to 10. This configuration will run 40
    /// permutations in parallel.
    #[arg(short, long, value_parser, default_value_t = 1)]
    threads: usize,

    /// Enable multi-node mode using MPI. This requires that the program has been built
    /// with the `mpi` feature enabled.
    #[arg(short = 'm', long, action = clap::ArgAction::SetTrue)]
    multi_node: bool,
}

fn main() {
    let cli = Cli::parse();

    let mut threads = cli.threads;

    if threads == 0 {
        eprintln!("Warning: Number of threads cannot be 0. Setting threads to 1.");
        threads = 1;
    }

    eprintln!("Ranked list 1: {}", cli.ranked_list1.display());
    eprintln!("Ranked list 2: {}", cli.ranked_list2.display());
    eprintln!("Permutations: {}", cli.permutations);
    eprintln!("Threads: {}", threads);
    eprintln!(
        "Multi-node mode: {}",
        if cli.multi_node {
            "enabled"
        } else {
            "disabled"
        }
    );

    let ranked_feature_list1 = read_ranked_feature_list_from_csv(
        cli.ranked_list1
            .to_str()
            .expect("Invalid file path for input1"),
    );
    let ranked_feature_list2 = read_ranked_feature_list_from_csv(
        cli.ranked_list2
            .to_str()
            .expect("Invalid file path for input2"),
    );

    eprintln!(
        "The product of the lengths of the threshold lists \
        (this describes the asymptotic runtime of a single job): {}",
        ranked_feature_list1.thresholds().len() * ranked_feature_list2.thresholds().len()
    );

    let background = cli.background.map(|background_file| {
        read_feature_list_from_file(
            background_file
                .to_str()
                .expect("Invalid file path for background"),
        )
    });

    let population_size = compute_population_size(
        &ranked_feature_list1,
        &ranked_feature_list2,
        background.as_ref(),
    );

    let mut tasks: Vec<Task> = vec![Task {
        id: 0,
        permute: false,
    }];
    tasks.extend((1..=cli.permutations).map(|id| Task { id, permute: true }));

    // If the mpi feature is used to build, and --multi-node is set, then run in
    // the multi-node mode. Otherwise, this will run in single-node mode.
    if cli.multi_node {
        #[cfg(feature = "mpi")]
        {
            let results = run_multi_node(
                tasks,
                cli.ranked_list1.to_str().unwrap(),
                cli.ranked_list2.to_str().unwrap(),
                population_size,
                threads,
            );

            // if results is not an empty vector, then get the final results and print
            // them out
            if !results.is_empty() {
                let final_result = empirical_pvalue(results);
                println!("{}", serde_json::to_string_pretty(&final_result).unwrap());
            }
            return;
        }

        #[cfg(not(feature = "mpi"))]
        {
            panic!("Multi-node mode requires MPI. Rebuild with the 'mpi' feature enabled.");
        }
    }

    let results = run_single_node(
        tasks,
        ranked_feature_list1,
        ranked_feature_list2,
        population_size,
        threads,
    );

    let final_result = empirical_pvalue(results);
    println!("{}", serde_json::to_string_pretty(&final_result).unwrap());
}
