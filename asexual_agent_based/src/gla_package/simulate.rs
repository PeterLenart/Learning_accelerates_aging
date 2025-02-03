use std::fs::File;
use csv::Writer;
use indicatif::{ProgressBar, ProgressStyle};

use crate::gla_package::agent_based::{
    get_death_population, get_population_b_stats, get_population_lmax_stats, get_reproduction_population,
    increment_age_population, initialize_population,
};

/// Struct to store the results of a simulation
/// 
/// # Fields
/// * `mean_b` - The mean value of the b parameter in the population
/// * `mean_lmax` - The mean value of the lmax parameter in the population
/// * `time` - The time current time of the simulation
/// * `replicate_id` - The ID of the replicate
#[derive(serde::Serialize)]
struct SimulationResult {
    mean_b: f64,
    mean_lmax: f64,
    time: f64,
    replicate_id: i32,
}

/// Run a simulation of the agent-based model
pub fn run_simulation(
    output_writer: &mut Writer<File>,
    population_cap: usize,
    simulation_time:usize,
    replicate_id: i32,
    aging_parameters: &[f64],
    learning_parameters: &[f64],
    growth_parameters: &[f64],
    initial_age_distribution: [f64; 2],
    initial_b_distribution: [f64; 2],
    initial_lmax_distribution: [f64; 2],
    time_step: f64,
    mutable_b: bool,
    mutable_lmax: bool,
    b_mutation_rate: f64,
    lmax_mutation_rate: f64,
    b_mutation_strength: f64,
    lmax_mutation_strength: f64,
    aging_intermediate_closure: impl Fn(f64, &[f64], &[f64], &[f64]) -> f64 + Send + Sync,
) {

    /// Initialize the population with the provided parameters
    // let mut wtr = Writer::from_path("foo.csv").unwrap();
    let mut population = initialize_population(
        population_cap,
        aging_parameters,
        learning_parameters,
        growth_parameters,
        initial_age_distribution,
        initial_b_distribution,
        initial_lmax_distribution,);

    /// Create a progress bar
    let bar = ProgressBar::new(simulation_time as u64);
    bar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:50.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-"),
    );

    /// Run the simulation
    for i in 0..simulation_time {
        /// Remove agents that die during the current time step
        get_death_population(&mut population, time_step, &aging_intermediate_closure);
        /// Create new agents through cloning until the population reaches the cap 
        get_reproduction_population(
            &mut population,
            population_cap,
            mutable_b,
            mutable_lmax,
            b_mutation_rate,
            lmax_mutation_rate,
            b_mutation_strength,
            lmax_mutation_strength,
        );
        /// Increment the age of all agents in the population
        increment_age_population(&mut population, time_step);
        /// Get the statistics for the current population
        let b_stats = get_population_b_stats(&population);
        let lmax_stats = get_population_lmax_stats(&population);

        /// Write the results to the output file
        let res = SimulationResult {
            mean_b: b_stats.0,
            mean_lmax: lmax_stats.0,
            time: (i as f64) * time_step,
            replicate_id: replicate_id,
        };
        let _ = output_writer.serialize(res);
        bar.inc(1);
    }
    bar.finish();
}
