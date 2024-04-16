use peroxide::fuga::GaussLegendre;
use peroxide::numerical::integral::integrate;
use rand::seq::SliceRandom;
use rand_distr::{Distribution, Normal};
use rayon::prelude::*;

/// Struct representing an agent in the simulation.
///
/// # Fields
/// * `age` - The current age of the agent in years.
/// * `aging_parameters` - A vector of floating-point numbers representing parameters for the aging function
/// * `learning_parameters` - A vector of floating-point numbers representing parameters for learning benefit function.
/// * `growth_parameters` - A vector of floating-point numbers that represent growth benefit function.
#[derive(Clone)]
pub struct Agent {
    pub age: f64,
    pub aging_parameters: Vec<f64>,
    pub learning_parameters: Vec<f64>,
    pub growth_parameters: Vec<f64>,
}

/// Initialize a population of agents with random parameters based on normal distributions.
///
/// # Arguments
/// * `initial_population_size` - The number of agents to create in the population.
/// * `aging_parameters` - A vector of floating-point numbers representing parameters for the aging function.
/// * `learning_parameters` - A vector of floating-point numbers representing parameters for the learning benefit function.
/// * `growth_parameters` - A vector of floating-point numbers that represent growth benefit function.
/// * `initial_age_distribution` - A two-element array representing the mean and standard deviation of the initial age distribution.
/// * `initial_b_distribution` - A two-element array representing the mean and standard deviation of the initial b distribution.
/// * `initial_lmax_distribution` - A two-element array representing the mean and standard deviation of the initial lmax distribution.
///
/// # Returns
/// Returns a vector of agents with random parameters based on the specified distributions.
pub fn initialize_population<'a, 'b>(
    initial_population_size: usize,
    aging_parameters: &[f64],
    learning_parameters: &[f64],
    growth_parameters: &[f64],
    initial_age_distribution: [f64; 2],
    initial_b_distribution: [f64; 2],
    initial_lmax_distribution: [f64; 2],
) -> Vec<Agent> {
    let mut population = Vec::with_capacity(initial_population_size);

    let age_dist = Normal::new(initial_age_distribution[0], initial_age_distribution[1]).unwrap();
    let b_dist = Normal::new(initial_b_distribution[0], initial_b_distribution[1]).unwrap();
    let lmax_dist =
        Normal::new(initial_lmax_distribution[0], initial_lmax_distribution[1]).unwrap();

    for _ in 0..initial_population_size {
        let age: f64 = age_dist.sample(&mut rand::thread_rng()).max(0.0).round();
        let b = b_dist.sample(&mut rand::thread_rng()).max(0.0);
        let lmax = lmax_dist.sample(&mut rand::thread_rng()).max(0.0);

        let mut agent_aging_parameters = aging_parameters.to_owned();
        agent_aging_parameters[1] = b;

        let mut agent_learning_parameters = learning_parameters.to_owned();
        agent_learning_parameters[0] = lmax;

        let agent_growth_parameters = growth_parameters.to_owned();

        let agent = Agent {
            age,
            aging_parameters: agent_aging_parameters,
            learning_parameters: agent_learning_parameters,
            growth_parameters: agent_growth_parameters,
        };

        population.push(agent);
    }
    population
}

/// Calculate the probability of death for an agent over a given time step using an aging model.
///
/// This function integrates an aging model over a specified time interval to estimate the probability of death.
/// It utilizes a custom closure that combines aging, learning, and growth parameters of the agent.
///
/// # Arguments
/// * `agent` - A reference to an Agent struct containing the agent's parameters and current age.
/// * `time_step` - The time interval over which to calculate the probability of death.
/// * `aging_intermediate_closure` - A closure that computes aging effects using three different parameter sets.
///
/// # Returns
/// Returns the estimated probability of death for the agent over the specified time step.
pub fn get_proba_of_death_agent(
    agent: &Agent,
    time_step: f64,
    aging_intermediate_closure: &dyn Fn(f64, &[f64], &[f64], &[f64]) -> f64,
) -> f64 {
    integrate(
        |x: f64| -> f64 {
            aging_intermediate_closure(
                x,
                &agent.aging_parameters,
                &agent.learning_parameters,
                &agent.growth_parameters,
            )
        },
        (agent.age, agent.age + time_step),
        GaussLegendre(5),
    )
}

/// Determine whether an agent dies during a given time step based on mortality probabilities.
///
/// # Arguments
/// * `agent` - A reference to the `Agent` struct, representing the agent whose death is being evaluated.
/// * `time_step` - The time interval over which to calculate the probability of death.
/// * `aging_intermediate_closure` - A closure that computes the aging effects using the agent's parameters.
///
/// # Returns
/// Returns `true` if the agent is dead by the end of the time step, otherwise `false`.
pub fn get_death_agent(
    agent: &Agent,
    time_step: f64,
    aging_intermediate_closure: &dyn Fn(f64, &[f64], &[f64], &[f64]) -> f64,
) -> bool {
    let proba_of_death = get_proba_of_death_agent(agent, time_step, aging_intermediate_closure);
    rand::random::<f64>() < proba_of_death
}

/// Determine which agents in a population die during a given time step and removes them from the population.
///
/// # Arguments
/// * `population` - A mutable reference to a vector of `Agent` structs representing the population of agents.
/// * `time_step` - The time interval over which to calculate the probability of death.
/// * `aging_intermediate_closure` - A closure that computes the aging effects using the agent's parameters.
pub fn get_death_population<F: Fn(f64, &[f64], &[f64], &[f64]) -> f64 + Send + Sync>(
    population: &mut Vec<Agent>,
    time_step: f64,
    aging_intermediate_closure: &F,
) {
    let death_test_parallel = population
        .par_iter()
        .map(|agent| get_death_agent(agent, time_step, aging_intermediate_closure))
        .collect::<Vec<_>>();
    let mut dead_agent_indexes: Vec<usize> = death_test_parallel
        .iter()
        .enumerate()
        .filter(|&(_, &value)| value)
        .map(|(index, _)| index)
        .collect();

    dead_agent_indexes.sort();
    dead_agent_indexes.reverse();

    for index in dead_agent_indexes.iter() {
        population.swap_remove(*index);
    }
}

/// Increment the age of all agents in the population by a specified time step.
///
/// # Arguments
/// * `population` - A mutable reference to a vector of `Agent` structs representing the population of agents.
/// * `time_step` - The time step by which to increment the age of all agents.
pub fn increment_age_population(population: &mut Vec<Agent>, time_step: f64) {
    for agent in population.iter_mut() {
        agent.age += time_step;
    }
}

/// Potentially mutate a parameter by adding a random value sampled from a normal distribution.
///
/// # Arguments
/// * `param` - A mutable reference to the parameter to mutate.
/// * `mutation_rate` - The probability of mutation for the parameter.
/// * `mutation_strength` - The standard deviation of the normal distribution used to sample the mutation value.
pub fn mutate_parameter(param: &mut f64, mutation_rate: f64, mutation_strength: f64) {
    if rand::random::<f64>() < mutation_rate {
        let mutation_dist = Normal::new(*param, mutation_strength).unwrap();
        *param = mutation_dist.sample(&mut rand::thread_rng()).max(0.0);
    }
}

/// Clone an agent and potentially mutate its parameters.
///
/// # Arguments
/// * `agent` - A reference to the `Agent` struct to clone.
/// * `mutable_b` - A boolean indicating whether the `b` parameter is mutable.
/// * `mutable_lmax` - A boolean indicating whether the `lmax` parameter is mutable.
/// * `b_mutation_rate` - The probability of mutation for the `b` parameter.
/// * `lmax_mutation_rate` - The probability of mutation for the `lmax` parameter.
/// * `b_mutation_strength` - The standard deviation of the normal distribution used to sample the mutation value for `b`.
/// * `lmax_mutation_strength` - The standard deviation of the normal distribution used to sample the mutation value for `lmax`.
///
/// # Returns
/// Returns a new `Agent` struct that is a clone of the input agent with potentially mutated parameters.
pub fn clone_agent(
    agent: &Agent,
    mutable_b: bool,
    mutable_lmax: bool,
    b_mutation_rate: f64,
    lmax_mutation_rate: f64,
    b_mutation_strength: f64,
    lmax_mutation_strength: f64,
) -> Agent {
    let mut b = agent.aging_parameters[1];
    if mutable_b {
        mutate_parameter(&mut b, b_mutation_rate, b_mutation_strength);
    }
    let mut lmax = agent.learning_parameters[0];
    if mutable_lmax {
        mutate_parameter(&mut lmax, lmax_mutation_rate, lmax_mutation_strength);
    }

    let mut agent_aging_parameters = agent.aging_parameters.to_owned();
    agent_aging_parameters[1] = b;

    let mut agent_learning_parameters = agent.learning_parameters.to_owned();
    agent_learning_parameters[0] = lmax;

    let agent_growth_parameters = agent.growth_parameters.to_owned();

    Agent {
        age: 0.0,
        aging_parameters: agent_aging_parameters,
        learning_parameters: agent_learning_parameters,
        growth_parameters: agent_growth_parameters,
    }
}

/// Generate new agents through asexual reproduction to maintain a constant population size.
///
/// # Arguments
/// * `population` - A mutable reference to a vector of `Agent` structs representing the population of agents.
/// * `population_cap` - The maximum population size.
/// * `mutable_b` - A boolean indicating whether the `b` parameter is mutable.
/// * `mutable_lmax` - A boolean indicating whether the `lmax` parameter is mutable.
/// * `b_mutation_rate` - The probability of mutation for the `b` parameter.
/// * `lmax_mutation_rate` - The probability of mutation for the `lmax` parameter.
/// * `b_mutation_strength` - The standard deviation of the normal distribution used to sample the mutation value for `b`.
/// * `lmax_mutation_strength` - The standard deviation of the normal distribution used to sample the mutation value for `lmax`.
pub fn get_reproduction_population(
    population: &mut Vec<Agent>,
    population_cap: usize,
    mutable_b: bool,
    mutable_lmax: bool,
    b_mutation_rate: f64,
    lmax_mutation_rate: f64,
    b_mutation_strength: f64,
    lmax_mutation_strength: f64,
) {
    let missing_population = population_cap as usize - population.len();
    let par_iter = (0..missing_population).into_par_iter();
    let new_babies: Vec<Agent> = par_iter
        .map(|_| {
            let parent = population.choose(&mut rand::thread_rng()).unwrap();
            clone_agent(
                parent,
                mutable_b,
                mutable_lmax,
                b_mutation_rate,
                lmax_mutation_rate,
                b_mutation_strength,
                lmax_mutation_strength,
            )
        })
        .collect();
    population.extend(new_babies);
}

/// Calculate the mean and variance of the b parameter in the population.
///
/// # Arguments
/// * `population` - A reference to a vector of `Agent` structs representing the population of agents.
///
/// # Returns
/// Returns a tuple containing the mean and variance of the b parameter in the population.
pub fn get_population_b_stats(population: &Vec<Agent>) -> (f64, f64) {
    let b_values = population
        .iter()
        .map(|agent| agent.aging_parameters[1])
        .collect::<Vec<_>>();

    let b_mean = b_values.par_iter().sum::<f64>() / b_values.len() as f64;

    let b_variance = b_values
        .par_iter()
        .map(|b| (b - b_mean).powi(2))
        .sum::<f64>()
        / b_values.len() as f64;

    (b_mean, b_variance)
}

/// Calculate the mean and variance of the lmax parameter in the population.
///
/// # Arguments
/// * `population` - A reference to a vector of `Agent` structs representing the population of agents.
///
/// # Returns
/// Returns a tuple containing the mean and variance of the lmax parameter in the population.
pub fn get_population_lmax_stats(population: &Vec<Agent>) -> (f64, f64) {
    let lmax_values = population
        .iter()
        .map(|agent| agent.learning_parameters[0])
        .collect::<Vec<_>>();

    let lmax_mean = lmax_values.par_iter().sum::<f64>() / lmax_values.len() as f64;

    let lmax_variance = lmax_values
        .par_iter()
        .map(|lmax| (lmax - lmax_mean).powi(2))
        .sum::<f64>()
        / lmax_values.len() as f64;

    (lmax_mean, lmax_variance)
}
