mod gla_package;
use csv::Writer;
use crate::gla_package::{gla::{
    aging_gompertz_makeham, fertility_brass_polynomial, find_maximum_fertility, gla_model,
    growth_function, learning_function, constant_fertility, mortality_improvement_function, toy_model, _aging_gompertz,step_fertility
}, simulate::run_simulation};

// Define a trait for the aging closure to handle different function signatures
trait AgingClosure: Send + Sync {
    fn call(&self, x: f64, aging_parameters: &[f64], secondary_parameters: &[f64], growth_parameters: &[f64]) -> f64;
}

// Implement for the GLA model closure
struct GlaClosure;
impl AgingClosure for GlaClosure {
    fn call(&self, x: f64, aging_parameters: &[f64], learning_parameters: &[f64], growth_parameters: &[f64]) -> f64 {
        gla_model(
            x,
            aging_gompertz_makeham as fn(f64, &[f64]) -> f64,
            learning_function,
            growth_function,
            aging_parameters,
            learning_parameters,
            growth_parameters,
            1e-5, // minimum_mortality
        )
    }
}

// Implement for the toy model closure
struct ToyClosure;
impl AgingClosure for ToyClosure {
    fn call(&self, x: f64, aging_parameters: &[f64], improvement_parameters: &[f64], _growth_parameters: &[f64]) -> f64 {
        toy_model(
            x,
            _aging_gompertz as fn(f64, &[f64]) -> f64,
            mortality_improvement_function,
            aging_parameters,
            improvement_parameters,
            1e-5, // minimum_mortality
        )
    }
}

// Define a trait for the fertility closure
trait FertilityClosure {
    fn call(&self, x: f64) -> f64;
}

struct NormalizedFertilityClosure {
    fertility_function: fn(f64, &[f64]) -> f64,
    fertility_parameters: Vec<f64>,
    maximum_fertility: f64,
}

impl FertilityClosure for NormalizedFertilityClosure {
    fn call(&self, x: f64) -> f64 {
        ((self.fertility_function)(x, &self.fertility_parameters) / self.maximum_fertility).min(1.0)
    }
}

// Constant fertility closure
struct ConstantFertilityClosure;

impl FertilityClosure for ConstantFertilityClosure {
    fn call(&self, _x: f64) -> f64 {
        1.0
    }
}

fn main() {
    // Define overall parameters
    let time_step = 1.0;
    let initial_female_proportion = 0.5;
    let minimum_mortality = 1e-5;
    let use_toy_model = false;
    let fertility_model = "brass_polynomial";

    // Define initial distributions
    let initial_age_distribution = [20.0, 10.0];
    let mut initial_b_distribution = [0.14, 0.005];
    let mut initial_lmax_distribution = [0.01606792505529796, 0.0];
    let mut initial_gmax_distribution = [0.05168141300917714, 0.0]; 

    // // Define initial distributions for the toy model
    // let initial_age_distribution = [20.0, 10.0];
    // let mut initial_b_distribution = [0.1, 0.0];
    // let mut initial_lmax_distribution = [5.0, 0.0];
    // let initial_gmax_distribution = [1.0, 0.0];

    // Define simulation parameters
    let population_cap = 10000;
    let simulation_time : usize = 600000;
    let replicate_number = 5;
    let assortative_mating = false;
    let remove_non_reproducing = true;
    let tradeoff = false;
    let start_b = initial_b_distribution[0];

    // Define mutation parameters
    let mutable_b = true;
    let mutable_lmax = false;
    let mutable_gmax = false;
    let b_mutation_rate: f64 = 0.02;
    let lmax_mutation_rate: f64 = 0.02;
    let gmax_mutation_rate: f64 = 0.02;
    let b_mutation_strength = 0.012;
    let lmax_mutation_strength = 0.012;
    let gmax_mutation_strength = 0.012;

    // Define parameters based on model type
    let (aging_parameters, learning_parameters, growth_parameters, aging_closure_box): 
    ([f64; 3], [f64; 3], [f64; 2], Box<dyn AgingClosure>) = 
    if !use_toy_model {
        println!("Running GLA model!");

        // Define GLA parameters
        let aging_params = [0.00275961297460256, 0.04326224872667336, 0.025201676835511704];
        let learning_params = [0.01606792505529796, 39.006865144958745, 0.11060749334680318];
        let growth_params: [f64; 2] = [0.05168141300917714, 0.08765165352033985];

        (aging_params, learning_params, growth_params, Box::new(GlaClosure))
    } else {
        println!("Running Toy model!");
        // Define toy model parameters
        let aging_params = [0.01, 0.1, 0.0]; // Pad to match [f64; 3]
        let mut learning_params = [5.0, f64::INFINITY, 0.2];
        let growth_params = [1.0, 1.0];

        (aging_params, learning_params, growth_params, Box::new(ToyClosure))
    };

    // Fertility model selection
    let (normalized_male_fertility_closure_box, normalized_female_fertility_closure_box, male_menopause, female_menopause) = 
        if fertility_model == "brass_polynomial" {
            println!("Using brass polynomial fertility model!");
            let female_fertility_parameters = [2.445e-5, 14.8, 32.836];
            // let male_fertility_parameters = [0.00000978, 14.8, 47.836];
            let male_fertility_parameters = female_fertility_parameters;
    
            let female_fertility_function = fertility_brass_polynomial;
            let male_fertility_function = fertility_brass_polynomial;
    
            let female_menopause = female_fertility_parameters[1] + female_fertility_parameters[2];
            let male_menopause = male_fertility_parameters[1] + male_fertility_parameters[2];
    
            let female_maximum_fertility = find_maximum_fertility(
                &female_fertility_function,
                &female_fertility_parameters,
                20.0,
            );
            let male_maximum_fertility =
                find_maximum_fertility(&male_fertility_function, &male_fertility_parameters, 20.0);

            let normalized_male_fertility_closure_box: Box<dyn FertilityClosure> = Box::new(NormalizedFertilityClosure {
                fertility_function: male_fertility_function,
                fertility_parameters: male_fertility_parameters.to_vec(),
                maximum_fertility: male_maximum_fertility,
            });

            let normalized_female_fertility_closure_box: Box<dyn FertilityClosure> = Box::new(NormalizedFertilityClosure {
                fertility_function: female_fertility_function,
                fertility_parameters: female_fertility_parameters.to_vec(),
                maximum_fertility: female_maximum_fertility,
            });

            (normalized_male_fertility_closure_box, normalized_female_fertility_closure_box, male_menopause, female_menopause)
        } else if fertility_model == "constant" {
            println!("Using constant fertility model!");
            let female_fertility_parameters = [1.0];
            let male_fertility_parameters = [1.0];

            let female_menopause : f64 = 0f64;
            let male_menopause : f64 = 0f64;

            let normalized_male_fertility_closure_box: Box<dyn FertilityClosure> = Box::new(ConstantFertilityClosure);
            let normalized_female_fertility_closure_box: Box<dyn FertilityClosure> = Box::new(ConstantFertilityClosure);

            (normalized_male_fertility_closure_box, normalized_female_fertility_closure_box, male_menopause, female_menopause)
        } else if fertility_model == "step" {
            let female_fertility_parameters = [1.0, 25.0, 100.0];
            let male_fertility_parameters = [1.0, 25.0, 100.0];

            let female_fertility_function = step_fertility;
            let male_fertility_function = step_fertility;

            let female_menopause = female_fertility_parameters[2];
            let male_menopause = male_fertility_parameters[2];

            let female_maximum_fertility = female_fertility_parameters[0];
            let male_maximum_fertility = male_fertility_parameters[0];

            let normalized_male_fertility_closure_box: Box<dyn FertilityClosure> = Box::new(NormalizedFertilityClosure {
                fertility_function: male_fertility_function,
                fertility_parameters: male_fertility_parameters.to_vec(),
                maximum_fertility: male_maximum_fertility,
            });

            let normalized_female_fertility_closure_box: Box<dyn FertilityClosure> = Box::new(NormalizedFertilityClosure {
                fertility_function: female_fertility_function,
                fertility_parameters: female_fertility_parameters.to_vec(),
                maximum_fertility: female_maximum_fertility,
            });

            (normalized_male_fertility_closure_box, normalized_female_fertility_closure_box, male_menopause, female_menopause)
        } else {
            panic!("Unknown fertility model!");
        };

    // Helper function to mimic numpy's linspace
    fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
        let step = (end - start) / (n - 1) as f64;
        (0..n).map(|i| start + step * i as f64).collect()
    }

    // let start_b_grid = linspace(1e-5, 0.1, 11);
    // // let start_improvement_grid = linspace(0.0, 60.0, 11);
    // // let start_improvement_grid = Vec::from([0.0]);
    // let improvement_strength_grid = linspace(0.0, 0.5, 11);
    // // let improvement_strength_grid = Vec::from([0.0]);

    // for &start_b in start_b_grid.iter(){
    //     // for &start_improvement in start_improvement_grid.iter(){
    //     for &improvement_strength in improvement_strength_grid.iter(){
    //         println!("######################################");
    //         // println!("###### Simulation with b = {} and improvement start = {} ######", start_b, start_improvement);
    //         println!("###### Simulation with b = {} and improvement strength = {} ######", start_b, improvement_strength);
    //         println!("######################################");
    //         initial_b_distribution = [start_b, 0.0];
    //         // initial_lmax_distribution = [start_improvement, 0.0];
    //         learning_parameters[2] = improvement_strength;

    //         // let base_name_part = format!("early_slope_toy_model_step_fertility_varying_start_and_b_{}_{}", initial_lmax_distribution[0], initial_b_distribution[0]);
    //         // let output_file_name = format!("./simulation_results/toy_model_simulations/varying_improvement_strength_late_40/{}_{}_{}_{}_{}_{}.csv", base_name_part, mating_name_part, learning_name_part, removal_name_part, tradeoff_name_part,learning_parameters[2]);

    //         let base_name_part = format!("toy_model_start_{}_b_{}", initial_lmax_distribution[0], initial_b_distribution[0]);
    //         let output_file_name = format!("./simulation_results/toy_model_simulations/varying_improvement_strength_early_5/{}_{}.csv", base_name_part, learning_parameters[2]);
    //         let mut wtr = Writer::from_path(output_file_name).unwrap();

    //         for i in 0..replicate_number{
    //             println!("Replicate : {}/{}", i+1, replicate_number);
    //             run_simulation(&mut wtr, population_cap, simulation_time, i, assortative_mating, &aging_parameters, &learning_parameters, &growth_parameters, initial_age_distribution, initial_b_distribution, initial_lmax_distribution, initial_gmax_distribution, initial_female_proportion, time_step, mutable_b, mutable_lmax, mutable_gmax, b_mutation_rate, lmax_mutation_rate, gmax_mutation_rate, b_mutation_strength, lmax_mutation_strength, gmax_mutation_strength, aging_intermediate_closure, &normalized_male_fertility_closure, &normalized_female_fertility_closure, tradeoff, start_b, remove_non_reproducing, male_menopause, female_menopause)
    //             }
    //         }
    //     }

    let start_b_grid = Vec::from([0.14]);
    let start_gmax_grid = linspace(0.055, 0.175, 5);

    let output_folder = "./simulation_results/gla_simulations/gmax_simulations/gmax_plateau/";
    std::fs::create_dir_all(output_folder).unwrap();

    let aging_closure = |x: f64,
                        aging_parameters: &[f64],
                        learning_parameters: &[f64],
                        growth_parameters: &[f64]|
    -> f64 {
        aging_closure_box.call(x, aging_parameters, learning_parameters, growth_parameters)
    };

    // let normalized_male_fertility_closure = |x| normalized_male_fertility_closure.call(x);
    let normalized_male_fertility_closure = Box::new(|x| normalized_male_fertility_closure_box.call(x));
    // let normalized_female_fertility_closure = |x| normalized_female_fertility_closure.call(x);
    let normalized_female_fertility_closure = Box::new(|x| normalized_female_fertility_closure_box.call(x));




    for &start_b in start_b_grid.iter(){
        for &start_gmax in start_gmax_grid.iter(){
            println!("######################################");
            println!("###### Simulation with b = {} and gmax = {} ######", start_b, start_gmax);
            println!("######################################");
            initial_b_distribution = [start_b, 0.005];
            initial_gmax_distribution = [start_gmax, 0.0];

            let base_name_part = format!("plateau_varying_gmax_and_b_{}_{}.csv", initial_gmax_distribution[0], initial_b_distribution[0]);
            let output_file_name = format!("{}/{}", output_folder, base_name_part);

            let mut wtr = Writer::from_path(output_file_name).unwrap();

            for i in 0..replicate_number{
                println!("Replicate : {}/{}", i+1, replicate_number);
                run_simulation(&mut wtr, population_cap, simulation_time, i, assortative_mating, &aging_parameters, &learning_parameters, &growth_parameters, initial_age_distribution, initial_b_distribution, initial_lmax_distribution, initial_gmax_distribution, initial_female_proportion, time_step, mutable_b, mutable_lmax, mutable_gmax, b_mutation_rate, lmax_mutation_rate, gmax_mutation_rate, b_mutation_strength, lmax_mutation_strength, gmax_mutation_strength, aging_closure, &normalized_male_fertility_closure, &normalized_female_fertility_closure, tradeoff, start_b, remove_non_reproducing, male_menopause, female_menopause)
                }
            }
        }
}

