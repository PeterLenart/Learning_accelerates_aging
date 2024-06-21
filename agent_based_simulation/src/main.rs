mod gla_package;
use csv::Writer;
use crate::gla_package::{gla::{
    aging_gompertz_makeham, fertility_brass_polynomial, find_maximum_fertility, gla_model,
    growth_function, learning_function, constant_fertility, mortality_improvement_function, toy_model, _aging_gompertz,step_fertility
}, simulate::run_simulation};

fn main() {

    // Define overall parameters
    let time_step = 1.0;
    let initial_female_proportion = 0.5;
    let minimum_mortality = 1e-5;

    // Define GLA parameters
    let aging_parameters = [0.00275961297460256,0.04326224872667336,0.025201676835511704] ;
    let learning_parameters = [0.01606792505529796,39.006865144958745,0.11060749334680318];
    let growth_parameters: [f64; 2] = [0.05168141300917714,0.08765165352033985];

    // // Define toy model parameters
    // let aging_parameters = [0.01, 0.1];
    // // Hijacking the learning parameters as the improvement parameters for the toy model
    // // let learning_parameters = [0.0, f64::INFINITY, 0.5];
    // let learning_parameters = [35.0, f64::INFINITY, 0.2];
    // // The growth parameters are not used in the toy model
    // let growth_parameters = [1.0, 1.0];


    // Define fertility parameters
    // let female_fertility_parameters = [1.0];
    // let male_fertility_parameters = [1.0];

    // let female_fertility_function = constant_fertility;
    // let male_fertility_function = constant_fertility;

    // let female_menopause : f64 = 0f64;
    // let male_menopause : f64 = 0f64;

    // let female_maximum_fertility = 1.0;
    // let male_maximum_fertility = 1.0;

    // let female_fertility_parameters = [1.0, 15.0, 100.0];
    // let male_fertility_parameters = [1.0, 15.0, 100.0];

    // let female_fertility_function = step_fertility;
    // let male_fertility_function = step_fertility;

    // let female_menopause = female_fertility_parameters[2];
    // let male_menopause = male_fertility_parameters[2];

    // let female_maximum_fertility = female_fertility_parameters[0];
    // let male_maximum_fertility = male_fertility_parameters[0];

    let female_fertility_parameters = [2.445e-5, 14.8, 32.836];
    let male_fertility_parameters = [2.445e-5, 14.8, 32.836];
    let male_fertility_parameters = [0.00000978, 14.8, 47.836];

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

    // Define normalized fertility closures
    let normalized_male_fertility_closure = Box::new(|x: f64| -> f64 {
        (male_fertility_function(x, &male_fertility_parameters) / male_maximum_fertility).min(1.0)
    });

    let normalized_female_fertility_closure = Box::new(|x: f64| -> f64{
        (female_fertility_function(x, &female_fertility_parameters) / female_maximum_fertility)
            .min(1.0)
    });

    // Define aging intermediate closure
    let aging_intermediate_closure = |x: f64,
                                      aging_parameters: &[f64],
                                      learning_parameters: &[f64],
                                      growth_parameters: &[f64]|
     -> f64 {
        gla_model(
            x,
            aging_gompertz_makeham as fn(f64, &[f64]) -> f64,
            learning_function,
            growth_function,
            &aging_parameters,
            &learning_parameters,
            &growth_parameters,
            minimum_mortality,
        )
    };

    // // Dirty hack to use toy model, the growth parameters are not used
    // let aging_intermediate_closure = |x: f64,
    //                                   aging_parameters: &[f64],
    //                                   improvement_parameters: &[f64],
    //                                   growth_parameters: &[f64]|
    //  -> f64 {
    //     toy_model(
    //         x,
    //         _aging_gompertz as fn(f64, &[f64]) -> f64,
    //         mortality_improvement_function,
    //         &aging_parameters,
    //         &improvement_parameters,
    //         minimum_mortality,
    //     )
    // };
    
    // #######################################################
    // ############# Simulation with learning ################
    // #######################################################

    // Define initial distributions
    let initial_age_distribution = [20.0, 10.0];
    let initial_b_distribution = [0.14, 0.005];
    let mut initial_lmax_distribution = [0.15, 0.0];
    let mut initial_gmax_distribution = [0.05168141300917714, 0.0]; 
    // let mut initial_gmax_distribution = [0.055, 0.0]

    // // Define initial distributions for the toy model
    // let initial_age_distribution = [20.0, 10.0];
    // let mut initial_b_distribution = [0.1, 0.0];
    // let mut initial_lmax_distribution = [35.0, 0.0];
    // let initial_gmax_distribution = [1.0, 0.0];

    // Define simulation parameters
    let population_cap = 10000;
    let simulation_time : usize = 1000;
    let replicate_number = 100;
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

    // Boilerplate to name the output file based on the simulation parameters
    // let base_name_part = "plateau_toy_model_step_fertility";
    // let base_name_part = "toy_model_constant_fertility_0_10";
    // Name the output file based on the simulation parameters
    // let base_name_part = format!("early_slope_toy_model_step_fertility_varying_start_and_b_{}_{}", initial_lmax_distribution[0], initial_b_distribution[0]);

    // // let mut learning_name_part = "with_learning";
    let mut learning_name_part = "with_improvement";
    // let mut learning_name_part = "no_improvement";
    let mut mating_name_part = "random_mating";
    let mut removal_name_part = "non_reproducing_kept";
    let mut tradeoff_name_part = "no_tradeoff";

    if tradeoff{
        tradeoff_name_part = "tradeoff";
    }

    if assortative_mating{
        mating_name_part = "assortative_mating";
    }

    if remove_non_reproducing{
        removal_name_part = "non_reproducing_removed";
    }

    // // Run the simulation
    // println!("######################################");
    // println!("###### Simulation with learning ######");
    // println!("######################################");

    // let output_file_name = format!("/mnt/towbin.data/shared/plenart/agent_based_simulations/{}_{}_{}_{}_{}_{}.csv", base_name_part, mating_name_part, learning_name_part, removal_name_part, tradeoff_name_part,initial_lmax_distribution[0]);
    // let mut wtr = Writer::from_path(output_file_name).unwrap();

    // // println!("Saving to file: {}", output_file_name)
    // for i in 0..replicate_number{
    //     println!("Replicate : {}/{}", i+1, replicate_number);
    //     run_simulation(&mut wtr, population_cap, simulation_time, i, assortative_mating, &aging_parameters, &learning_parameters, &growth_parameters, initial_age_distribution, initial_b_distribution, initial_lmax_distribution, initial_gmax_distribution, initial_female_proportion, time_step, mutable_b, mutable_lmax, mutable_gmax, b_mutation_rate, lmax_mutation_rate, gmax_mutation_rate, b_mutation_strength, lmax_mutation_strength, gmax_mutation_strength, aging_intermediate_closure, &normalized_male_fertility_closure, &normalized_female_fertility_closure, tradeoff, start_b, remove_non_reproducing, male_menopause, female_menopause)
    // }

    // Helper function to mimic numpy's linspace
    fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
        let step = (end - start) / (n - 1) as f64;
        (0..n).map(|i| start + step * i as f64).collect()
    }


    // let start_b_grid = linspace(0.02, 0.1, 11);
    // // let start_improvement_grid = linspace(0.0, 60.0, 11);
    // let start_improvement_grid = Vec::from([0.0]);

    // for &start_b in start_b_grid.iter(){
    //     for &start_improvement in start_improvement_grid.iter(){
    //         println!("######################################");
    //         println!("###### Simulation with b = {} and improvement start = {} ######", start_b, start_improvement);
    //         println!("######################################");
    //         initial_b_distribution = [start_b, 0.0];
    //         initial_lmax_distribution = [start_improvement, 0.0];

    //         let base_name_part = format!("early_slope_toy_model_step_fertility_varying_start_and_b_{}_{}", initial_lmax_distribution[0], initial_b_distribution[0]);
    //         let output_file_name = format!("./simulation_results/toy_model_simulations/varying_start_no_improvement/{}_{}_{}_{}_{}_{}.csv", base_name_part, mating_name_part, learning_name_part, removal_name_part, tradeoff_name_part,learning_parameters[2]);

    //         let mut wtr = Writer::from_path(output_file_name).unwrap();

    //         for i in 0..replicate_number{
    //             println!("Replicate : {}/{}", i+1, replicate_number);
    //             run_simulation(&mut wtr, population_cap, simulation_time, i, assortative_mating, &aging_parameters, &learning_parameters, &growth_parameters, initial_age_distribution, initial_b_distribution, initial_lmax_distribution, initial_gmax_distribution, initial_female_proportion, time_step, mutable_b, mutable_lmax, mutable_gmax, b_mutation_rate, lmax_mutation_rate, gmax_mutation_rate, b_mutation_strength, lmax_mutation_strength, gmax_mutation_strength, aging_intermediate_closure, &normalized_male_fertility_closure, &normalized_female_fertility_closure, tradeoff, start_b, remove_non_reproducing, male_menopause, female_menopause)
    //             }
    //         }
        // }

    let start_b_grid = linspace(1e-5, 0.14, 11);
    // let start_improvement_grid = linspace(0.0, 60.0, 11);
    let start_lmax_grid = Vec::from([0.0]);

    for &start_b in start_b_grid.iter(){
        for &start_lmax in start_lmax_grid.iter(){
            println!("######################################");
            println!("###### Simulation with b = {} and lmax = {} ######", start_b, start_lmax);
            println!("######################################");
            initial_b_distribution = [start_b, 0.0];
            initial_lmax_distribution = [start_lmax, 0.0];

            let base_name_part = format!("early_slope_toy_model_step_fertility_varying_lmax_and_b_{}_{}", initial_lmax_distribution[0], initial_b_distribution[0]);
            let output_file_name = format!("./simulation_results/gla_simulations/lmax_simulations/varying_lmax_no_learning/{}_{}_{}_{}_{}.csv", base_name_part, mating_name_part, learning_name_part, removal_name_part, tradeoff_name_part);

            let mut wtr = Writer::from_path(output_file_name).unwrap();

            for i in 0..replicate_number{
                println!("Replicate : {}/{}", i+1, replicate_number);
                run_simulation(&mut wtr, population_cap, simulation_time, i, assortative_mating, &aging_parameters, &learning_parameters, &growth_parameters, initial_age_distribution, initial_b_distribution, initial_lmax_distribution, initial_gmax_distribution, initial_female_proportion, time_step, mutable_b, mutable_lmax, mutable_gmax, b_mutation_rate, lmax_mutation_rate, gmax_mutation_rate, b_mutation_strength, lmax_mutation_strength, gmax_mutation_strength, aging_intermediate_closure, &normalized_male_fertility_closure, &normalized_female_fertility_closure, tradeoff, start_b, remove_non_reproducing, male_menopause, female_menopause)
                }
            }
        }

    // println!("#########################################");
    // println!("###### Simulation without learning ######");
    // println!("#########################################");
    // initial_lmax_distribution = [0.0, 0.0];
    // // learning_name_part = "no_learning";
    // learning_name_part = "no_improvement";

    // let output_file_name = format!("./simulation_results/toy_model_simulations/{}_{}_{}_{}_{}.csv", base_name_part, mating_name_part, learning_name_part, removal_name_part, tradeoff_name_part);
    // let mut wtr = Writer::from_path(output_file_name).unwrap();

    // for i in 0..replicate_number{
    //     println!("Replicate : {}/{}", i+1, replicate_number);
    //     run_simulation(&mut wtr, population_cap, simulation_time, i, assortative_mating, &aging_parameters, &learning_parameters, &growth_parameters, initial_age_distribution, initial_b_distribution, initial_lmax_distribution, initial_gmax_distribution, initial_female_proportion, time_step, mutable_b, mutable_lmax, mutable_gmax, b_mutation_rate, lmax_mutation_rate, gmax_mutation_rate, b_mutation_strength, lmax_mutation_strength, gmax_mutation_strength, aging_intermediate_closure, &normalized_male_fertility_closure, &normalized_female_fertility_closure, tradeoff, start_b, remove_non_reproducing, male_menopause, female_menopause)
    // }

//     // ######################################################
//     // ############### Simulation with growth ###############
//     // ######################################################

//     // Define initial distributions
//     let initial_age_distribution = [20.0, 10.0];
//     let initial_b_distribution = [0.14, 0.005];
//     // let mut initial_lmax_distribution = [0.15, 0.0];
//     let initial_lmax_distribution = [0.01606792505529796, 0.0];
//     // let mut initial_gmax_distribution = [0.05168141300917714, 0.0]; 
//     let mut initial_gmax_distribution = [0.175, 0.0];

//     // Define simulation parameters
//     let population_cap = 10000;
//     let simulation_time : usize = 1000;
//     let replicate_number = 500;
//     let assortative_mating = false;
//     let remove_non_reproducing = true;
//     let tradeoff = false;
//     let start_b = initial_b_distribution[0];

//     // Define mutation parameters
//     let mutable_b = true;
//     let mutable_lmax = false;
//     let mutable_gmax = false;
//     let b_mutation_rate: f64 = 0.02;
//     let lmax_mutation_rate: f64 = 0.02;
//     let gmax_mutation_rate: f64 = 0.02;
//     let b_mutation_strength = 0.012;
//     let lmax_mutation_strength = 0.012;
//     let gmax_mutation_strength = 0.012;

//     // Boilerplate to name the output file based on the simulation parameters
//     // let base_name_part = "plateau_brass_polynomial_equal_both";
//     let base_name_part = "early_slope_brass_polynomial_equal_both";
//     let mut growth_name_part = "high_growth";
//     let mut mating_name_part = "random_mating";
//     let mut removal_name_part = "non_reproducing_kept";
//     let mut tradeoff_name_part = "no_tradeoff";

//     if tradeoff{
//         tradeoff_name_part = "tradeoff";
//     }

//     if assortative_mating{
//         mating_name_part = "assortative_mating";
//     }

//     if remove_non_reproducing{
//         removal_name_part = "non_reproducing_removed";
//     }

//     // Run the simulation
//     println!("######################################");
//     println!("###### Simulation with high Gmax ######");
//     println!("######################################");

//     let output_file_name = format!("./simulation_results/gmax_simulations/{}_{}_{}_{}_{}_{}.csv", base_name_part, mating_name_part, growth_name_part, removal_name_part, tradeoff_name_part,initial_gmax_distribution[0]);
//     let mut wtr = Writer::from_path(output_file_name).unwrap();

//     for i in 0..replicate_number{
//         println!("Replicate : {}/{}", i+1, replicate_number);
//         run_simulation(&mut wtr, population_cap, simulation_time, i, assortative_mating, &aging_parameters, &learning_parameters, &growth_parameters, initial_age_distribution, initial_b_distribution, initial_lmax_distribution, initial_gmax_distribution, initial_female_proportion, time_step, mutable_b, mutable_lmax, mutable_gmax, b_mutation_rate, lmax_mutation_rate, gmax_mutation_rate, b_mutation_strength, lmax_mutation_strength, gmax_mutation_strength, aging_intermediate_closure, &normalized_male_fertility_closure, &normalized_female_fertility_closure, tradeoff, start_b, remove_non_reproducing, male_menopause, female_menopause)
//     }

//     // println!("#########################################");
//     // println!("###### Simulation with low Gmax ######");
//     // println!("#########################################");
//     // initial_gmax_distribution = [0.055, 0.0];
//     // growth_name_part = "low_growth";

//     // let output_file_name = format!("./simulation_results/gmax_simulations/{}_{}_{}_{}_{}_{}.csv", base_name_part, mating_name_part, growth_name_part, removal_name_part, tradeoff_name_part, initial_gmax_distribution[0]);
//     // let mut wtr = Writer::from_path(output_file_name).unwrap();

//     // for i in 0..replicate_number{
//     //     println!("Replicate : {}/{}", i+1, replicate_number);
//     //     run_simulation(&mut wtr, population_cap, simulation_time, i, assortative_mating, &aging_parameters, &learning_parameters, &growth_parameters, initial_age_distribution, initial_b_distribution, initial_lmax_distribution, initial_gmax_distribution, initial_female_proportion, time_step, mutable_b, mutable_lmax, mutable_gmax, b_mutation_rate, lmax_mutation_rate, gmax_mutation_rate, b_mutation_strength, lmax_mutation_strength, gmax_mutation_strength, aging_intermediate_closure, &normalized_male_fertility_closure, &normalized_female_fertility_closure, tradeoff, start_b, remove_non_reproducing, male_menopause, female_menopause)
//     // }
}