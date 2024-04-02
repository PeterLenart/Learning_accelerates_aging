mod gla_package;
use crate::gla_package::{
    gla::{aging_gompertz_makeham, gla_model, growth_function, learning_function},
    simulate::run_simulation,
};
use csv::Writer;

// use easybench::bench;

// use peroxide::fuga::{GaussLegendre, G7K15R};
// use peroxide::numerical::integral::{gauss_kronrod_quadrature, integrate};

fn main() {
    let time_step = 1.0;
    let initial_female_proportion = 0.5;
    let minimum_mortality = 1e-5;
    let aging_parameters = [0.00275961297460256,0.04326224872667336,0.025201676835511704] ;
    let learning_parameters = [0.01606792505529796,39.006865144958745,0.11060749334680318];
    let growth_parameters: [f64; 2] = [0.05168141300917714,0.08765165352033985];

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

    let initial_age_distribution = [20.0, 10.0];
    let initial_b_distribution = [0.14, 0.005];
    let mut initial_lmax_distribution = [2.0, 0.0];

    let population_cap = 10000;
    let simulation_time : usize = 10000;
    let replicate_number = 100;
    let assortative_mating = false;

    let mutable_b = true;
    let mutable_lmax = false;
    let b_mutation_rate: f64 = 0.02;
    let lmax_mutation_rate: f64 = 0.02;
    let b_mutation_strength = 0.012;
    let lmax_mutation_strength = 0.012;

    let base_name_part = "asexual_test";
    let mut learning_name_part = "with_learning";
    let mut mating_name_part = "random_mating";

    if assortative_mating{
        mating_name_part = "assortative_mating";
    }

    println!("######################################");
    println!("###### Simulation with learning ######");
    println!("######################################");

    let output_file_name = format!("./simulation_results/{}_{}_{}_{}.csv", base_name_part, mating_name_part, learning_name_part,initial_lmax_distribution[0]);
    let mut wtr = Writer::from_path(output_file_name).unwrap();


    for i in 0..replicate_number {
        println!("Replicate : {}/{}", i + 1, replicate_number);
        run_simulation(
            &mut wtr,
            population_cap,
            simulation_time,
            i,
            &aging_parameters,
            &learning_parameters,
            &growth_parameters,
            initial_age_distribution,
            initial_b_distribution,
            initial_lmax_distribution,
            time_step,
            mutable_b,
            mutable_lmax,
            b_mutation_rate,
            lmax_mutation_rate,
            b_mutation_strength,
            lmax_mutation_strength,
            aging_intermediate_closure,
        )
    }

    // println!("#########################################");
    // println!("###### Simulation without learning ######");
    // println!("#########################################");
    // initial_lmax_distribution = [0.0, 0.0];
    // learning_name_part = "no_learning";

    // let output_file_name = format!("./simulation_results/{}_{}_{}.csv", base_name_part, mating_name_part, learning_name_part);
    // let mut wtr = Writer::from_path(output_file_name).unwrap();


    // for i in 0..replicate_number {
    //     println!("Replicate : {}/{}", i + 1, replicate_number);
    //     run_simulation(
    //         &mut wtr,
    //         population_cap,
    //         simulation_time,
    //         i,
    //         &aging_parameters,
    //         &learning_parameters,
    //         &growth_parameters,
    //         initial_age_distribution,
    //         initial_b_distribution,
    //         initial_lmax_distribution,
    //         time_step,
    //         mutable_b,
    //         mutable_lmax,
    //         b_mutation_rate,
    //         lmax_mutation_rate,
    //         b_mutation_strength,
    //         lmax_mutation_strength,
    //         aging_intermediate_closure,
    //     )
    // }
}
