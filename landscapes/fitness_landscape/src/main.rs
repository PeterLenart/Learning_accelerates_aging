mod gla_package;
use crate::gla_package::{gla::{
    aging_gompertz_makeham, fertility_brass_polynomial, constant_fertility, gla_model,
    growth_function, learning_function, _aging_gompertz, toy_model, mortality_improvement_function, step_fertility,
}, landscape::{generate_grids, compute_fitness_grid_parallel_lmax, compute_fitness_grid_parallel_gmax, compute_fitness_grid_point_lmax, compute_fitness_grid_point_gmax}};
use ndarray::{Array2, Array1, stack, Axis, s};
use std::error::Error;
use std::fs::File;
use csv::WriterBuilder;

fn main() -> Result<(), Box<dyn Error>>{

    // LMAX GRID

    // Generate two regular grids of b and lmax values, a base grid and a grid where b is reduced by some factor
    let dim = 500;
    
    // let (base_b_grid, reduced_b_grid) = generate_grids(
    //     "regular",
    //     1e-5, 
    //     0.14, 
    //     0.0, 
    //     0.20,
    //     dim,
    //     0.95,
    // );

    let (base_b_grid, reduced_b_grid) = generate_grids(
        "regular",
        1e-5, 
        0.1, 
        0.0, 
        1.0,
        dim,
        0.95,
    );

    // Define GLA parameters
    let minimum_mortality = 1e-5f64;
    // let aging_parameters = [0.00275961297460256,0.04326224872667336,0.025201676835511704] ;
    // let learning_parameters = [0.01606792505529796,39.006865144958745,0.11060749334680318];
    // // let learning_parameters = [0.01606792505529796, 35.0,0.11060749334680318];
    // let growth_parameters: [f64; 2] = [0.05168141300917714,0.08765165352033985];

    // Define toy model parameters
    let aging_parameters = [0.01, 0.02];
    // Hijacking the learning parameters as the improvement parameters for the toy model
    let learning_parameters = [0.5, 0.0, 10.0];
    // The growth parameters are not used in the toy model
    let growth_parameters = [1.0, 1.0];

    // Define intermediate aging closure
    // let aging_intermediate_closure = |x: f64,
    //                                   aging_parameters: &[f64],
    //                                   learning_parameters: &[f64],
    //                                   growth_parameters: &[f64]|
    //  -> f64 {
    //     gla_model(
    //         x,
    //         aging_gompertz_makeham as fn(f64, &[f64]) -> f64,
    //         learning_function,
    //         growth_function,
    //         &aging_parameters,
    //         &learning_parameters,
    //         &growth_parameters,
    //         minimum_mortality,
    //     )
    // };

    // Dirty hack to use toy model, the growth parameters are not used
    let aging_intermediate_closure = |x: f64,
                                      aging_parameters: &[f64],
                                      improvement_parameters: &[f64],
                                      growth_parameters: &[f64]|
     -> f64 {
        toy_model(
            x,
            _aging_gompertz as fn(f64, &[f64]) -> f64,
            mortality_improvement_function,
            &aging_parameters,
            &improvement_parameters,
            minimum_mortality,
        )
    };

    // Define fertility parameters
    // let fertility_parameters = [2.445e-5, 14.8, 32.836];
    // let fertility_function = fertility_brass_polynomial;

    // let fertility_parameters = [1f64];
    // let fertility_function = constant_fertility;

    let fertility_parameters = [1f64, 15.0, 100.0];
    let fertility_function = step_fertility;

    // let fertility_parameters = [0.1f64, 0.1];
    // let fertility_function = linear_fertility;

    // Define fertility closure
    let fertility_closure = |x: f64| -> f64 {
        fertility_function(x, &fertility_parameters)
    };

    // Compute the fitness landscape for the base b grid
    let fitness_base_b = compute_fitness_grid_parallel_lmax(&base_b_grid, &aging_intermediate_closure, &aging_parameters, &learning_parameters, &growth_parameters, &fertility_closure);
    
    // Save result
    // let output_file_name = "./output/toy_model_landscapes/csv/fitness_40_50.csv";
    // Name the output file based on the parameters
    let output_file_name = format!("./output/toy_model_landscapes/csv/fitness_{}_{}.csv", learning_parameters[1], learning_parameters[2]);

    // Combine the grid and fitness difference into one array
    let fitness_base_b_reshaped = fitness_base_b.clone().insert_axis(Axis(1));
    let combined_data = stack![Axis(1), base_b_grid.clone(), fitness_base_b_reshaped];

    // Write to CSV
    let file = File::create(output_file_name)?;
    let mut writer = WriterBuilder::new().has_headers(true).from_writer(file);

    // Write headers
    let _ = writer.write_record(&["b", "lmax", "fitness"]);

    // Write data
    for row in combined_data.genrows() {
        let _ = writer.serialize(row.to_vec());
    }

    let _ = writer.flush();

    // Compute the fitness landscape for the reduced b grid
    let fitness_reduced_b = compute_fitness_grid_parallel_lmax(&reduced_b_grid, &aging_intermediate_closure, &aging_parameters, &learning_parameters, &growth_parameters, &fertility_closure);

    // Compute the difference in fitness between the base and reduced b grids
    let fitness_difference = fitness_reduced_b - fitness_base_b;

    // Save result
    // let output_file_name = "./output/toy_model_landscapes/csv/fitness_difference_40_50.csv";
    // Name the output file based on the parameters
    let output_file_name = format!("./output/toy_model_landscapes/csv/fitness_difference_{}_{}.csv", learning_parameters[1], learning_parameters[2]);

    // Combine the grid and fitness difference into one array
    let fitness_difference_reshaped = fitness_difference.insert_axis(Axis(1));
     let combined_data = stack![Axis(1), base_b_grid, fitness_difference_reshaped];

    // Write to CSV
    let file = File::create(output_file_name)?;
    let mut writer = WriterBuilder::new().has_headers(true).from_writer(file);

    // Write headers
    let _ = writer.write_record(&["b", "lmax", "fitness_difference"]);

    // Write data
    for row in combined_data.genrows() {
        let _ = writer.serialize(row.to_vec());
    }

    let _ = writer.flush();

    Ok(())


    // // GMAX GRID

    // // Generate two regular grids of b and gmax values, a base grid and a grid where b is reduced by some factor
    // let dim = 500;
    
    // let (base_b_grid, reduced_b_grid) = generate_grids(
    //     "regular",
    //     1e-5, 
    //     0.14, 
    //     0.0, 
    //     0.20,
    //     dim,
    //     0.95,
    // );

    // // Define GLA parameters
    // let minimum_mortality = 1e-5f64;
    // let aging_parameters = [0.00275961297460256,0.04326224872667336,0.025201676835511704] ;
    // let learning_parameters = [0.01606792505529796,39.006865144958745,0.11060749334680318];
    // let growth_parameters: [f64; 2] = [0.05168141300917714,0.08765165352033985];

    // // Define intermediate aging closure
    // let aging_intermediate_closure = |x: f64,
    //                                   aging_parameters: &[f64],
    //                                   learning_parameters: &[f64],
    //                                   growth_parameters: &[f64]|
    //  -> f64 {
    //     gla_model(
    //         x,
    //         aging_gompertz_makeham as fn(f64, &[f64]) -> f64,
    //         learning_function,
    //         growth_function,
    //         &aging_parameters,
    //         &learning_parameters,
    //         &growth_parameters,
    //         minimum_mortality,
    //     )
    // };

    // // Define fertility parameters
    // let fertility_parameters = [2.445e-5, 14.8, 32.836];
    // let fertility_function = fertility_brass_polynomial;

    // // let fertility_parameters = [1f64];
    // // let fertility_function = constant_fertility;

    // // let fertility_parameters = [1f64, 0.02];
    // // let fertility_function = linear_fertility;

    // // Define fertility closure
    // let fertility_closure = |x: f64| -> f64 {
    //     fertility_function(x, &fertility_parameters)
    // };

    // // Compute the fitness landscape for the base b grid
    // let fitness_base_b = compute_fitness_grid_parallel_gmax(&base_b_grid, &aging_intermediate_closure, &aging_parameters, &learning_parameters, &growth_parameters, &fertility_closure);

    // // Save result
    // let output_file_name = "../output/gmax_fitness_landscapes/csv/gmax_fitness_brass_polynomial_g7k15.csv";

    // // Combine the grid and fitness difference into one array
    // let fitness_base_b_reshaped = fitness_base_b.clone().insert_axis(Axis(1));
    // let combined_data = stack![Axis(1), base_b_grid.clone(), fitness_base_b_reshaped];

    // // Write to CSV
    // let file = File::create(output_file_name)?;
    // let mut writer = WriterBuilder::new().has_headers(true).from_writer(file);

    // // Write headers
    // let _ = writer.write_record(&["b", "gmax", "fitness"]);

    // // Write data
    // for row in combined_data.genrows() {
    //     let _ = writer.serialize(row.to_vec());
    // }

    // let _ = writer.flush();

    // // Compute the fitness landscape for the reduced b grid
    // let fitness_reduced_b = compute_fitness_grid_parallel_gmax(&reduced_b_grid, &aging_intermediate_closure, &aging_parameters, &learning_parameters, &growth_parameters, &fertility_closure);

    // // Compute the difference in fitness between the base and reduced b grids
    // let fitness_difference = fitness_reduced_b - fitness_base_b;

    // // Save result
    // let output_file_name = "../output/gmax_fitness_landscapes/csv/gmax_fitness_difference_brass_polynomial_g7k15.csv";

    // // Combine the grid and fitness difference into one array
    // let fitness_difference_reshaped = fitness_difference.insert_axis(Axis(1));
    //  let combined_data = stack![Axis(1), base_b_grid, fitness_difference_reshaped];

    // // Write to CSV
    // let file = File::create(output_file_name)?;
    // let mut writer = WriterBuilder::new().has_headers(true).from_writer(file);

    // // Write headers
    // let _ = writer.write_record(&["b", "gmax", "fitness_difference"]);

    // // Write data
    // for row in combined_data.genrows() {
    //     let _ = writer.serialize(row.to_vec());
    // }

    // let _ = writer.flush();

    // Ok(())

}