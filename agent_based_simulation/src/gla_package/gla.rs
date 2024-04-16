use ndarray::{Array, ArrayView1};
use optimize::{Minimizer, NelderMeadBuilder};

/// Calculate mortality based on the Gompertz model. 
///
/// # Arguments
/// * `x` - The age at which to calculate the mortality rate.
/// * `aging_parameters` - An array containing the parameters for the Gompertz function:
///     - `aging_parameters[0]` (a): The initial mortality rate (intercept).
///     - `aging_parameters[1]` (b): The rate of increase in mortality rate with age (slope).
///
/// # Returns
/// Returns the mortality rate at age `x`, calculated using the provided `a` and `b` parameters.
///
/// # Examples
/// ```
/// let age = 50.0;
/// let parameters = [0.01, 0.02];
/// let mortality_rate = _aging_gompertz(age, &parameters);
/// assert_eq!(mortality_rate, 0.01 * (50.0 * 0.02).exp());
/// ```
pub fn _aging_gompertz(x: f64, aging_parameters: &[f64]) -> f64 {
    let (a, b) = (aging_parameters[0], aging_parameters[1]);
    a * (x * b).exp()
}

/// Calculate mortality based on the Gompertz-Makeham model.
///
/// # Arguments
/// * `x` - The age at which to calculate the mortality rate.
/// * `aging_parameters` - An array containing the parameters for the Gompertz-Makeham function:
///     - `aging_parameters[0]` (a): The initial mortality rate (intercept).
///     - `aging_parameters[1]` (b): The rate of increase in mortality rate with age (slope).
///     - `aging_parameters[2]` (c): Constant age-independent mortality factor.
///
/// # Returns
/// Returns the mortality rate at age `x`, calculated using the provided `a`, `b`, and `c` parameters.
///
/// # Examples
/// ```
/// let age = 50.0;
/// let parameters = [0.01, 0.02, 0.005];
/// let mortality_rate = aging_gompertz_makeham(age, &parameters);
/// assert_eq!(mortality_rate, 0.005 + 0.01 * (50.0 * 0.02).exp());
/// ```
pub fn aging_gompertz_makeham(x: f64, aging_parameters: &[f64]) -> f64 {
    let (a, b, c) = (
        aging_parameters[0],
        aging_parameters[1],
        aging_parameters[2],
    );
    c + a * (x * b).exp()
}

/// Calculate benefit of learning on mortality based on a modified logistic model.
///
/// # Arguments
/// * `x` - The age at which to calculate the learning benefit.
/// * `learning_parameters` - An array containing the parameters for the learning curve:
///     - `learning_parameters[0]` (lmax): The maximum learning benefit.
///     - `learning_parameters[1]` (k): The inflection point of the curve.
///     - `learning_parameters[2]` (n): The steepness of the curve.
///
/// # Returns
/// Returns the benefit of learning on mortality at the input `x`, calculated using the provided `lmax`, `k`, and `n`.
///
pub fn learning_function(x: f64, learning_parameters: &[f64]) -> f64 {
    let (lmax, k, n) = (
        learning_parameters[0],
        learning_parameters[1],
        learning_parameters[2],
    );
    lmax * ((1_f64 / (1_f64 + (n * (x - k)).exp())) - 1_f64)
}

pub fn growth_function(x: f64, growth_parameters: &[f64]) -> f64 {
    let (gmax, growth_rate) = (growth_parameters[0], growth_parameters[1]);
    gmax * ((1_f64 / (1_f64 + x.powf(growth_rate))) - 1_f64)
}

pub fn fertility_brass_polynomial(x: f64, fertility_parameters: &[f64]) -> f64 {
    let (c, d, w) = (
        fertility_parameters[0],
        fertility_parameters[1],
        fertility_parameters[2],
    );
    if (x > d) && (x < (d + w)) {
        c * (x - d) * ((d + w - x).powi(2))
    } else {
        0_f64
    }
}

pub fn constant_fertility(_x: f64, fertility_parameters: &[f64]) -> f64 {
    let c = fertility_parameters[0];
    c
}



pub fn gla_model<T>(
    x: f64,
    aging_func: T,
    learning_func: T,
    growth_func: T,
    aging_parameters: &[f64],
    learning_parameters: &[f64],
    growth_parameters: &[f64],
    minimum_mortality: f64,
) -> f64
where
    T: Fn(f64, &[f64]) -> f64,
{
    let aging_result = aging_func(x, aging_parameters);
    let learning_result = learning_func(x, learning_parameters);
    let growth_result = growth_func(x, growth_parameters);

    let gla_result = aging_result + learning_result + growth_result;

    if gla_result < minimum_mortality {
        minimum_mortality
    } else {
        gla_result
    }
}

pub fn find_maximum_fertility<T>(
    fertility_function: &T,
    fertility_parameters: &[f64],
    first_guess: f64,
) -> f64
where
    T: Fn(f64, &[f64]) -> f64,
{
    let fertility_cost_function =
        |x: ArrayView1<f64>| -fertility_function(x[0], fertility_parameters);

    // Create a minimizer using the builder pattern. If some of the parameters are not given, default values are used.
    let minimizer = NelderMeadBuilder::default()
        .xtol(1e-6f64)
        .ftol(1e-6f64)
        .maxiter(50000)
        .build()
        .unwrap();

    // Set the starting guess
    let args = Array::from_vec(vec![first_guess]);

    // Run the optimization
    let ans = minimizer.minimize(&fertility_cost_function, args.view());

    fertility_function(ans[0], fertility_parameters)
}
