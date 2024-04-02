pub fn _aging_gompertz(x: f64, aging_parameters: &[f64]) -> f64 {
    let (a, b) = (aging_parameters[0], aging_parameters[1]);
    a * (x * b).exp()
}

pub fn aging_gompertz_makeham(x: f64, aging_parameters: &[f64]) -> f64 {
    let (a, b, c) = (
        aging_parameters[0],
        aging_parameters[1],
        aging_parameters[2],
    );
    c + a * (x * b).exp()
}

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
