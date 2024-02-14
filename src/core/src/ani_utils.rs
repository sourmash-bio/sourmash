// Streamlined set of utils for containment --> ANI estimation
// Equations based off of: https://github.com/KoslickiLab/mutation-rate-ci-calculator
// Reference: https://doi.org/10.1101/2022.01.11.475870

use crate::Error;
use roots::{find_root_brent, SimpleConvergency};
use statrs::distribution::{ContinuousCDF, Normal};

fn exp_n_mutated(l: f64, k: f64, r1: f64) -> f64 {
    let q = r1_to_q(k, r1);
    l * q
}

fn var_n_mutated(l: f64, k: f64, r1: f64, q: Option<f64>) -> Result<f64, Error> {
    if r1 == 0.0 {
        return Ok(0.0);
    }

    let q = q.unwrap_or_else(|| r1_to_q(k, r1));

    let var_n = l * (1.0 - q) * (q * (2.0 * k + (2.0 / r1) - 1.0) - 2.0 * k)
        + k * (k - 1.0) * (1.0 - q).powi(2)
        + (2.0 * (1.0 - q) / (r1.powi(2))) * ((1.0 + (k - 1.0) * (1.0 - q)) * r1 - q);

    if var_n < 0.0 {
        Err(Error::ANIEstimationError {
            message: "varN is less than 0.0".into(),
        })
    } else {
        Ok(var_n)
    }
}

fn exp_n_mutated_squared(l: f64, k: f64, p: f64) -> Result<f64, Error> {
    let var_n = var_n_mutated(l, k, p, None)?;
    let exp_n_squared = exp_n_mutated(l, k, p).powi(2);
    Ok(var_n + exp_n_squared)
}

fn probit(p: f64) -> f64 {
    Normal::new(0.0, 1.0).unwrap().inverse_cdf(p)
}

fn r1_to_q(k: f64, r1: f64) -> f64 {
    1.0 - (1.0 - r1).powi(k as i32)
}

fn get_exp_probability_nothing_common(
    ani_estimate: f64,
    ksize: f64,
    f_scaled: f64,
    n_unique_kmers: f64,
) -> Result<f64, Error> {
    if ani_estimate == 0.0 || ani_estimate == 1.0 {
        Ok(1.0 - ani_estimate)
    } else {
        // Calculate the expected log probability.
        let exp_nmut = exp_n_mutated(n_unique_kmers, ksize, 1.0 - ani_estimate);
        let mut expected_log_probability = (n_unique_kmers - exp_nmut) * (1.0 - f_scaled).ln();

        if expected_log_probability.is_infinite() {
            expected_log_probability = f64::NEG_INFINITY;
        }
        // Return the exponential of the expected log probability.
        Ok(expected_log_probability.exp())
    }
}

/// Streamlined function for ANI from containment.
/// todo: report ANI as % in 5.0?
pub fn ani_from_containment(containment: f64, ksize: f64) -> f64 {
    if containment == 0.0 {
        1.0
    } else if containment == 1.0 {
        0.0
    } else {
        1.0 - (1.0 - containment.powf(1.0 / ksize))
    }
}

// Calculate containment to ANI with confidence intervals
pub fn ani_from_containment_ci(
    containment: f64,
    ksize: f64,
    scaled: u64,
    n_unique_kmers: u64,
    confidence: Option<f64>,
    prob_threshold: Option<f64>,
) -> Result<(f64, f64, f64, f64), Error> {
    let confidence = confidence.unwrap_or(0.95);
    let prob_threshold = prob_threshold.unwrap_or(1e-3);

    let point_estimate = ani_from_containment(containment, ksize);

    // conversions needed throughout
    let scaled_f64 = scaled as f64;
    let f_scaled = 1.0 / scaled_f64;
    let n_unique_kmers = n_unique_kmers as f64;
    let alpha = 1.0 - confidence;

    let z_alpha = probit(1.0 - alpha / 2.0);
    let bias_factor = 1.0 - (1.0 - f_scaled).powi(n_unique_kmers as i32);
    let term_1 = (1.0 - f_scaled) / (f_scaled * (n_unique_kmers).powi(3) * bias_factor.powi(2));
    let term_2 = |pest: f64| {
        n_unique_kmers * exp_n_mutated(n_unique_kmers, ksize, pest)
            - exp_n_mutated_squared(n_unique_kmers, ksize, pest).unwrap_or(0.0)
    };
    let term_3 = |pest: f64| {
        var_n_mutated(n_unique_kmers, ksize, pest, None).unwrap_or(0.0)
            / (n_unique_kmers as f64).powi(2)
    };

    let var_direct = |pest: f64| term_1 * term_2(pest) + term_3(pest);

    let f1 = |pest: f64| {
        (1.0 - pest).powi(ksize as i32) + z_alpha * var_direct(pest).sqrt() - containment
    };
    let f2 = |pest: f64| {
        (1.0 - pest).powi(ksize as i32) - z_alpha * var_direct(pest).sqrt() - containment
    };

    let mut convergency = SimpleConvergency {
        eps: 1e-15,
        max_iter: 1000,
    };
    // check this:: what is the default? 0?
    let sol1 = find_root_brent(0.0000001, 0.9999999, &f1, &mut convergency).unwrap_or_default();
    let sol2 = find_root_brent(0.0000001, 0.9999999, &f2, &mut convergency).unwrap_or_default();

    // Actually, we really only want to calculate this if ANI is 0, since in all other cases
    // there's already something in common!
    // .. meaning it's not worth calculating for prefetch or gather results
    // do we want to create a separate version of this function excluding this, or just keep it?
    let prob_nothing_in_common =
        get_exp_probability_nothing_common(n_unique_kmers, ksize, point_estimate, f_scaled)?;

    Ok((point_estimate, sol1, sol2, prob_nothing_in_common))
}
