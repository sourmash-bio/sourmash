// Streamlined set of utils for containment --> ANI estimation
// Equations based off of: https://github.com/KoslickiLab/mutation-rate-ci-calculator
// Reference: https://doi.org/10.1101/2022.01.11.475870

use roots::{find_root_brent, SimpleConvergency};
use statrs::distribution::{ContinuousCDF, Normal};

use crate::{Error, ScaledType};

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

// prior versions of ani estimation also returned this value by default.
// BUT, it's not really something we need to calculate for prefetch/gather,
// since there will always be something in common (or comparison would not be happening)
// I'm not sure where the right place to put this back in is.. perhaps add minhash.contained_ani()
// and return this there.
// Usage:
// let prob_nothing_in_common =
//     get_exp_probability_nothing_common(n_unique_kmers, ksize, point_estimate, f_scaled)?;
// let prob_threshold = prob_threshold.unwrap_or(1e-3);
#[allow(dead_code)]
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
        0.0
    } else if containment == 1.0 {
        1.0
    } else {
        1.0 - (1.0 - containment.powf(1.0 / ksize))
    }
}

// Calculate containment to ANI with confidence intervals
pub fn ani_ci_from_containment(
    containment: f64,
    ksize: f64,
    scaled: ScaledType,
    n_unique_kmers: u64,
    confidence: Option<f64>,
) -> Result<(f64, f64), Error> {
    if containment == 0.0 {
        return Ok((0.0, 0.0));
    } else if containment == 1.0 {
        return Ok((1.0, 1.0));
    }
    let confidence = confidence.unwrap_or(0.95);

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
        var_n_mutated(n_unique_kmers, ksize, pest, None).unwrap_or(0.0) / (n_unique_kmers).powi(2)
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

    let dist_sol1 =
        find_root_brent(0.0000001, 0.9999999, &f1, &mut convergency).unwrap_or_default();
    let dist_sol2 =
        find_root_brent(0.0000001, 0.9999999, &f2, &mut convergency).unwrap_or_default();

    Ok((1.0 - dist_sol1, 1.0 - dist_sol2))
}

#[cfg(test)]
mod tests {

    use super::*;
    use std::f64::EPSILON;

    #[test]
    fn test_containment_to_ani_zero() {
        let contain = 0.0;
        let ksize = 21;
        let scaled = 10;
        let n_unique_kmers = 100;
        let confidence = Some(0.95);
        let res = ani_from_containment(contain, ksize as f64);
        assert_eq!(res, 0.0);
        let (ci_low, ci_high) =
            ani_ci_from_containment(contain, ksize as f64, scaled, n_unique_kmers, confidence)
                .unwrap();

        eprintln!("{}", ci_low);
        eprintln!("{}", ci_high);
        assert_eq!(ci_low, 0.0);
        assert_eq!(ci_high, 0.0);
    }

    #[test]
    fn test_containment_to_ani_one() {
        let contain = 1.0;
        let ksize = 21;
        let scaled = 10;
        let n_unique_kmers = 100;
        let confidence = None;
        let res = ani_from_containment(contain, ksize as f64);
        assert_eq!(res, 1.0);
        let (ci_low, ci_high) =
            ani_ci_from_containment(contain, ksize as f64, scaled, n_unique_kmers, confidence)
                .unwrap();
        assert_eq!(ci_low, 1.0);
        assert_eq!(ci_high, 1.0);
    }

    #[test]
    fn test_containment_to_ani_scaled1() {
        let contain = 0.5;
        let ksize = 21;
        let scaled = 1;
        let n_unique_kmers = 10000;
        let confidence = None;
        let ani = ani_from_containment(contain, ksize as f64);
        assert!((ani - 0.9675317785238916) < EPSILON);
        let (ci_low, ci_high) =
            ani_ci_from_containment(contain, ksize as f64, scaled, n_unique_kmers, confidence)
                .unwrap();
        assert!((ci_low - 0.9635213980271021) < EPSILON);
        assert!((ci_high - 0.9712900870335944) < EPSILON);
    }

    #[test]
    fn test_containment_to_ani_scaled100() {
        let contain = 0.1;
        let ksize = 31;
        let scaled = 100;
        let n_unique_kmers = 10000;
        let confidence = None;
        let ani = ani_from_containment(contain, ksize as f64);
        assert!((ani - 0.9284145445194744) < EPSILON);
        let (ci_low, ci_high) =
            ani_ci_from_containment(contain, ksize as f64, scaled, n_unique_kmers, confidence)
                .unwrap();
        assert!((ci_low - 0.9094445232754665) < EPSILON);
        assert!((ci_high - 0.9467922076143345) < EPSILON);
    }

    #[test]
    fn test_containment_to_ani_scaled100_2() {
        let contain = 0.5;
        let ksize = 21;
        let scaled = 100;
        let n_unique_kmers = 10000;
        let confidence = None;
        let ani = ani_from_containment(contain, ksize as f64);
        assert!((ani - 0.9675317785238916) < EPSILON);
        let (ci_low, ci_high) =
            ani_ci_from_containment(contain, ksize as f64, scaled, n_unique_kmers, confidence)
                .unwrap();
        assert!((ci_low - 0.9569003945603415) < EPSILON);
        assert!((ci_high - 0.9762879360833708) < EPSILON);
    }

    #[test]
    fn test_var_n_mutated_zero() {
        let r = 0.0;
        let ksize = 31;
        let nkmers = 200;
        let var_n_mut = var_n_mutated(nkmers as f64, ksize as f64, r, None).unwrap(); // Assuming the function returns a Result
        assert_eq!(var_n_mut, 0.0, "Expected variance to be 0 for r=0");
    }

    #[test]
    fn test_var_n_mutated_value_error() {
        let r = 10.0;
        let ksize = 31;
        let nkmers = 200;
        match var_n_mutated(nkmers as f64, ksize as f64, r, None) {
            Err(e) => assert_eq!(
                e.to_string(),
                "error while calculating ANI confidence intervals: varN is less than 0.0",
                "Unexpected error message"
            ),
            Ok(_) => panic!("Expected an error, but got Ok"),
        }
    }

    #[test]
    fn test_var_n_mutated_success() {
        let r = 0.4;
        let ksize = 31;
        let nkmers = 200_000;
        let var_n_mut = var_n_mutated(nkmers as f64, ksize as f64, r, None).unwrap(); // Assuming the function returns a Result
        let expected = 0.10611425440741508;
        assert!(
            (var_n_mut - expected).abs() < f64::EPSILON,
            "Variance did not match expected value"
        );
    }

    #[test]
    fn test_r1_to_q() {
        let k = 2.0;
        let r1 = 0.5;
        let result = r1_to_q(k, r1);
        let expected = 0.75;

        assert!(
            (result - expected).abs() < EPSILON,
            "The result of r1_to_q({}, {}) was {}, but {} was expected",
            k,
            r1,
            result,
            expected
        );
    }

    #[test]
    fn test_exp_n_mutated() {
        let l = 100.0;
        let k = 2.0;
        let r1 = 0.5;

        // Calculate the expected result based on the inputs
        let expected_q = r1_to_q(k, r1);
        let expected_result = l * expected_q;

        let result = exp_n_mutated(l, k, r1);

        assert!(
            (result - expected_result).abs() < EPSILON,
            "The result of exp_n_mutated({}, {}, {}) was {}, but {} was expected",
            l,
            k,
            r1,
            result,
            expected_result
        );
    }

    #[test]
    fn test_get_exp_probability_nothing_common_ani_zero() {
        let ani_estimate = 0.0;
        let ksize = 31.0;
        let f_scaled = 0.1;
        let n_unique_kmers = 1000.0;

        let result =
            get_exp_probability_nothing_common(ani_estimate, ksize, f_scaled, n_unique_kmers)
                .unwrap();
        assert_eq!(
            result, 1.0,
            "Expected probability for ani_estimate of 0 to be 1.0"
        );
    }

    #[test]
    fn test_get_exp_probability_nothing_common_ani_one() {
        let ani_estimate = 1.0;
        let ksize = 31.0;
        let f_scaled = 0.1;
        let n_unique_kmers = 1000.0;

        let result =
            get_exp_probability_nothing_common(ani_estimate, ksize, f_scaled, n_unique_kmers)
                .unwrap();
        assert_eq!(
            result, 0.0,
            "Expected probability for ani_estimate of 1 to be 0.0"
        );
    }

    #[test]
    fn test_get_exp_probability_nothing_common() {
        let contain = 0.1;
        let ksize = 31 as f64;
        let scaled = 10;
        let f_scaled = 1.0 / scaled as f64;
        let n_unique_kmers = 1000;

        let ani = ani_from_containment(contain, ksize);
        let result =
            get_exp_probability_nothing_common(ani, ksize, f_scaled, n_unique_kmers as f64)
                .unwrap();
        assert!(
            result >= 0.0 && result <= 1.0,
            "The result should be a valid probability"
        );
        assert!((result - 0.000026561398887587855) < EPSILON);
    }
}
