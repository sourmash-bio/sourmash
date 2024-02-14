use roots::{find_root_brent, SimpleConvergency};
use statrs::distribution::{ContinuousCDF, Normal};
use std::error::Error;

#[derive(Debug)]
#[allow(dead_code)]
pub struct CiAniResult {
    pub point_estimate: f64,
    prob_nothing_in_common: f64,
    dist_low: Option<f64>,
    dist_high: Option<f64>,
    p_threshold: f64,
}

fn exp_n_mutated(l: u64, k: u32, r1: f64) -> f64 {
    let q = r1_to_q(k, r1);
    l as f64 * q
}

fn var_n_mutated(l: u64, k: u32, r1: f64, q: Option<f64>) -> Result<f64, Box<dyn Error>> {
    if r1 == 0.0 {
        return Ok(0.0);
    }

    let q = q.unwrap_or_else(|| r1_to_q(k, r1));

    let var_n = l as f64 * (1.0 - q) * (q * (2.0 * k as f64 + (2.0 / r1) - 1.0) - 2.0 * k as f64)
        + k as f64 * (k as f64 - 1.0) * (1.0 - q).powi(2)
        + (2.0 * (1.0 - q) / (r1.powi(2))) * ((1.0 + (k as f64 - 1.0) * (1.0 - q)) * r1 - q);

    if var_n < 0.0 {
        Err("Error: varN < 0.0!".into())
    } else {
        Ok(var_n)
    }
}

fn exp_n_mutated_squared(l: u64, k: u32, p: f64) -> Result<f64, Box<dyn Error>> {
    let var_n = var_n_mutated(l, k, p, None)?;
    Ok(var_n + exp_n_mutated(l, k, p).powi(2))
}

fn probit(p: f64) -> f64 {
    Normal::new(0.0, 1.0).unwrap().inverse_cdf(p)
}

fn r1_to_q(k: u32, r1: f64) -> f64 {
    1.0 - (1.0 - r1).powi(k as i32)
}

fn get_exp_probability_nothing_common(
    mutation_rate: f64,
    kmer_size: u32,
    scaled: u64,
    n_unique_kmers: u64,
) -> f64 {
    // Inverse of the scale factor, used in probability calculation.
    let inverse_scaled = 1.0 / scaled as f64;

    // Handle special cases for mutation rate.
    if mutation_rate == 1.0 {
        1.0
    } else if mutation_rate == 0.0 {
        0.0
    } else {
        // Calculate the expected log probability.
        let expected_log_probability =
            get_expected_log_probability(n_unique_kmers, kmer_size, mutation_rate, inverse_scaled);
        // Return the exponential of the expected log probability.
        expected_log_probability.exp()
    }
}

fn get_expected_log_probability(
    n_unique_kmers: u64,
    ksize: u32,
    mutation_rate: f64,
    scaled_fraction: f64,
) -> f64 {
    let exp_nmut = exp_n_mutated(n_unique_kmers, ksize, mutation_rate);
    let result = (n_unique_kmers as f64 - exp_nmut) * (1.0 - scaled_fraction).ln();

    if result.is_infinite() {
        f64::NEG_INFINITY
    } else {
        result
    }
}

fn handle_seqlen_nkmers(
    ksize: u32,
    sequence_len_bp: Option<u64>,
    n_unique_kmers: Option<u64>,
) -> Result<u64, Box<dyn Error>> {
    match n_unique_kmers {
        Some(n) => Ok(n),
        None => match sequence_len_bp {
            Some(len) => Ok(len.saturating_sub(ksize as u64 - 1)),
            None => Err(
                "Error: distance estimation requires 'sequence_len_bp' or 'n_unique_kmers'".into(),
            ),
        },
    }
}

fn check_distance(dist: f64) -> Result<f64, Box<dyn Error>> {
    if dist >= 0.0 && dist <= 1.0 {
        Ok(dist)
    } else {
        Err(format!("Error: distance value {:.4} is not between 0 and 1!", dist).into())
    }
}

impl CiAniResult {
    fn new(
        point_estimate: f64,
        prob_nothing_in_common: f64,
        dist_low: Option<f64>,
        dist_high: Option<f64>,
        p_threshold: f64,
    ) -> Result<Self, Box<dyn Error>> {
        let dist_low_checked = dist_low.map_or(Ok(None), |d| check_distance(d).map(Some))?;
        let dist_high_checked = dist_high.map_or(Ok(None), |d| check_distance(d).map(Some))?;

        Ok(Self {
            point_estimate,
            prob_nothing_in_common,
            dist_low: dist_low_checked,
            dist_high: dist_high_checked,
            p_threshold,
        })
    }
}

pub fn containment_to_distance(
    containment: f64,
    ksize: u32,
    scaled: u64,
    n_unique_kmers: Option<u64>,
    sequence_len_bp: Option<u64>,
    confidence: Option<f64>,
    estimate_ci: Option<bool>,
    prob_threshold: Option<f64>,
) -> Result<CiAniResult, Box<dyn Error>> {
    let scaled_f64 = scaled as f64;
    let n_unique_kmers = handle_seqlen_nkmers(ksize, sequence_len_bp, n_unique_kmers)?;

    let confidence = confidence.unwrap_or(0.95);
    let estimate_ci = estimate_ci.unwrap_or(false);
    let prob_threshold = prob_threshold.unwrap_or(1e-3);

    let point_estimate = if containment == 0.0 {
        1.0
    } else if containment == 1.0 {
        0.0
    } else {
        1.0 - containment.powf(1.0 / ksize as f64)
    };

    let mut sol1 = None;
    let mut sol2 = None;

    if estimate_ci {
        let alpha = 1.0 - confidence;
        let z_alpha = probit(1.0 - alpha / 2.0);
        let f_scaled = 1.0 / scaled_f64;
        let bias_factor = 1.0 - (1.0 - f_scaled).powi(n_unique_kmers as i32);
        let term_1 =
            (1.0 - f_scaled) / (f_scaled * (n_unique_kmers as f64).powi(3) * bias_factor.powi(2));
        let term_2 = |pest: f64| {
            n_unique_kmers as f64 * exp_n_mutated(n_unique_kmers, ksize, pest)
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
        sol1 = find_root_brent(0.0000001, 0.9999999, &f1, &mut convergency).ok();
        sol2 = find_root_brent(0.0000001, 0.9999999, &f2, &mut convergency).ok();
    }

    let prob_nothing_in_common =
        get_exp_probability_nothing_common(point_estimate, ksize, scaled, n_unique_kmers);

    CiAniResult::new(
        point_estimate,
        prob_nothing_in_common,
        sol1,
        sol2,
        prob_threshold,
    )
}

/* I calculate the distance, then I do 1 - distance to get the ANI identity.

fn main() {
    let containment = 0.1;
    let ksize = 31;
    let scaled = 1000; // u64
    let n_unique_kmers = 1000;
    let sequence_len_bp = n_unique_kmers * scaled;
    let confidence = Some(0.95);
    let estimate_ci = Some(true);
    let prob_threshold = Some(1e-3);

    let result = containment_to_distance(
        containment,
        ksize,
        scaled,
        Some(n_unique_kmers),
        Some(sequence_len_bp),
        confidence,
        estimate_ci,
        prob_threshold,
    );


    match result {
        Ok(ci_result) => {
            let ani_result = 1.0 - ci_result.point_estimate;
            println!("ANI: {}", ani_result);
        }

        Err(e) => println!("Error occurred: {:?}", e),
    }
}
*/
