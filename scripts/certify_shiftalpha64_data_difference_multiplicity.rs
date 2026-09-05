//! Certify the maximum multiplicity of pairwise differences among data coefficients.

use std::env;
use std::time::Instant;

const REDUCTION: u64 = 0x1b;
const DEFAULT_BLOCKS: usize = 16_384;

#[inline(always)]
fn multiply_x(value: u64) -> u64 {
    (value << 1) ^ (REDUCTION & 0u64.wrapping_sub(value >> 63))
}

fn main() {
    let blocks = env::args()
        .nth(1)
        .map(|value| value.parse::<usize>().expect("invalid block count"))
        .unwrap_or(DEFAULT_BLOCKS);
    assert!((2..=DEFAULT_BLOCKS).contains(&blocks));

    let mut alpha_zero = 1u64;
    for _ in 0..64 {
        alpha_zero = multiply_x(alpha_zero);
    }
    let mut coefficients = vec![0u64; blocks];
    coefficients[0] = alpha_zero;
    for index in 1..blocks {
        coefficients[index] = multiply_x(coefficients[index - 1]);
    }
    assert!(coefficients.windows(2).all(|pair| pair[0] != pair[1]));

    let pair_count = blocks * (blocks - 1) / 2;
    let started = Instant::now();
    let mut differences = Vec::with_capacity(pair_count);
    for second in 1..blocks {
        let right = coefficients[second];
        for first in 0..second {
            differences.push(coefficients[first] ^ right);
        }
    }
    assert_eq!(differences.len(), pair_count);
    assert!(differences.iter().all(|difference| *difference != 0));
    let generation_seconds = started.elapsed().as_secs_f64();

    let sort_started = Instant::now();
    differences.sort_unstable();
    let sort_seconds = sort_started.elapsed().as_secs_f64();
    let mut distinct = 0u64;
    let mut repeated_records = 0u64;
    let mut maximum_multiplicity = 0u64;
    let mut runs_with_repetitions = 0u64;
    let mut index = 0usize;
    while index < differences.len() {
        let mut end = index + 1;
        while end < differences.len() && differences[end] == differences[index] {
            end += 1;
        }
        let multiplicity = (end - index) as u64;
        distinct += 1;
        maximum_multiplicity = maximum_multiplicity.max(multiplicity);
        if multiplicity > 1 {
            runs_with_repetitions += 1;
            repeated_records += multiplicity - 1;
        }
        index = end;
    }
    let elapsed = started.elapsed().as_secs_f64();

    println!("{{");
    println!("  \"schema\": \"shiftalpha64-data-difference-multiplicity-v1\",");
    println!("  \"data_blocks\": {},", blocks);
    println!("  \"unordered_data_pairs\": {},", pair_count);
    println!("  \"distinct_nonzero_differences\": {},", distinct);
    println!("  \"maximum_difference_multiplicity\": {},", maximum_multiplicity);
    println!("  \"runs_with_repetitions\": {},", runs_with_repetitions);
    println!("  \"repeated_records\": {},", repeated_records);
    println!("  \"generation_seconds\": {:.12},", generation_seconds);
    println!("  \"sort_seconds\": {:.12},", sort_seconds);
    println!("  \"elapsed_seconds\": {:.12},", elapsed);
    println!("  \"validation\": {{");
    println!("    \"data_coefficients_are_distinct\": true,");
    println!("    \"all_pairwise_differences_are_nonzero\": true,");
    println!("    \"sorted_run_mass_equals_pair_count\": true");
    println!("  }}");
    println!("}}");
}
