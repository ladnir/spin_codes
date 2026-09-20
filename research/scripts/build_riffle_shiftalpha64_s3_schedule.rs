//! Build the exact S3-canonical finite-support schedule table.

use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

const REDUCTION: u64 = 0x1b;
const DEFAULT_BLOCKS: usize = 16_384;
const P0_TAG: u64 = 1u64 << 63;

#[derive(Clone, Copy)]
struct TaggedRecord {
    key: u64,
    tagged_multiplicity: u64,
}

#[target_feature(enable = "pclmulqdq")]
#[inline]
unsafe fn multiply(left: u64, right: u64) -> u64 {
    use std::arch::x86_64::{
        _mm_clmulepi64_si128, _mm_cvtsi128_si64, _mm_set_epi64x, _mm_srli_si128,
    };
    let a = _mm_set_epi64x(0, left as i64);
    let b = _mm_set_epi64x(0, right as i64);
    let product = _mm_clmulepi64_si128::<0x00>(a, b);
    let low = _mm_cvtsi128_si64(product) as u64;
    let high = _mm_cvtsi128_si64(_mm_srli_si128::<8>(product)) as u64;
    let overflow = (high >> 63) ^ (high >> 61) ^ (high >> 60);
    let mut reduced = low ^ high ^ (high << 1) ^ (high << 3) ^ (high << 4);
    reduced ^= overflow ^ (overflow << 1) ^ (overflow << 3) ^ (overflow << 4);
    reduced
}

#[inline(always)]
fn field_multiply(left: u64, right: u64) -> u64 {
    unsafe { multiply(left, right) }
}

#[inline(always)]
fn multiply_x(value: u64) -> u64 {
    (value << 1) ^ (REDUCTION & 0u64.wrapping_sub(value >> 63))
}

fn multiply_scalar(mut left: u64, mut right: u64) -> u64 {
    let mut result = 0u64;
    while right != 0 {
        if right & 1 != 0 {
            result ^= left;
        }
        right >>= 1;
        left = multiply_x(left);
    }
    result
}

fn power(mut base: u64, mut exponent: u64) -> u64 {
    let mut result = 1u64;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = field_multiply(result, base);
        }
        exponent >>= 1;
        if exponent != 0 {
            base = field_multiply(base, base);
        }
    }
    result
}

fn inverse(value: u64) -> u64 {
    assert_ne!(value, 0);
    power(value, u64::MAX - 1)
}

#[inline(always)]
fn invariant_from_factors(
    first: u64,
    second: u64,
    third: u64,
    inverse_first: u64,
    inverse_second: u64,
    inverse_third: u64,
) -> u64 {
    assert_eq!(first ^ second ^ third, 0);
    let first_squared = field_multiply(first, first);
    let second_squared = field_multiply(second, second);
    let numerator_base = first_squared
        ^ field_multiply(first, second)
        ^ second_squared;
    let numerator = field_multiply(
        field_multiply(numerator_base, numerator_base),
        numerator_base,
    );
    let inverse_product = field_multiply(
        field_multiply(inverse_first, inverse_second),
        inverse_third,
    );
    let inverse_denominator = field_multiply(inverse_product, inverse_product);
    field_multiply(numerator, inverse_denominator)
}

fn choose2(value: u64) -> u64 {
    value * (value - 1) / 2
}

fn choose3(value: u64) -> u64 {
    value * (value - 1) * (value - 2) / 6
}

fn main() {
    if !std::is_x86_feature_detected!("pclmulqdq") {
        panic!("pclmulqdq is required");
    }
    let output_path = env::args()
        .nth(1)
        .expect("usage: schedule-builder output.bin [blocks]");
    let blocks = env::args()
        .nth(2)
        .map(|value| value.parse::<usize>().expect("invalid block count"))
        .unwrap_or(DEFAULT_BLOCKS);
    assert!((3..=DEFAULT_BLOCKS).contains(&blocks));

    for probe in 0..1000u64 {
        let left = 0x0123456789abcdefu64
            .wrapping_add(probe.wrapping_mul(0x9e3779b97f4a7c15));
        let right = 0xfedcba9876543210u64
            .wrapping_add(probe.wrapping_mul(0xd1342543de82ef95));
        assert_eq!(field_multiply(left, right), multiply_scalar(left, right));
    }

    let started = Instant::now();
    let mut gamma = vec![1u64; blocks];
    let mut gamma_inverse = vec![1u64; blocks];
    let inverse_gamma = inverse(2);
    for index in 1..blocks {
        gamma[index] = multiply_x(gamma[index - 1]);
        gamma_inverse[index] = field_multiply(gamma_inverse[index - 1], inverse_gamma);
    }
    let mut inverse_one_plus = vec![0u64; blocks];
    for difference in 1..blocks {
        inverse_one_plus[difference] = inverse(1 ^ gamma[difference]);
    }

    let exact_entries = (blocks - 1) * (blocks - 2) / 2 + blocks - 1;
    let mut records = Vec::<TaggedRecord>::with_capacity(exact_entries);
    let mut data_mass = 0u64;
    for third_difference in 2..blocks {
        let third_coefficient = gamma[third_difference];
        let multiplicity = (blocks - third_difference) as u64;
        for second_difference in 1..third_difference {
            let second_coefficient = gamma[second_difference];
            let first = second_coefficient ^ third_coefficient;
            let second = 1 ^ third_coefficient;
            let third = 1 ^ second_coefficient;
            let inverse_first = field_multiply(
                gamma_inverse[second_difference],
                inverse_one_plus[third_difference - second_difference],
            );
            let key = invariant_from_factors(
                first,
                second,
                third,
                inverse_first,
                inverse_one_plus[third_difference],
                inverse_one_plus[second_difference],
            );
            records.push(TaggedRecord {
                key,
                tagged_multiplicity: multiplicity,
            });
            data_mass += multiplicity;
        }
    }
    let mut p0_mass = 0u64;
    for difference in 1..blocks {
        let value = gamma[difference];
        let multiplicity = (blocks - difference) as u64;
        let first = 1 ^ value;
        let second = value;
        let third = 1u64;
        let key = invariant_from_factors(
            first,
            second,
            third,
            inverse_one_plus[difference],
            gamma_inverse[difference],
            1,
        );
        records.push(TaggedRecord {
            key,
            tagged_multiplicity: P0_TAG | multiplicity,
        });
        p0_mass += multiplicity;
    }
    assert_eq!(records.len(), exact_entries);
    assert_eq!(data_mass, choose3(blocks as u64));
    assert_eq!(p0_mass, choose2(blocks as u64));

    records.sort_unstable_by_key(|record| record.key);
    let mut writer = BufWriter::with_capacity(
        16 * 1024 * 1024,
        File::create(&output_path).expect("failed to create canonical schedule"),
    );
    let mut distinct = 0u64;
    let mut maximum_data = 0u64;
    let mut maximum_p0 = 0u64;
    let mut index = 0usize;
    while index < records.len() {
        let key = records[index].key;
        let mut data = 0u64;
        let mut p0 = 0u64;
        while index < records.len() && records[index].key == key {
            let tagged = records[index].tagged_multiplicity;
            if tagged & P0_TAG == 0 {
                data += tagged;
            } else {
                p0 += tagged & !P0_TAG;
            }
            index += 1;
        }
        writer.write_all(&key.to_le_bytes()).expect("schedule write failed");
        writer.write_all(&data.to_le_bytes()).expect("schedule write failed");
        writer.write_all(&p0.to_le_bytes()).expect("schedule write failed");
        distinct += 1;
        maximum_data = maximum_data.max(data);
        maximum_p0 = maximum_p0.max(p0);
    }
    writer.flush().expect("failed to flush canonical schedule");
    let elapsed = started.elapsed().as_secs_f64();
    println!("{{");
    println!("  \"schema\": \"riffle-shiftalpha64-s3-invariant-schedule-v1\",");
    println!("  \"data_blocks\": {},", blocks);
    println!("  \"raw_gap_entries\": {},", exact_entries);
    println!("  \"distinct_invariant_keys\": {},", distinct);
    println!("  \"data_support_mass\": {},", data_mass);
    println!("  \"p0_support_mass\": {},", p0_mass);
    println!("  \"maximum_data_multiplicity\": {},", maximum_data);
    println!("  \"maximum_p0_multiplicity\": {},", maximum_p0);
    println!("  \"binary_output\": \"{}\",", output_path.replace('\\', "\\\\"));
    println!("  \"binary_record\": \"uint64 key, uint64 data multiplicity, uint64 p0 multiplicity\",");
    println!("  \"elapsed_seconds\": {:.12},", elapsed);
    println!("  \"validation\": {{");
    println!("    \"pclmul_matches_scalar_1000_probes\": true,");
    println!("    \"data_multiplicity_sums_to_choose_blocks_3\": true,");
    println!("    \"p0_multiplicity_sums_to_choose_blocks_2\": true,");
    println!("    \"invariant_records_are_sorted_and_unique\": true");
    println!("  }}");
    println!("}}");
}
