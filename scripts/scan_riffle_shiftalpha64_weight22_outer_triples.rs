//! Match exact EBCH weight-22 additive ratios against finite outer supports.

use std::env;
use std::convert::TryInto;
use std::fs;
use std::time::Instant;

const REDUCTION: u64 = 0x1b;
const DEFAULT_BLOCKS: usize = 16_384;

#[derive(Clone, Copy, Default)]
struct Entry {
    key: u64,
    count: u64,
}

struct RatioTable {
    entries: Vec<Entry>,
    mask: usize,
}

#[inline(always)]
fn mix(mut value: u64) -> u64 {
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

impl RatioTable {
    fn new(records: &[(u64, u64)]) -> Self {
        let size = (2 * records.len()).next_power_of_two();
        let mut result = Self {
            entries: vec![Entry::default(); size],
            mask: size - 1,
        };
        for &(key, count) in records {
            assert_ne!(key, 0, "zero ratio cannot occur in a minimum triple");
            let mut index = mix(key) as usize & result.mask;
            while result.entries[index].key != 0 {
                assert_ne!(result.entries[index].key, key, "duplicate ratio record");
                index = (index + 1) & result.mask;
            }
            result.entries[index] = Entry { key, count };
        }
        result
    }

    #[inline(always)]
    fn get(&self, key: u64) -> u64 {
        let mut index = mix(key) as usize & self.mask;
        loop {
            let entry = self.entries[index];
            if entry.key == key {
                return entry.count;
            }
            if entry.key == 0 {
                return 0;
            }
            index = (index + 1) & self.mask;
        }
    }
}

fn read_records(path: &str) -> Vec<(u64, u64)> {
    let bytes = fs::read(path).expect("failed to read ratio histogram");
    assert_eq!(bytes.len() % 16, 0, "partial ratio record");
    let mut records = Vec::with_capacity(bytes.len() / 16);
    for record in bytes.chunks_exact(16) {
        let ratio = u64::from_le_bytes(record[..8].try_into().unwrap());
        let count = u64::from_le_bytes(record[8..].try_into().unwrap());
        records.push((ratio, count));
    }
    assert!(records.windows(2).all(|pair| pair[0].0 < pair[1].0));
    records
}

#[inline(always)]
fn multiply_x(value: u64) -> u64 {
    (value << 1) ^ (REDUCTION & 0u64.wrapping_sub(value >> 63))
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

fn choose3(value: u64) -> u64 {
    value * (value - 1) * (value - 2) / 6
}

fn main() {
    if !std::is_x86_feature_detected!("pclmulqdq") {
        panic!("pclmulqdq is required");
    }
    let path = env::args().nth(1).expect("usage: scanner ratios.bin [blocks]");
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

    let records = read_records(&path);
    let ratio_mass: u64 = records.iter().map(|record| record.1).sum();
    assert_eq!(ratio_mass, 1_365_504);
    let table = RatioTable::new(&records);

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

    let started = Instant::now();
    let mut data_difference_pairs_with_hits = 0u64;
    let mut data_supports_with_hits = 0u64;
    let mut data_outer_words = 0u64;
    let mut p0_differences_with_hits = 0u64;
    let mut p0_supports_with_hits = 0u64;
    let mut p0_outer_words = 0u64;
    let mut maximum_used_ratio_count = 0u64;

    for third_difference in 2..blocks {
        let numerator = 1 ^ gamma[third_difference];
        let multiplicity = (blocks - third_difference) as u64;
        for second_difference in 1..third_difference {
            let ratio = field_multiply(
                field_multiply(numerator, gamma_inverse[second_difference]),
                inverse_one_plus[third_difference - second_difference],
            );
            let count = table.get(ratio);
            if count != 0 {
                data_difference_pairs_with_hits += 1;
                data_supports_with_hits += multiplicity;
                data_outer_words += multiplicity * count;
                maximum_used_ratio_count = maximum_used_ratio_count.max(count);
            }
        }
    }

    for difference in 1..blocks {
        let ratio = field_multiply(gamma[difference], inverse_one_plus[difference]);
        let count = table.get(ratio);
        if count != 0 {
            let multiplicity = (blocks - difference) as u64;
            p0_differences_with_hits += 1;
            p0_supports_with_hits += multiplicity;
            p0_outer_words += multiplicity * count;
            maximum_used_ratio_count = maximum_used_ratio_count.max(count);
        }
    }

    let data_supports = choose3(blocks as u64);
    let p0_supports = blocks as u64 * (blocks as u64 - 1) / 2;
    let elapsed = started.elapsed().as_secs_f64();
    println!("{{");
    println!("  \"schema\": \"riffle-shiftalpha64-weight22-outer-triples-v1\",");
    println!("  \"data_blocks\": {},", blocks);
    println!("  \"ratio_records\": {},", records.len());
    println!("  \"ratio_ordered_pair_mass\": {},", ratio_mass);
    println!("  \"data_only\": {{");
    println!("    \"supports\": {},", data_supports);
    println!("    \"difference_pairs_with_hits\": {},", data_difference_pairs_with_hits);
    println!("    \"supports_with_hits\": {},", data_supports_with_hits);
    println!("    \"exact_weight22_triple_outer_words\": {}", data_outer_words);
    println!("  }},");
    println!("  \"p0_and_two_data\": {{");
    println!("    \"supports\": {},", p0_supports);
    println!("    \"differences_with_hits\": {},", p0_differences_with_hits);
    println!("    \"supports_with_hits\": {},", p0_supports_with_hits);
    println!("    \"exact_weight22_triple_outer_words\": {}", p0_outer_words);
    println!("  }},");
    println!("  \"exact_finite_weight22_triple_outer_words\": {},", data_outer_words + p0_outer_words);
    println!("  \"maximum_compatible_scalars_on_one_used_ratio\": {},", maximum_used_ratio_count);
    println!("  \"elapsed_seconds\": {:.12},", elapsed);
    println!("  \"validation\": {{");
    println!("    \"pclmul_matches_scalar_1000_probes\": true,");
    println!("    \"ratio_histogram_mass_is_1365504\": true");
    println!("  }}");
    println!("}}");
}
