//! Match exact EBCH (22,22,24) ratios against every finite outer support.

use std::convert::TryInto;
use std::env;
use std::fs;
use std::time::Instant;

const REDUCTION: u64 = 0x1b;
const DEFAULT_BLOCKS: usize = 16_384;
const DEFAULT_EXPECTED_RATIO_MASS: usize = 27_765_248;

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
            assert_ne!(key, 0);
            let mut index = mix(key) as usize & result.mask;
            while result.entries[index].key != 0 {
                assert_ne!(result.entries[index].key, key);
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

fn read_raw_ratios(path: &str, expected_ratio_mass: usize) -> Vec<u64> {
    let bytes = fs::read(path).expect("failed to read raw ratio stream");
    assert_eq!(bytes.len() % 8, 0, "partial raw ratio record");
    let mut ratios = Vec::with_capacity(bytes.len() / 8);
    for record in bytes.chunks_exact(8) {
        ratios.push(u64::from_le_bytes(record.try_into().unwrap()));
    }
    assert_eq!(ratios.len(), expected_ratio_mass);
    assert!(ratios.iter().all(|&ratio| ratio != 0 && ratio != 1));
    ratios
}

fn compress_ratios(mut ratios: Vec<u64>) -> Vec<(u64, u64)> {
    ratios.sort_unstable();
    let mut records = Vec::new();
    let mut begin = 0usize;
    while begin < ratios.len() {
        let ratio = ratios[begin];
        let mut end = begin + 1;
        while end < ratios.len() && ratios[end] == ratio {
            end += 1;
        }
        records.push((ratio, (end - begin) as u64));
        begin = end;
    }
    records
}

fn write_histogram(path: &str, records: &[(u64, u64)]) {
    let mut bytes = Vec::with_capacity(16 * records.len());
    for &(ratio, count) in records {
        bytes.extend_from_slice(&ratio.to_le_bytes());
        bytes.extend_from_slice(&count.to_le_bytes());
    }
    fs::write(path, bytes).expect("failed to write ratio histogram");
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

fn accumulate_hit(
    table: &RatioTable,
    ratio: u64,
    multiplicity: u64,
    ratio_hits: &mut u64,
    profile_supports: &mut u64,
    outer_words: &mut u64,
    maximum_count: &mut u64,
) {
    let count = table.get(ratio);
    if count != 0 {
        *ratio_hits += 1;
        *profile_supports += multiplicity;
        *outer_words += multiplicity * count;
        *maximum_count = (*maximum_count).max(count);
    }
}

fn main() {
    if !std::is_x86_feature_detected!("pclmulqdq") {
        panic!("pclmulqdq is required");
    }
    let path = env::args()
        .nth(1)
        .expect("usage: scanner raw-ratios.bin [blocks] [histogram-output] [expected-mass] [profile]");
    let blocks = env::args()
        .nth(2)
        .map(|value| value.parse::<usize>().expect("invalid block count"))
        .unwrap_or(DEFAULT_BLOCKS);
    let histogram_output = env::args().nth(3).filter(|value| value != "-");
    let expected_ratio_mass = env::args()
        .nth(4)
        .map(|value| value.parse::<usize>().expect("invalid expected ratio mass"))
        .unwrap_or(DEFAULT_EXPECTED_RATIO_MASS);
    let profile = env::args().nth(5).unwrap_or_else(|| "22_22_24".to_owned());
    let (schema, unique_weight, exact_key) = match profile.as_str() {
        "22_22_24" => (
            "riffle-shiftalpha64-weight222224-outer-v1",
            24,
            "exact_finite_weight222224_outer_words",
        ),
        "22_22_26" => (
            "riffle-shiftalpha64-weight222226-outer-v1",
            26,
            "exact_finite_weight222226_outer_words",
        ),
        "24_24_22" => (
            "riffle-shiftalpha64-weight242422-outer-v1",
            22,
            "exact_finite_weight242422_outer_words",
        ),
        _ => panic!("unsupported equal-pair profile"),
    };
    assert!((3..=DEFAULT_BLOCKS).contains(&blocks));

    for probe in 0..1000u64 {
        let left = 0x0123456789abcdefu64
            .wrapping_add(probe.wrapping_mul(0x9e3779b97f4a7c15));
        let right = 0xfedcba9876543210u64
            .wrapping_add(probe.wrapping_mul(0xd1342543de82ef95));
        assert_eq!(field_multiply(left, right), multiply_scalar(left, right));
    }

    let preparation_started = Instant::now();
    let records = compress_ratios(read_raw_ratios(&path, expected_ratio_mass));
    let ratio_mass: u64 = records.iter().map(|record| record.1).sum();
    assert_eq!(ratio_mass, expected_ratio_mass as u64);
    if let Some(output) = histogram_output.as_deref() {
        write_histogram(output, &records);
    }
    let distinct_ratios = records.len();
    let maximum_ratio_count = records.iter().map(|record| record.1).max().unwrap();
    let table = RatioTable::new(&records);
    let preparation_seconds = preparation_started.elapsed().as_secs_f64();

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

    let scan_started = Instant::now();
    let mut data_ratio_hits = [0u64; 3];
    let mut data_profile_supports = [0u64; 3];
    let mut data_outer_words = [0u64; 3];
    let mut p0_ratio_hits = [0u64; 3];
    let mut p0_profile_supports = [0u64; 3];
    let mut p0_outer_words = [0u64; 3];
    let mut maximum_used_ratio_count = 0u64;

    for third_difference in 2..blocks {
        let third = gamma[third_difference];
        let multiplicity = (blocks - third_difference) as u64;
        for second_difference in 1..third_difference {
            let second = gamma[second_difference];
            let common = field_multiply(
                gamma_inverse[second_difference],
                inverse_one_plus[third_difference - second_difference],
            );
            let ratios = [
                field_multiply(1 ^ second, inverse_one_plus[third_difference]),
                field_multiply(1 ^ second, common),
                field_multiply(1 ^ third, common),
            ];
            for role in 0..3 {
                accumulate_hit(
                    &table,
                    ratios[role],
                    multiplicity,
                    &mut data_ratio_hits[role],
                    &mut data_profile_supports[role],
                    &mut data_outer_words[role],
                    &mut maximum_used_ratio_count,
                );
            }
        }
    }

    for difference in 1..blocks {
        let value = gamma[difference];
        let multiplicity = (blocks - difference) as u64;
        let ratios = [
            gamma_inverse[difference],
            inverse_one_plus[difference],
            field_multiply(value, inverse_one_plus[difference]),
        ];
        for role in 0..3 {
            accumulate_hit(
                &table,
                ratios[role],
                multiplicity,
                &mut p0_ratio_hits[role],
                &mut p0_profile_supports[role],
                &mut p0_outer_words[role],
                &mut maximum_used_ratio_count,
            );
        }
    }

    let data_outer_total: u64 = data_outer_words.iter().sum();
    let p0_outer_total: u64 = p0_outer_words.iter().sum();
    let scan_seconds = scan_started.elapsed().as_secs_f64();
    println!("{{");
    println!("  \"schema\": \"{}\",", schema);
    println!("  \"profile\": \"{}\",", profile);
    println!("  \"data_blocks\": {},", blocks);
    println!("  \"raw_ratio_mass\": {},", ratio_mass);
    println!("  \"distinct_ratios\": {},", distinct_ratios);
    println!("  \"maximum_global_ratio_multiplicity\": {},", maximum_ratio_count);
    println!("  \"data_only_supports\": {},", choose3(blocks as u64));
    println!("  \"data_ratio_hits_by_weight{}_role\": {:?},", unique_weight, data_ratio_hits);
    println!("  \"data_profile_supports_by_weight{}_role\": {:?},", unique_weight, data_profile_supports);
    println!("  \"data_outer_words_by_weight{}_role\": {:?},", unique_weight, data_outer_words);
    println!("  \"p0_and_two_data_supports\": {},", blocks as u64 * (blocks as u64 - 1) / 2);
    println!("  \"p0_ratio_hits_by_weight{}_role\": {:?},", unique_weight, p0_ratio_hits);
    println!("  \"p0_profile_supports_by_weight{}_role\": {:?},", unique_weight, p0_profile_supports);
    println!("  \"p0_outer_words_by_weight{}_role\": {:?},", unique_weight, p0_outer_words);
    println!("  \"{}\": {},", exact_key, data_outer_total + p0_outer_total);
    println!("  \"maximum_used_ratio_count\": {},", maximum_used_ratio_count);
    println!("  \"preparation_seconds\": {:.12},", preparation_seconds);
    println!("  \"scan_seconds\": {:.12},", scan_seconds);
    println!("  \"validation\": {{");
    println!("    \"pclmul_matches_scalar_1000_probes\": true,");
    println!("    \"all_three_weight{}_roles_scanned\": true,", unique_weight);
    println!("    \"raw_ratio_mass_matches_exact_relation_count\": true");
    println!("  }}");
    println!("}}");
}
