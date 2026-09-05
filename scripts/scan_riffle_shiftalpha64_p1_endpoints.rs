//! Exact BCH-weight histograms for endpoint triples containing p1.
//!
//! A p1 support contains two finite projective points with coefficients a,b.
//! Its two finite field values are equal.  If a finite coordinate is
//! normalized to the all-one BCH message e, the p1 value is e(a+b).  If p1
//! is normalized to e, both finite values are e/(a+b).

use std::arch::x86_64::{
    _mm_clmulepi64_si128, _mm_cvtsi128_si64, _mm_set_epi64x, _mm_srli_si128,
};
use std::env;
use std::time::Instant;

const REDUCTION: u64 = 0x1b;
const ALL_ONE_MESSAGE: u64 = 0x8e63cf44efd4fa21;
const DEFAULT_DATA_BLOCKS: usize = 16_384;

#[derive(Clone, Copy, Default)]
struct Codeword {
    low: u64,
    high: u64,
}

fn load_bch_rows() -> Vec<Codeword> {
    let source = include_str!("scan_riffle_shiftalpha64_endpoint_triples.cpp");
    let start = source.find("constexpr std::array<Codeword, 64> kRows{{")
        .expect("missing BCH row table");
    let tail = &source[start..];
    let end = tail.find("}};").expect("unterminated BCH row table");
    let table = &tail[..end];
    let mut values = Vec::new();
    let mut rest = table;
    while let Some(position) = rest.find("0x") {
        let hexadecimal = &rest[position + 2..];
        let digits: String = hexadecimal
            .chars()
            .take_while(|character| character.is_ascii_hexdigit())
            .collect();
        values.push(u64::from_str_radix(&digits, 16).expect("invalid BCH row"));
        rest = &hexadecimal[digits.len()..];
    }
    assert_eq!(values.len(), 128, "wrong BCH row count");
    values
        .chunks_exact(2)
        .map(|pair| Codeword { low: pair[0], high: pair[1] })
        .collect()
}

fn initialize_tables(rows: &[Codeword]) -> Vec<Vec<Codeword>> {
    let mut tables = (0..4)
        .map(|_| vec![Codeword::default(); 1 << 16])
        .collect::<Vec<_>>();
    for chunk in 0..4 {
        for value in 1usize..(1 << 16) {
            let bit = value.trailing_zeros() as usize;
            let previous = value & (value - 1);
            tables[chunk][value] = Codeword {
                low: tables[chunk][previous].low ^ rows[16 * chunk + bit].low,
                high: tables[chunk][previous].high ^ rows[16 * chunk + bit].high,
            };
        }
    }
    tables
}

#[inline(always)]
fn bch_weight(message: u64, tables: &[Vec<Codeword>]) -> usize {
    let mut low = 0u64;
    let mut high = 0u64;
    for chunk in 0..4 {
        let row = tables[chunk][((message >> (16 * chunk)) & 0xffff) as usize];
        low ^= row.low;
        high ^= row.high;
    }
    (low.count_ones() + high.count_ones()) as usize
}

#[inline(always)]
fn multiply_x(value: u64) -> u64 {
    (value << 1) ^ (REDUCTION & 0u64.wrapping_sub(value >> 63))
}

#[target_feature(enable = "pclmulqdq")]
#[inline]
unsafe fn multiply(left: u64, right: u64) -> u64 {
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

#[inline]
fn field_multiply(left: u64, right: u64) -> u64 {
    // The executable checks the CPU feature before entering the scan.
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
    assert_ne!(value, 0, "inverse of zero");
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

fn print_histogram(name: &str, histogram: &[u64; 129], trailing_comma: bool) {
    println!("  \"{}\": {{", name);
    let entries = histogram
        .iter()
        .enumerate()
        .filter(|(_, count)| **count != 0)
        .collect::<Vec<_>>();
    for (index, (weight, count)) in entries.iter().enumerate() {
        let comma = if index + 1 == entries.len() { "" } else { "," };
        println!("    \"{}\": {}{}", weight, count, comma);
    }
    println!("  }}{}", if trailing_comma { "," } else { "" });
}

fn main() {
    if !std::is_x86_feature_detected!("pclmulqdq") {
        panic!("pclmulqdq is required");
    }
    let blocks = env::args()
        .nth(1)
        .map(|value| value.parse::<usize>().expect("invalid block count"))
        .unwrap_or(DEFAULT_DATA_BLOCKS);
    assert!((2..=DEFAULT_DATA_BLOCKS).contains(&blocks));

    let rows = load_bch_rows();
    let tables = initialize_tables(&rows);
    for probe in 0..1000u64 {
        let left = 0x0123456789abcdefu64
            .wrapping_add(probe.wrapping_mul(0x9e3779b97f4a7c15));
        let right = 0xfedcba9876543210u64
            .wrapping_add(probe.wrapping_mul(0xd1342543de82ef95));
        assert_eq!(field_multiply(left, right), multiply_scalar(left, right));
    }
    assert_eq!(field_multiply(2, inverse(2)), 1);
    assert_eq!(bch_weight(ALL_ONE_MESSAGE, &tables), 128);

    let mut gamma = vec![1u64; blocks];
    let mut gamma_inverse = vec![1u64; blocks];
    let inverse_gamma = inverse(2);
    for index in 1..blocks {
        gamma[index] = multiply_x(gamma[index - 1]);
        gamma_inverse[index] = field_multiply(gamma_inverse[index - 1], inverse_gamma);
    }
    let alpha_zero = power(2, 64);
    let inverse_alpha_zero = inverse(alpha_zero);
    let mut alpha = vec![0u64; blocks];
    let mut alpha_inverse = vec![0u64; blocks];
    alpha[0] = alpha_zero;
    alpha_inverse[0] = inverse_alpha_zero;
    for index in 1..blocks {
        alpha[index] = multiply_x(alpha[index - 1]);
        alpha_inverse[index] = field_multiply(alpha_inverse[index - 1], inverse_gamma);
    }

    let mut direct_factor = vec![0u64; blocks];
    let mut inverse_factor = vec![0u64; blocks];
    for difference in 1..blocks {
        let one_plus = 1 ^ gamma[difference];
        direct_factor[difference] = field_multiply(ALL_ONE_MESSAGE, one_plus);
        inverse_factor[difference] = field_multiply(ALL_ONE_MESSAGE, inverse(one_plus));
    }

    let started = Instant::now();
    let mut finite_histogram = [0u64; 129];
    let mut p1_histogram = [0u64; 129];
    let mut overlap_supports = 0u64;

    // The finite pair consists of p0 and one data point.
    for index in 0..blocks {
        let finite_value = field_multiply(ALL_ONE_MESSAGE, alpha[index]);
        let p1_value = field_multiply(ALL_ONE_MESSAGE, alpha_inverse[index]);
        finite_histogram[bch_weight(finite_value, &tables)] += 1;
        p1_histogram[bch_weight(p1_value, &tables)] += 1;
        overlap_supports += u64::from(finite_value == ALL_ONE_MESSAGE);
    }

    // The finite pair consists of two data points.  Factoring out the first
    // exponent leaves one precomputed multiplier for each difference.
    for difference in 1..blocks {
        let count = blocks - difference;
        let forward = direct_factor[difference];
        let backward = inverse_factor[difference];
        for index in 0..count {
            let finite_value = field_multiply(forward, alpha[index]);
            let p1_value = field_multiply(backward, alpha_inverse[index]);
            finite_histogram[bch_weight(finite_value, &tables)] += 1;
            p1_histogram[bch_weight(p1_value, &tables)] += 1;
            overlap_supports += u64::from(finite_value == ALL_ONE_MESSAGE);
        }
    }

    let supports = (blocks as u64 + 1) * blocks as u64 / 2;
    assert_eq!(finite_histogram.iter().sum::<u64>(), supports);
    assert_eq!(p1_histogram.iter().sum::<u64>(), supports);
    let elapsed = started.elapsed().as_secs_f64();

    println!("{{");
    println!("  \"schema\": \"riffle-shiftalpha64-p1-endpoint-histograms-v1\",");
    println!("  \"data_blocks\": {},", blocks);
    println!("  \"second_parity_supports\": {},", supports);
    println!("  \"all_three_all_one_overlap_supports\": {},", overlap_supports);
    print_histogram("finite_coordinate_normalization_histogram", &finite_histogram, true);
    print_histogram("p1_coordinate_normalization_histogram", &p1_histogram, true);
    println!("  \"elapsed_seconds\": {:.12},", elapsed);
    println!("  \"validation\": {{");
    println!("    \"pclmul_matches_scalar_1000_probes\": true,");
    println!("    \"all_one_bch_weight_is_128\": true,");
    println!("    \"each_histogram_sums_to_support_count\": true");
    println!("  }}");
    println!("}}");
}
