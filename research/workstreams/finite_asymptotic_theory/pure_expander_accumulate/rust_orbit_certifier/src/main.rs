use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, ToPrimitive, Zero};
use rayon::prelude::*;
use serde::Serialize;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Clone, Copy, Debug)]
struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    fn point(x: f64) -> Self {
        Self { lo: x, hi: x }
    }

    fn zero() -> Self {
        Self::point(0.0)
    }

    fn one() -> Self {
        Self::point(1.0)
    }

    fn add(self, rhs: Self) -> Self {
        Self {
            lo: next_down(self.lo + rhs.lo),
            hi: next_up(self.hi + rhs.hi),
        }
    }

    fn sub(self, rhs: Self) -> Self {
        Self {
            lo: next_down(self.lo - rhs.hi),
            hi: next_up(self.hi - rhs.lo),
        }
    }

    fn mul(self, rhs: Self) -> Self {
        let products = [
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        ];
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for value in products {
            lo = lo.min(value);
            hi = hi.max(value);
        }
        Self {
            lo: next_down(lo),
            hi: next_up(hi),
        }
    }

    fn nonnegative(self) -> Self {
        assert!(self.hi >= 0.0 && self.lo.is_finite() && self.hi.is_finite());
        Self {
            lo: self.lo.max(0.0),
            hi: self.hi.max(0.0),
        }
    }

    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    fn scale_pow2(self, exponent: i32) -> Self {
        Self {
            lo: scalbn(self.lo, exponent),
            hi: scalbn(self.hi, exponent),
        }
    }

    fn pow_nonnegative(self, mut exponent: usize) -> Self {
        let mut result = Self::one();
        let mut factor = self;
        while exponent != 0 {
            if exponent & 1 != 0 {
                result = result.mul(factor);
            }
            exponent >>= 1;
            if exponent != 0 {
                factor = factor.mul(factor);
            }
        }
        result
    }

    fn divide_by_positive(self, denominator: Self) -> Self {
        assert!(denominator.lo > 0.0);
        let lo_denominator = if self.lo >= 0.0 {
            denominator.hi
        } else {
            denominator.lo
        };
        let hi_denominator = if self.hi >= 0.0 {
            denominator.lo
        } else {
            denominator.hi
        };
        Self {
            lo: next_down(self.lo / lo_denominator),
            hi: next_up(self.hi / hi_denominator),
        }
    }

    fn finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite()
    }

    fn maximum_absolute(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }
}

fn next_up(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x == -0.0 {
        return f64::from_bits(1);
    }
    let bits = x.to_bits();
    f64::from_bits(if x >= 0.0 { bits + 1 } else { bits - 1 })
}

fn next_down(x: f64) -> f64 {
    if x.is_nan() || x == f64::NEG_INFINITY {
        return x;
    }
    if x == 0.0 {
        return -f64::from_bits(1);
    }
    let bits = x.to_bits();
    f64::from_bits(if x > 0.0 { bits - 1 } else { bits + 1 })
}

fn scalbn(x: f64, exponent: i32) -> f64 {
    x * (2.0f64).powi(exponent)
}

fn add_upper(left: f64, right: f64) -> f64 {
    next_up(left + right)
}

fn multiply_upper(left: f64, right: f64) -> f64 {
    assert!(left >= 0.0 && right >= 0.0);
    next_up(left * right)
}

fn add_lower(left: f64, right: f64) -> f64 {
    next_down(left + right).max(0.0)
}

fn delta_product(left: Interval, right: Interval) -> Interval {
    left.add(right).add(left.mul(right))
}

fn delta_power(delta: Interval, mut exponent: usize) -> Interval {
    let mut result = Interval::zero();
    let mut factor = delta;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = delta_product(result, factor);
        }
        exponent >>= 1;
        if exponent != 0 {
            factor = delta_product(factor, factor);
        }
    }
    result
}

fn binomial(n: usize, k: usize) -> BigUint {
    if k > n {
        return BigUint::zero();
    }
    let k = k.min(n - k);
    let mut value = BigUint::one();
    for index in 0..k {
        value *= BigUint::from(n - index);
        value /= BigUint::from(index + 1);
    }
    value
}

fn krawtchouk(n: usize, degree: usize, weight: usize) -> BigInt {
    let low = degree.saturating_sub(n - weight);
    let high = degree.min(weight);
    let mut total = BigInt::zero();
    for overlap in low..=high {
        let term = BigInt::from_biguint(
            Sign::Plus,
            binomial(weight, overlap) * binomial(n - weight, degree - overlap),
        );
        if overlap & 1 == 0 {
            total += term;
        } else {
            total -= term;
        }
    }
    total
}

fn f64_as_ratio(value: f64) -> (BigInt, BigUint) {
    assert!(value.is_finite());
    if value == 0.0 {
        return (BigInt::zero(), BigUint::one());
    }
    let bits = value.to_bits();
    let negative = bits >> 63 != 0;
    let exponent_bits = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1u64 << 52) - 1);
    let (mantissa, exponent) = if exponent_bits == 0 {
        (fraction, -1074)
    } else {
        ((1u64 << 52) | fraction, exponent_bits - 1023 - 52)
    };
    let mut numerator = BigInt::from(mantissa);
    let mut denominator = BigUint::one();
    if exponent >= 0 {
        numerator <<= exponent as usize;
    } else {
        denominator <<= (-exponent) as usize;
    }
    if negative {
        numerator = -numerator;
    }
    (numerator, denominator)
}

fn compare_float_ratio(value: f64, numerator: &BigInt, denominator: &BigUint) -> Ordering {
    let (float_numerator, float_denominator) = f64_as_ratio(value);
    let left = float_numerator * BigInt::from_biguint(Sign::Plus, denominator.clone());
    let right = numerator * BigInt::from_biguint(Sign::Plus, float_denominator);
    left.cmp(&right)
}

fn rational_interval(numerator: &BigInt, denominator: &BigUint) -> Interval {
    assert!(!denominator.is_zero());
    if numerator.is_zero() {
        return Interval::zero();
    }
    let approximation = numerator.to_f64().expect("numerator exceeds f64 range")
        / denominator.to_f64().expect("denominator exceeds f64 range");
    assert!(approximation.is_finite());
    match compare_float_ratio(approximation, numerator, denominator) {
        Ordering::Equal => Interval::point(approximation),
        Ordering::Less => {
            let mut lo = approximation;
            let mut hi = next_up(lo);
            while compare_float_ratio(hi, numerator, denominator) == Ordering::Less {
                lo = hi;
                hi = next_up(hi);
            }
            Interval { lo, hi }
        }
        Ordering::Greater => {
            let mut hi = approximation;
            let mut lo = next_down(hi);
            while compare_float_ratio(lo, numerator, denominator) == Ordering::Greater {
                hi = lo;
                lo = next_down(lo);
            }
            Interval { lo, hi }
        }
    }
}

fn integer_interval(value: &BigUint) -> Interval {
    rational_interval(&BigInt::from_biguint(Sign::Plus, value.clone()), &BigUint::one())
}

fn pair_type_orbit(row: [usize; 4]) -> Vec<[usize; 4]> {
    let [a, b, c, d] = row;
    let complements = [[a, b, c, d], [c, d, a, b], [b, a, d, c], [d, c, b, a]];
    let mut orbit = BTreeSet::new();
    for item in complements {
        orbit.insert(item);
        orbit.insert([item[0], item[2], item[1], item[3]]);
    }
    orbit.into_iter().collect()
}

#[derive(Clone)]
struct PairType {
    first_one: Interval,
    second_one: Interval,
    determinant: Interval,
    scale: Interval,
}

fn make_pair_types(message_bits: usize, right_degree: usize) -> (Vec<PairType>, usize) {
    let denominator = binomial(message_bits, right_degree);
    let denominator_bigint = BigInt::from_biguint(Sign::Plus, denominator.clone());
    let determinant_denominator = BigUint::from(4u8) * &denominator * &denominator;
    let bias_numerators: Vec<BigInt> = (0..=message_bits)
        .map(|weight| krawtchouk(message_bits, right_degree, weight))
        .collect();
    let bias_intervals: Vec<Interval> = bias_numerators
        .iter()
        .map(|value| rational_interval(value, &denominator))
        .collect();
    let mut result = Vec::new();
    let mut covered_types = 0usize;
    for n11 in 0..=message_bits {
        for n01 in 0..=message_bits - n11 {
            for n10 in 0..=message_bits - n11 - n01 {
                let n00 = message_bits - n11 - n01 - n10;
                let row = [n00, n01, n10, n11];
                let orbit = pair_type_orbit(row);
                if row != orbit[0] {
                    continue;
                }
                covered_types += orbit.len();
                let first_weight = n10 + n11;
                let second_weight = n01 + n11;
                let difference_weight = n01 + n10;
                let first_one = Interval::one()
                    .sub(bias_intervals[first_weight])
                    .scale_pow2(-1);
                let second_one = Interval::one()
                    .sub(bias_intervals[second_weight])
                    .scale_pow2(-1);
                let determinant_numerator =
                    &bias_numerators[difference_weight] * &denominator_bigint
                        - &bias_numerators[first_weight] * &bias_numerators[second_weight];
                let determinant =
                    rational_interval(&determinant_numerator, &determinant_denominator);
                let multiplicity = BigUint::from(orbit.len())
                    * binomial(message_bits, n11)
                    * binomial(message_bits - n11, n01)
                    * binomial(message_bits - n11 - n01, n10);
                let scale = multiplicity << message_bits.saturating_sub(2);
                result.push(PairType {
                    first_one,
                    second_one,
                    determinant,
                    scale: integer_interval(&scale),
                });
            }
        }
    }
    (result, covered_types)
}

fn coefficient_intervals(degree: usize, weight: usize) -> Vec<Interval> {
    let low = (2 * weight).saturating_sub(degree);
    (low..=weight)
        .map(|intersection| {
            let cross = weight - intersection;
            integer_interval(
                &(binomial(weight, intersection) * binomial(degree - weight, cross)),
            )
        })
        .collect()
}

fn probability_from_marginals(
    first_one: Interval,
    second_one: Interval,
    determinant: Interval,
) -> [Interval; 4] {
    let first_zero = Interval::one().sub(first_one);
    let second_zero = Interval::one().sub(second_one);
    let a0 = first_zero.mul(second_zero).nonnegative();
    let b0 = first_zero.mul(second_one).nonnegative();
    let c0 = first_one.mul(second_zero).nonnegative();
    let d0 = first_one.mul(second_one).nonnegative();
    [
        a0.add(determinant).nonnegative(),
        b0.sub(determinant).nonnegative(),
        c0.sub(determinant).nonnegative(),
        d0.add(determinant).nonnegative(),
    ]
}

fn symmetric_diagonal(
    probability: [Interval; 4],
    degree: usize,
    weight: usize,
    coefficients: &[Interval],
) -> Interval {
    let [a, b, c, d] = probability;
    let bc = b.mul(c).nonnegative();
    let low = (2 * weight).saturating_sub(degree);
    let mut result = Interval::zero();
    for (offset, coefficient) in coefficients.iter().enumerate() {
        let intersection = low + offset;
        let cross = weight - intersection;
        let zero_zero = degree - 2 * weight + intersection;
        result = result.add(
            coefficient
                .mul(a.pow_nonnegative(zero_zero))
                .mul(bc.pow_nonnegative(cross))
                .mul(d.pow_nonnegative(intersection)),
        );
    }
    result
}

fn centered_diagonal(
    first_one: Interval,
    second_one: Interval,
    determinant: Interval,
    output_bits: usize,
    level: usize,
    coefficients: &[Interval],
) -> Interval {
    let first_zero = Interval::one().sub(first_one);
    let second_zero = Interval::one().sub(second_one);
    let a0 = first_zero.mul(second_zero).nonnegative();
    let b0 = first_zero.mul(second_one).nonnegative();
    let c0 = first_one.mul(second_zero).nonnegative();
    let d0 = first_one.mul(second_one).nonnegative();
    let a = a0.add(determinant).nonnegative();
    let b = b0.sub(determinant).nonnegative();
    let c = c0.sub(determinant).nonnegative();
    let d = d0.add(determinant).nonnegative();
    let bc = b.mul(c).nonnegative();
    let bc0 = b0.mul(c0).nonnegative();
    let bc_difference = determinant.mul(determinant.sub(b0.add(c0)));
    let low = (2 * level).saturating_sub(output_bits);
    let mut result = Interval::zero();
    for (offset, coefficient) in coefficients.iter().enumerate() {
        let intersection = low + offset;
        let cross = level - intersection;
        let zero_zero = output_bits - 2 * level + intersection;
        let baseline = coefficient
            .mul(a0.pow_nonnegative(zero_zero))
            .mul(bc0.pow_nonnegative(cross))
            .mul(d0.pow_nonnegative(intersection));
        let current = coefficient
            .mul(a.pow_nonnegative(zero_zero))
            .mul(bc.pow_nonnegative(cross))
            .mul(d.pow_nonnegative(intersection));
        let mut term = current.sub(baseline);
        let regular = (zero_zero == 0 || a0.lo > 0.0)
            && (cross == 0 || bc0.lo > 0.0)
            && (intersection == 0 || d0.lo > 0.0);
        if regular {
            let mut ratios = Vec::with_capacity(3);
            if zero_zero != 0 {
                ratios.push((determinant.divide_by_positive(a0), zero_zero));
            }
            if cross != 0 {
                ratios.push((bc_difference.divide_by_positive(bc0), cross));
            }
            if intersection != 0 {
                ratios.push((determinant.divide_by_positive(d0), intersection));
            }
            let stable = ratios
                .iter()
                .all(|(ratio, _)| ratio.finite() && ratio.maximum_absolute() <= 0.5);
            if stable {
                let mut ratio_delta = Interval::zero();
                for (ratio, exponent) in ratios {
                    ratio_delta = delta_product(ratio_delta, delta_power(ratio, exponent));
                }
                term = baseline.mul(ratio_delta);
            }
        }
        result = result.add(term);
    }
    result
}

fn sector_zero_orbit(
    pair: &PairType,
    output_bits: usize,
    level: usize,
    coefficients: &[Interval],
) -> Interval {
    let first_zero = Interval::one().sub(pair.first_one);
    let second_zero = Interval::one().sub(pair.second_one);
    let negative = pair.determinant.neg();
    let transforms = [
        (pair.first_one, pair.second_one, pair.determinant),
        (first_zero, pair.second_one, negative),
        (pair.first_one, second_zero, negative),
        (first_zero, second_zero, pair.determinant),
    ];
    let mut orbit = Interval::zero();
    for (first, second, determinant) in transforms {
        orbit = orbit.add(centered_diagonal(
            first,
            second,
            determinant,
            output_bits,
            level,
            coefficients,
        ));
    }
    pair.scale.mul(orbit)
}

fn sector_orbit(
    pair: &PairType,
    output_bits: usize,
    sector: usize,
    level: usize,
    coefficients: &[Interval],
) -> Interval {
    if sector == 0 {
        return sector_zero_orbit(pair, output_bits, level, coefficients);
    }
    let degree = output_bits - 2 * sector;
    let weight = level - sector;
    let first_zero = Interval::one().sub(pair.first_one);
    let second_zero = Interval::one().sub(pair.second_one);
    let negative = pair.determinant.neg();
    let transforms = [
        (pair.first_one, pair.second_one, pair.determinant),
        (first_zero, pair.second_one, negative),
        (pair.first_one, second_zero, negative),
        (first_zero, second_zero, pair.determinant),
    ];
    let mut orbit = Interval::zero();
    for (first, second, determinant) in transforms {
        let determinant_power = match sector {
            1 => determinant,
            2 => determinant.mul(determinant).nonnegative(),
            _ => unreachable!(),
        };
        let diagonal = symmetric_diagonal(
            probability_from_marginals(first, second, determinant),
            degree,
            weight,
            coefficients,
        );
        orbit = orbit.add(determinant_power.mul(diagonal));
    }
    pair.scale.mul(orbit)
}

#[derive(Clone, Copy)]
struct CertifiedBounds {
    positive_lower: f64,
    positive_upper: f64,
    signed: Interval,
}

fn certify_entry(
    pair_types: &[PairType],
    output_bits: usize,
    sector: usize,
    level: usize,
) -> CertifiedBounds {
    let degree = output_bits - 2 * sector;
    let weight = level - sector;
    let coefficients = coefficient_intervals(degree, weight);
    let chunk_sums: Vec<CertifiedBounds> = pair_types
        .par_chunks(512)
        .map(|chunk| {
            let mut result = CertifiedBounds {
                positive_lower: 0.0,
                positive_upper: 0.0,
                signed: Interval::zero(),
            };
            for pair in chunk {
                let contribution = sector_orbit(
                    pair,
                    output_bits,
                    sector,
                    level,
                    &coefficients,
                );
                result.positive_lower = add_lower(
                    result.positive_lower,
                    contribution.lo.max(0.0),
                );
                result.positive_upper = add_upper(
                    result.positive_upper,
                    contribution.hi.max(0.0),
                );
                result.signed = result.signed.add(contribution);
            }
            result
        })
        .collect();
    chunk_sums.into_iter().fold(
        CertifiedBounds {
            positive_lower: 0.0,
            positive_upper: 0.0,
            signed: Interval::zero(),
        },
        |mut result, chunk| {
            result.positive_lower = add_lower(result.positive_lower, chunk.positive_lower);
            result.positive_upper = add_upper(result.positive_upper, chunk.positive_upper);
            result.signed = result.signed.add(chunk.signed);
            result
        },
    )
}

#[derive(Serialize)]
struct Parameters {
    message_bits: usize,
    output_bits: usize,
    right_degree: usize,
}

#[derive(Serialize)]
struct Entry {
    sector: usize,
    level: usize,
    scaled_diagonal_upper: f64,
    positive_orbit_lower: f64,
    positive_orbit_upper: f64,
    signed_orbit_lower: f64,
    signed_orbit_upper: f64,
    selected_upper_source: &'static str,
}

#[derive(Serialize)]
struct Receipt {
    schema: &'static str,
    status: &'static str,
    parameters: Parameters,
    entries: Vec<Entry>,
    enumerated_pair_types: usize,
    covered_pair_types: usize,
    canonical_complement_swap_orbits: bool,
    worker_threads: usize,
    elapsed_seconds: f64,
    scope: Vec<&'static str>,
}

fn parse_entry(text: &str) -> (usize, usize) {
    let (sector, level) = text
        .split_once(':')
        .unwrap_or_else(|| panic!("entry must have form SECTOR:LEVEL"));
    (sector.parse().unwrap(), level.parse().unwrap())
}

fn parse_args() -> (usize, usize, usize, Vec<(usize, usize)>, PathBuf, usize) {
    let mut message_bits = 256usize;
    let mut output_bits = 512usize;
    let mut right_degree = 33usize;
    let mut entries = Vec::new();
    let mut output = PathBuf::from("sector_zero_complement_orbits_rust_outward.json");
    let mut threads = 1usize;
    let mut arguments = env::args().skip(1);
    while let Some(argument) = arguments.next() {
        let value = arguments.next().unwrap_or_else(|| panic!("missing value after {argument}"));
        match argument.as_str() {
            "--message-bits" => message_bits = value.parse().unwrap(),
            "--output-bits" => output_bits = value.parse().unwrap(),
            "--right-degree" => right_degree = value.parse().unwrap(),
            "--level" => entries.push((0, value.parse().unwrap())),
            "--entry" => entries.push(parse_entry(&value)),
            "--output" => output = PathBuf::from(value),
            "--threads" => threads = value.parse().unwrap(),
            _ => panic!("unknown argument {argument}"),
        }
    }
    assert!(!entries.is_empty());
    entries.sort_unstable();
    entries.dedup();
    (message_bits, output_bits, right_degree, entries, output, threads)
}

fn main() {
    let (message_bits, output_bits, right_degree, requested_entries, output, threads) = parse_args();
    assert!(right_degree % 2 == 1);
    assert!(message_bits >= 2);
    assert!(requested_entries.iter().all(|(sector, level)| {
        *sector <= 2 && *sector <= *level && *level <= output_bits - *sector
    }));
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let started = Instant::now();
    let (pair_types, covered_types) = make_pair_types(message_bits, right_degree);
    eprintln!("prepared {} representatives covering {} types in {:.3}s", pair_types.len(), covered_types, started.elapsed().as_secs_f64());
    let entries: Vec<Entry> = requested_entries
        .iter()
        .map(|(sector, level)| {
            let bounds = certify_entry(&pair_types, output_bits, *sector, *level);
            let signed_upper = bounds.signed.hi.max(0.0);
            let (value, source) = if signed_upper <= bounds.positive_upper {
                (signed_upper, "signed_global_orbit_sum")
            } else {
                (bounds.positive_upper, "positive_orbit_sum")
            };
            eprintln!("sector,{sector},level,{level},upper,{value:.17e},positive_lower,{:.17e},positive_upper,{:.17e},signed_lower,{:.17e},signed_upper,{:.17e},elapsed_seconds,{:.3}", bounds.positive_lower, bounds.positive_upper, bounds.signed.lo, bounds.signed.hi, started.elapsed().as_secs_f64());
            Entry {
                sector: *sector,
                level: *level,
                scaled_diagonal_upper: value,
                positive_orbit_lower: bounds.positive_lower,
                positive_orbit_upper: bounds.positive_upper,
                signed_orbit_lower: bounds.signed.lo,
                signed_orbit_upper: bounds.signed.hi,
                selected_upper_source: source,
            }
        })
        .collect();
    let receipt = Receipt {
        schema: "pure-ea-sector-zero-complement-orbits-rust-outward-v1",
        status: "OUTWARD_BINARY64_UPPER_BOUND",
        parameters: Parameters { message_bits, output_bits, right_degree },
        entries,
        enumerated_pair_types: pair_types.len(),
        covered_pair_types: covered_types,
        canonical_complement_swap_orbits: true,
        worker_threads: threads,
        elapsed_seconds: started.elapsed().as_secs_f64(),
        scope: vec![
            "Complement-orbit averaging is exact for odd row degree.",
            "Exact big-integer comparisons bracket every rational bias, determinant, coefficient, and multiplicity in binary64.",
            "Every arithmetic operation expands its result by one adjacent binary64 value.",
            "Each thread performs upward summation; the main thread combines chunk sums upward in deterministic order.",
            "The bound sums positive orbit upper endpoints and does not assume orbit positivity.",
            "The receipt bounds only the listed sector-zero diagonals.",
        ],
    };
    let encoded = serde_json::to_string_pretty(&receipt).unwrap() + "\n";
    fs::write(&output, encoded).unwrap();
    println!("receipt,{},elapsed_seconds,{:.3}", output.display(), receipt.elapsed_seconds);
}
