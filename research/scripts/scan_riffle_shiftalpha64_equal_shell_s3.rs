//! Count exact equal-shell three-block outer words through the AGL/S3 quotient.

use std::collections::HashSet;
use std::convert::TryInto;
use std::env;
use std::fs::{self, File};
use std::io::{BufReader, Read};
use std::time::Instant;

const GROUP_SIZE: usize = 127 * 128;
const REDUCTION: u64 = 0x1b;
const DECODER_INVERSE: u64 = 0x92a451cd307d0e63;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Codeword {
    low: u64,
    high: u64,
}

const GENERATOR_ROWS: [Codeword; 64] = [
    Codeword{low:0xf4845518b9582a1f,high:0x8000000000000000},Codeword{low:0xe908aa3172b0543e,high:0x8000000000000001},
    Codeword{low:0xd2115462e560a87c,high:0x8000000000000003},Codeword{low:0xa422a8c5cac150f8,high:0x8000000000000007},
    Codeword{low:0x4845518b9582a1f0,high:0x800000000000000f},Codeword{low:0x908aa3172b0543e0,high:0x800000000000001e},
    Codeword{low:0x2115462e560a87c0,high:0x800000000000003d},Codeword{low:0x422a8c5cac150f80,high:0x800000000000007a},
    Codeword{low:0x845518b9582a1f00,high:0x80000000000000f4},Codeword{low:0x08aa3172b0543e00,high:0x80000000000001e9},
    Codeword{low:0x115462e560a87c00,high:0x80000000000003d2},Codeword{low:0x22a8c5cac150f800,high:0x80000000000007a4},
    Codeword{low:0x45518b9582a1f000,high:0x8000000000000f48},Codeword{low:0x8aa3172b0543e000,high:0x8000000000001e90},
    Codeword{low:0x15462e560a87c000,high:0x8000000000003d21},Codeword{low:0x2a8c5cac150f8000,high:0x8000000000007a42},
    Codeword{low:0x5518b9582a1f0000,high:0x800000000000f484},Codeword{low:0xaa3172b0543e0000,high:0x800000000001e908},
    Codeword{low:0x5462e560a87c0000,high:0x800000000003d211},Codeword{low:0xa8c5cac150f80000,high:0x800000000007a422},
    Codeword{low:0x518b9582a1f00000,high:0x80000000000f4845},Codeword{low:0xa3172b0543e00000,high:0x80000000001e908a},
    Codeword{low:0x462e560a87c00000,high:0x80000000003d2115},Codeword{low:0x8c5cac150f800000,high:0x80000000007a422a},
    Codeword{low:0x18b9582a1f000000,high:0x8000000000f48455},Codeword{low:0x3172b0543e000000,high:0x8000000001e908aa},
    Codeword{low:0x62e560a87c000000,high:0x8000000003d21154},Codeword{low:0xc5cac150f8000000,high:0x8000000007a422a8},
    Codeword{low:0x8b9582a1f0000000,high:0x800000000f484551},Codeword{low:0x172b0543e0000000,high:0x800000001e908aa3},
    Codeword{low:0x2e560a87c0000000,high:0x800000003d211546},Codeword{low:0x5cac150f80000000,high:0x800000007a422a8c},
    Codeword{low:0xb9582a1f00000000,high:0x80000000f4845518},Codeword{low:0x72b0543e00000000,high:0x80000001e908aa31},
    Codeword{low:0xe560a87c00000000,high:0x80000003d2115462},Codeword{low:0xcac150f800000000,high:0x80000007a422a8c5},
    Codeword{low:0x9582a1f000000000,high:0x8000000f4845518b},Codeword{low:0x2b0543e000000000,high:0x8000001e908aa317},
    Codeword{low:0x560a87c000000000,high:0x8000003d2115462e},Codeword{low:0xac150f8000000000,high:0x8000007a422a8c5c},
    Codeword{low:0x582a1f0000000000,high:0x800000f4845518b9},Codeword{low:0xb0543e0000000000,high:0x800001e908aa3172},
    Codeword{low:0x60a87c0000000000,high:0x800003d2115462e5},Codeword{low:0xc150f80000000000,high:0x800007a422a8c5ca},
    Codeword{low:0x82a1f00000000000,high:0x80000f4845518b95},Codeword{low:0x0543e00000000000,high:0x80001e908aa3172b},
    Codeword{low:0x0a87c00000000000,high:0x80003d2115462e56},Codeword{low:0x150f800000000000,high:0x80007a422a8c5cac},
    Codeword{low:0x2a1f000000000000,high:0x8000f4845518b958},Codeword{low:0x543e000000000000,high:0x8001e908aa3172b0},
    Codeword{low:0xa87c000000000000,high:0x8003d2115462e560},Codeword{low:0x50f8000000000000,high:0x8007a422a8c5cac1},
    Codeword{low:0xa1f0000000000000,high:0x800f4845518b9582},Codeword{low:0x43e0000000000000,high:0x801e908aa3172b05},
    Codeword{low:0x87c0000000000000,high:0x803d2115462e560a},Codeword{low:0x0f80000000000000,high:0x807a422a8c5cac15},
    Codeword{low:0x1f00000000000000,high:0x80f4845518b9582a},Codeword{low:0x3e00000000000000,high:0x81e908aa3172b054},
    Codeword{low:0x7c00000000000000,high:0x83d2115462e560a8},Codeword{low:0xf800000000000000,high:0x87a422a8c5cac150},
    Codeword{low:0xf000000000000000,high:0x8f4845518b9582a1},Codeword{low:0xe000000000000000,high:0x9e908aa3172b0543},
    Codeword{low:0xc000000000000000,high:0xbd2115462e560a87},Codeword{low:0x8000000000000000,high:0xfa422a8c5cac150f},
];

#[derive(Clone)]
struct Orbit {
    representative_index: usize,
    size: usize,
    stabilizers: Vec<u16>,
}

#[derive(Clone, Copy)]
struct Action {
    multiplier: u8,
    translation: u8,
}

#[inline(always)]
fn gf7_multiply(mut left: u8, mut right: u8) -> u8 {
    let mut result = 0u8;
    for _ in 0..7 {
        if right & 1 != 0 { result ^= left; }
        right >>= 1;
        let top = left & 0x40;
        left = (left << 1) & 0x7f;
        if top != 0 { left ^= 0x03; }
    }
    result
}

fn gf7_inverse(value: u8) -> u8 {
    assert_ne!(value, 0);
    (1..128).find(|&candidate| gf7_multiply(value, candidate) == 1).unwrap()
}

#[inline(always)]
fn action_from_index(index: usize) -> Action {
    Action { multiplier: (index / 128 + 1) as u8, translation: (index % 128) as u8 }
}

#[inline(always)]
fn action_index(action: Action) -> usize {
    (action.multiplier as usize - 1) * 128 + action.translation as usize
}

#[inline(always)]
fn compose(left: Action, right: Action) -> Action {
    Action {
        multiplier: gf7_multiply(left.multiplier, right.multiplier),
        translation: gf7_multiply(left.multiplier, right.translation) ^ left.translation,
    }
}

#[inline(always)]
fn inverse_action(action: Action) -> Action {
    let inverse = gf7_inverse(action.multiplier);
    Action { multiplier: inverse, translation: gf7_multiply(inverse, action.translation) }
}

fn encode(mut message: u64) -> Codeword {
    let mut result = Codeword { low: 0, high: 0 };
    while message != 0 {
        let bit = message.trailing_zeros() as usize;
        result.low ^= GENERATOR_ROWS[bit].low;
        result.high ^= GENERATOR_ROWS[bit].high;
        message &= message - 1;
    }
    result
}

fn read_messages(path: &str, expected: usize, weight: u32) -> (Vec<u64>, Vec<Codeword>) {
    let bytes = fs::read(path).expect("failed to read shell messages");
    assert_eq!(bytes.len(), expected * 8, "shell file has the wrong length");
    let mut messages = Vec::with_capacity(expected);
    for record in bytes.chunks_exact(8) {
        messages.push(u64::from_le_bytes(record.try_into().unwrap()));
    }
    assert!(messages.windows(2).all(|pair| pair[0] < pair[1]));
    let codewords: Vec<_> = messages.iter().map(|&message| encode(message)).collect();
    assert!(codewords.iter().all(|word| word.low.count_ones() + word.high.count_ones() == weight));
    (messages, codewords)
}

fn coordinate_maps() -> ([u8; 128], [u8; 128]) {
    let mut coordinate_to_point = [0u8; 128];
    let mut point = 1u8;
    for coordinate in 0..127 {
        coordinate_to_point[coordinate] = point;
        point = gf7_multiply(point, 2);
    }
    coordinate_to_point[127] = 0;
    let mut point_to_coordinate = [0u8; 128];
    for coordinate in 0..128 {
        point_to_coordinate[coordinate_to_point[coordinate] as usize] = coordinate as u8;
    }
    (coordinate_to_point, point_to_coordinate)
}

fn carryless_multiply_low(mut left: u64, mut right: u64) -> u64 {
    let mut result = 0u64;
    while right != 0 {
        if right & 1 != 0 { result ^= left; }
        right >>= 1;
        left <<= 1;
    }
    result
}

fn build_byte_action_table() -> Vec<u64> {
    let (coordinate_to_point, point_to_coordinate) = coordinate_maps();
    let mut decode_basis = [0u64; 128];
    for coordinate in 0..64 {
        decode_basis[coordinate] = carryless_multiply_low(1u64 << coordinate, DECODER_INVERSE);
    }
    let mut basis = vec![0u64; 128 * GROUP_SIZE];
    for coordinate in 0..128 {
        let point = coordinate_to_point[coordinate];
        for action_index_value in 0..GROUP_SIZE {
            let action = action_from_index(action_index_value);
            let destination_point = gf7_multiply(action.multiplier, point) ^ action.translation;
            let destination = point_to_coordinate[destination_point as usize] as usize;
            basis[coordinate * GROUP_SIZE + action_index_value] = decode_basis[destination];
        }
    }
    let mut table = vec![0u64; 16 * 256 * GROUP_SIZE];
    for byte_position in 0..16 {
        for value in 1usize..256 {
            let bit = value.trailing_zeros() as usize;
            let previous = value & (value - 1);
            let source = (byte_position * 8 + bit) * GROUP_SIZE;
            let previous_offset = (byte_position * 256 + previous) * GROUP_SIZE;
            let output_offset = (byte_position * 256 + value) * GROUP_SIZE;
            for action in 0..GROUP_SIZE {
                table[output_offset + action] = table[previous_offset + action] ^ basis[source + action];
            }
        }
    }
    table
}

#[inline(always)]
fn codeword_bytes(codeword: Codeword) -> [u8; 16] {
    let mut result = [0u8; 16];
    result[..8].copy_from_slice(&codeword.low.to_le_bytes());
    result[8..].copy_from_slice(&codeword.high.to_le_bytes());
    result
}

fn transform_all(codeword: Codeword, table: &[u64], output: &mut [u64]) {
    output.fill(0);
    for (byte_position, &value) in codeword_bytes(codeword).iter().enumerate() {
        if value == 0 { continue; }
        let offset = (byte_position * 256 + value as usize) * GROUP_SIZE;
        let row = &table[offset..offset + GROUP_SIZE];
        for action in 0..GROUP_SIZE {
            output[action] ^= row[action];
        }
    }
}

#[inline(always)]
fn transform_one(codeword: Codeword, action: usize, table: &[u64]) -> u64 {
    let mut result = 0u64;
    for (byte_position, &value) in codeword_bytes(codeword).iter().enumerate() {
        result ^= table[(byte_position * 256 + value as usize) * GROUP_SIZE + action];
    }
    result
}

fn classify_orbits(
    messages: &[u64],
    codewords: &[Codeword],
    table: &[u64],
    weight: u32,
) -> (Vec<Orbit>, Vec<u16>, Vec<u16>) {
    let mut covered = vec![false; messages.len()];
    let mut orbit_ids = vec![u16::MAX; messages.len()];
    let mut sections = vec![u16::MAX; messages.len()];
    let mut orbits = Vec::new();
    let mut transformed = vec![0u64; GROUP_SIZE];
    let mut cursor = 0usize;
    while cursor < messages.len() {
        while cursor < messages.len() && covered[cursor] { cursor += 1; }
        if cursor == messages.len() { break; }
        transform_all(codewords[cursor], table, &mut transformed);
        let mut images: Vec<(u64, u16)> = transformed.iter().enumerate()
            .map(|(action, &message)| (message, action as u16)).collect();
        images.sort_unstable();
        let orbit_id = orbits.len() as u16;
        let mut stabilizers = Vec::new();
        let mut image_index = 0usize;
        let mut orbit_size = 0usize;
        while image_index < images.len() {
            let message = images[image_index].0;
            let first_action = images[image_index].1;
            let mut end = image_index + 1;
            while end < images.len() && images[end].0 == message { end += 1; }
            let shell_index = messages.binary_search(&message).expect("affine image left shell");
            assert!(!covered[shell_index], "affine orbit overlap");
            covered[shell_index] = true;
            orbit_ids[shell_index] = orbit_id;
            sections[shell_index] = first_action;
            if message == messages[cursor] {
                stabilizers.extend(images[image_index..end].iter().map(|row| row.1));
            }
            orbit_size += 1;
            image_index = end;
        }
        assert_eq!(orbit_size * stabilizers.len(), GROUP_SIZE);
        orbits.push(Orbit { representative_index: cursor, size: orbit_size, stabilizers });
    }
    assert!(covered.iter().all(|&value| value));
    let size_4064 = orbits.iter().filter(|orbit| orbit.size == 4064).count();
    let size_16256 = orbits.iter().filter(|orbit| orbit.size == GROUP_SIZE).count();
    match weight {
        22 => assert_eq!((orbits.len(), size_4064, size_16256), (15, 0, 15)),
        24 => assert_eq!((orbits.len(), size_4064, size_16256), (532, 147, 385)),
        _ => unreachable!(),
    }
    (orbits, orbit_ids, sections)
}

fn sorted_triple(mut values: [u64; 3]) -> [u64; 3] {
    values.sort_unstable();
    values
}

fn canonical_seed(
    indices: [usize; 3],
    messages: &[u64],
    codewords: &[Codeword],
    orbits: &[Orbit],
    orbit_ids: &[u16],
    sections: &[u16],
    table: &[u64],
) -> [u64; 3] {
    let minimum_representative = indices.iter().map(|&index| {
        messages[orbits[orbit_ids[index] as usize].representative_index]
    }).min().unwrap();
    let mut best = [u64::MAX; 3];
    for &vertex in &indices {
        let orbit = &orbits[orbit_ids[vertex] as usize];
        if messages[orbit.representative_index] != minimum_representative { continue; }
        let inverse_section = inverse_action(action_from_index(sections[vertex] as usize));
        for &stabilizer in &orbit.stabilizers {
            let action = compose(action_from_index(stabilizer as usize), inverse_section);
            let action = action_index(action);
            let candidate = sorted_triple([
                transform_one(codewords[indices[0]], action, table),
                transform_one(codewords[indices[1]], action, table),
                transform_one(codewords[indices[2]], action, table),
            ]);
            best = best.min(candidate);
        }
    }
    assert_eq!(best[0], minimum_representative);
    best
}

fn build_seed_orbits(
    messages: &[u64],
    codewords: &[Codeword],
    orbits: &[Orbit],
    orbit_ids: &[u16],
    sections: &[u16],
    table: &[u64],
    weight: u32,
    expected_ordered_relations: u64,
) -> (Vec<[u64; 3]>, u64) {
    let mut seeds = HashSet::<[u64; 3]>::with_capacity(50_000);
    let mut ordered_relations = 0u64;
    for (orbit_number, orbit) in orbits.iter().enumerate() {
        let x_index = orbit.representative_index;
        let x_message = messages[x_index];
        let x_codeword = codewords[x_index];
        let mut partner_count = 0u64;
        for y_index in 0..messages.len() {
            let distance = (x_codeword.low ^ codewords[y_index].low).count_ones()
                + (x_codeword.high ^ codewords[y_index].high).count_ones();
            if distance != weight { continue; }
            partner_count += 1;
            let y_message = messages[y_index];
            let z_message = x_message ^ y_message;
            if y_message >= z_message { continue; }
            let z_index = messages.binary_search(&z_message).expect("XOR partner left shell");
            let x_rep = messages[orbit.representative_index];
            let y_rep = messages[orbits[orbit_ids[y_index] as usize].representative_index];
            let z_rep = messages[orbits[orbit_ids[z_index] as usize].representative_index];
            if x_rep > y_rep || x_rep > z_rep { continue; }
            seeds.insert(canonical_seed(
                [x_index, y_index, z_index], messages, codewords, orbits,
                orbit_ids, sections, table,
            ));
        }
        assert_eq!(partner_count & 1, 0);
        ordered_relations += orbit.size as u64 * partner_count;
        if (orbit_number + 1) % 16 == 0 || orbit_number + 1 == orbits.len() {
            eprintln!("seed scan: {}/{} orbits, {} canonical triple orbits",
                orbit_number + 1, orbits.len(), seeds.len());
        }
    }
    assert_eq!(ordered_relations, expected_ordered_relations);
    let mut result: Vec<_> = seeds.into_iter().collect();
    result.sort_unstable();
    (result, ordered_relations)
}

fn triple_stabilizer(
    seed: [u64; 3],
    messages: &[u64],
    codewords: &[Codeword],
    orbits: &[Orbit],
    orbit_ids: &[u16],
    sections: &[u16],
    table: &[u64],
) -> Vec<u16> {
    let indices = seed.map(|message| messages.binary_search(&message).unwrap());
    let mut stabilizers = HashSet::<u16>::new();
    for source in 0..3 {
        let source_orbit_id = orbit_ids[indices[source]];
        let inverse_source = inverse_action(action_from_index(sections[indices[source]] as usize));
        for target in 0..3 {
            if orbit_ids[indices[target]] != source_orbit_id { continue; }
            let target_section = action_from_index(sections[indices[target]] as usize);
            for &word_stabilizer in &orbits[source_orbit_id as usize].stabilizers {
                let action = compose(
                    target_section,
                    compose(action_from_index(word_stabilizer as usize), inverse_source),
                );
                let action_value = action_index(action);
                let image = sorted_triple([
                    transform_one(codewords[indices[0]], action_value, table),
                    transform_one(codewords[indices[1]], action_value, table),
                    transform_one(codewords[indices[2]], action_value, table),
                ]);
                if image == seed { stabilizers.insert(action_value as u16); }
            }
        }
    }
    assert!(stabilizers.contains(&0));
    assert_eq!(GROUP_SIZE % stabilizers.len(), 0);
    let mut result: Vec<_> = stabilizers.into_iter().collect();
    result.sort_unstable();
    result
}

#[target_feature(enable = "pclmulqdq")]
#[inline]
unsafe fn multiply(left: u64, right: u64) -> u64 {
    use std::arch::x86_64::{_mm_clmulepi64_si128, _mm_cvtsi128_si64, _mm_set_epi64x, _mm_srli_si128};
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
fn field_multiply(left: u64, right: u64) -> u64 { unsafe { multiply(left, right) } }

fn power(mut base: u64, mut exponent: u64) -> u64 {
    let mut result = 1u64;
    while exponent != 0 {
        if exponent & 1 != 0 { result = field_multiply(result, base); }
        exponent >>= 1;
        if exponent != 0 { base = field_multiply(base, base); }
    }
    result
}

fn field_inverse(value: u64) -> u64 {
    assert_ne!(value, 0);
    power(value, u64::MAX - 1)
}

fn batch_invert(values: &mut [u64], prefix: &mut Vec<u64>) {
    prefix.clear();
    prefix.reserve(values.len().saturating_sub(prefix.capacity()));
    let mut product = 1u64;
    for &value in values.iter() {
        assert_ne!(value, 0);
        prefix.push(product);
        product = field_multiply(product, value);
    }
    let mut inverse = field_inverse(product);
    for index in (0..values.len()).rev() {
        let value = values[index];
        values[index] = field_multiply(inverse, prefix[index]);
        inverse = field_multiply(inverse, value);
    }
}

fn invariant_components(first: u64, second: u64) -> (u64, u64) {
    let third = first ^ second;
    let first_squared = field_multiply(first, first);
    let second_squared = field_multiply(second, second);
    let product_first_second = field_multiply(first, second);
    let numerator_base = first_squared ^ product_first_second ^ second_squared;
    let numerator = field_multiply(field_multiply(numerator_base, numerator_base), numerator_base);
    let product = field_multiply(product_first_second, third);
    let denominator = field_multiply(product, product);
    (numerator, denominator)
}

fn expand_relation_keys(
    seeds: &[[u64; 3]],
    messages: &[u64],
    codewords: &[Codeword],
    orbits: &[Orbit],
    orbit_ids: &[u16],
    sections: &[u16],
    table: &[u64],
    expected_unordered: usize,
) -> (Vec<u64>, usize, usize) {
    let mut keys = Vec::<u64>::with_capacity(expected_unordered);
    let mut first_images = vec![0u64; GROUP_SIZE];
    let mut second_images = vec![0u64; GROUP_SIZE];
    let mut numerators = Vec::<u64>::with_capacity(GROUP_SIZE);
    let mut denominators = Vec::<u64>::with_capacity(GROUP_SIZE);
    let mut inverse_prefix = Vec::<u64>::with_capacity(GROUP_SIZE);
    let mut nontrivial_stabilizers = 0usize;
    let mut special_keys = 0usize;
    for (seed_number, &seed) in seeds.iter().enumerate() {
        let stabilizer = triple_stabilizer(
            seed, messages, codewords, orbits, orbit_ids, sections, table,
        );
        if stabilizer.len() != 1 { nontrivial_stabilizers += 1; }
        transform_all(encode(seed[0]), table, &mut first_images);
        transform_all(encode(seed[1]), table, &mut second_images);
        numerators.clear();
        denominators.clear();
        for action in 0..GROUP_SIZE {
            if stabilizer.len() != 1 && !stabilizer.iter().all(|&right| {
                action <= action_index(compose(action_from_index(action), action_from_index(right as usize)))
            }) { continue; }
            let (numerator, denominator) = invariant_components(first_images[action], second_images[action]);
            numerators.push(numerator);
            denominators.push(denominator);
        }
        assert_eq!(numerators.len(), GROUP_SIZE / stabilizer.len());
        batch_invert(&mut denominators, &mut inverse_prefix);
        for index in 0..numerators.len() {
            let key = field_multiply(numerators[index], denominators[index]);
            if key == 0 { special_keys += 1; }
            keys.push(key);
        }
        if (seed_number + 1) % 64 == 0 || seed_number + 1 == seeds.len() {
            eprintln!("relation expansion: {}/{} triple orbits, {} unordered triples",
                seed_number + 1, seeds.len(), keys.len());
        }
    }
    assert_eq!(keys.len(), expected_unordered);
    (keys, nontrivial_stabilizers, special_keys)
}

fn bucketed_sort(keys: &mut [u64]) {
    const BUCKETS: usize = 1 << 16;
    let mut counts = vec![0usize; BUCKETS];
    for &key in keys.iter() { counts[(key >> 48) as usize] += 1; }
    let mut starts = vec![0usize; BUCKETS + 1];
    for bucket in 0..BUCKETS { starts[bucket + 1] = starts[bucket] + counts[bucket]; }
    let mut next = starts[..BUCKETS].to_vec();
    for bucket in 0..BUCKETS {
        while next[bucket] < starts[bucket + 1] {
            let position = next[bucket];
            let target = (keys[position] >> 48) as usize;
            if target == bucket {
                next[bucket] += 1;
            } else {
                let target_position = next[target];
                keys.swap(position, target_position);
                next[target] += 1;
            }
        }
    }
    for bucket in 0..BUCKETS {
        keys[starts[bucket]..starts[bucket + 1]].sort_unstable();
    }
    assert!(keys.windows(2).all(|pair| pair[0] <= pair[1]));
}

fn merge_schedule(keys: &[u64], schedule_path: &str) -> (u128, u128, u64, u64, u64, u64) {
    let file = File::open(schedule_path).expect("failed to open invariant schedule");
    let bytes = file.metadata().unwrap().len() as usize;
    assert_eq!(bytes % 24, 0);
    let records = bytes / 24;
    let mut reader = BufReader::with_capacity(16 * 1024 * 1024, file);
    let mut buffer = vec![0u8; 24 * 65_536];
    let mut relation_index = 0usize;
    let mut data_outer = 0u128;
    let mut p0_outer = 0u128;
    let mut matched_keys = 0u64;
    let mut matched_relations = 0u64;
    let mut distinct_relations = 0u64;
    let mut maximum_relation_multiplicity = 0u64;
    let mut previous_schedule = None;
    let mut records_read = 0usize;
    while records_read < records {
        let take = (records - records_read).min(65_536);
        reader.read_exact(&mut buffer[..take * 24]).expect("partial schedule record");
        for record in buffer[..take * 24].chunks_exact(24) {
            let key = u64::from_le_bytes(record[0..8].try_into().unwrap());
            let data = u64::from_le_bytes(record[8..16].try_into().unwrap());
            let p0 = u64::from_le_bytes(record[16..24].try_into().unwrap());
            if let Some(previous) = previous_schedule { assert!(previous < key); }
            previous_schedule = Some(key);
            while relation_index < keys.len() && keys[relation_index] < key {
                let value = keys[relation_index];
                let begin = relation_index;
                while relation_index < keys.len() && keys[relation_index] == value { relation_index += 1; }
                let count = (relation_index - begin) as u64;
                distinct_relations += 1;
                maximum_relation_multiplicity = maximum_relation_multiplicity.max(count);
            }
            if relation_index < keys.len() && keys[relation_index] == key {
                let begin = relation_index;
                while relation_index < keys.len() && keys[relation_index] == key { relation_index += 1; }
                let count = (relation_index - begin) as u64;
                let automorphisms = if key == 0 { 3u128 } else { 1u128 };
                data_outer += count as u128 * data as u128 * automorphisms;
                p0_outer += count as u128 * p0 as u128 * automorphisms;
                matched_keys += 1;
                matched_relations += count;
                distinct_relations += 1;
                maximum_relation_multiplicity = maximum_relation_multiplicity.max(count);
            }
        }
        records_read += take;
    }
    while relation_index < keys.len() {
        let value = keys[relation_index];
        let begin = relation_index;
        while relation_index < keys.len() && keys[relation_index] == value { relation_index += 1; }
        let count = (relation_index - begin) as u64;
        distinct_relations += 1;
        maximum_relation_multiplicity = maximum_relation_multiplicity.max(count);
    }
    (data_outer, p0_outer, matched_keys, matched_relations, distinct_relations, maximum_relation_multiplicity)
}

fn multiply_scalar(mut left: u64, mut right: u64) -> u64 {
    let mut result = 0u64;
    while right != 0 {
        if right & 1 != 0 { result ^= left; }
        right >>= 1;
        left = (left << 1) ^ (REDUCTION & 0u64.wrapping_sub(left >> 63));
    }
    result
}

fn main() {
    if !std::is_x86_feature_detected!("pclmulqdq") { panic!("pclmulqdq is required"); }
    let arguments: Vec<String> = env::args().collect();
    if arguments.len() != 4 { panic!("usage: scanner shell-messages.bin invariant-schedule.bin weight"); }
    let weight: u32 = arguments[3].parse().expect("invalid shell weight");
    let (expected_shell, expected_ordered, expected_unordered) = match weight {
        22 => (243_840usize, 1_365_504u64, 227_584usize),
        24 => (6_855_968usize, 4_018_336_896u64, 669_722_816usize),
        _ => panic!("only weight 22 and 24 shells are supported"),
    };
    for probe in 0..1000u64 {
        let left = 0x0123456789abcdefu64.wrapping_add(probe.wrapping_mul(0x9e3779b97f4a7c15));
        let right = 0xfedcba9876543210u64.wrapping_add(probe.wrapping_mul(0xd1342543de82ef95));
        assert_eq!(field_multiply(left, right), multiply_scalar(left, right));
    }

    let started = Instant::now();
    let (messages, codewords) = read_messages(&arguments[1], expected_shell, weight);
    let loaded_seconds = started.elapsed().as_secs_f64();
    let table = build_byte_action_table();
    let action_table_seconds = started.elapsed().as_secs_f64() - loaded_seconds;
    assert_eq!(transform_one(codewords[0], 0, &table), messages[0]);
    let (orbits, orbit_ids, sections) = classify_orbits(&messages, &codewords, &table, weight);
    let classified_seconds = started.elapsed().as_secs_f64() - loaded_seconds - action_table_seconds;
    let seed_started = Instant::now();
    let (seeds, ordered_relations) = build_seed_orbits(
        &messages, &codewords, &orbits, &orbit_ids, &sections, &table,
        weight, expected_ordered,
    );
    let seed_seconds = seed_started.elapsed().as_secs_f64();
    let expansion_started = Instant::now();
    let (mut keys, nontrivial_triple_stabilizers, special_relation_keys) = expand_relation_keys(
        &seeds, &messages, &codewords, &orbits, &orbit_ids, &sections, &table,
        expected_unordered,
    );
    let expansion_seconds = expansion_started.elapsed().as_secs_f64();
    drop(table);
    let sort_started = Instant::now();
    bucketed_sort(&mut keys);
    let sort_seconds = sort_started.elapsed().as_secs_f64();
    let merge_started = Instant::now();
    let (data_outer, p0_outer, matched_keys, matched_relations, distinct_relations, maximum_relation_multiplicity) =
        merge_schedule(&keys, &arguments[2]);
    let merge_seconds = merge_started.elapsed().as_secs_f64();
    let total_seconds = started.elapsed().as_secs_f64();
    println!("{{");
    println!("  \"schema\": \"riffle-shiftalpha64-equal-shell-s3-scan-v1\",");
    println!("  \"shell_weight\": {},", weight);
    println!("  \"shell_size\": {},", messages.len());
    println!("  \"affine_word_orbits\": {},", orbits.len());
    println!("  \"affine_triple_orbits\": {},", seeds.len());
    println!("  \"triple_orbits_with_nontrivial_stabilizer\": {},", nontrivial_triple_stabilizers);
    println!("  \"ordered_additive_relations\": {},", ordered_relations);
    println!("  \"unordered_additive_relations\": {},", keys.len());
    println!("  \"distinct_relation_invariant_keys\": {},", distinct_relations);
    println!("  \"maximum_relation_invariant_multiplicity\": {},", maximum_relation_multiplicity);
    println!("  \"special_J_zero_relations\": {},", special_relation_keys);
    println!("  \"matched_invariant_keys\": {},", matched_keys);
    println!("  \"matched_unordered_relations\": {},", matched_relations);
    println!("  \"exact_data_only_outer_words\": {},", data_outer);
    println!("  \"exact_p0_outer_words\": {},", p0_outer);
    println!("  \"exact_finite_outer_words\": {},", data_outer + p0_outer);
    println!("  \"timing_seconds\": {{");
    println!("    \"load_and_encode\": {:.12},", loaded_seconds);
    println!("    \"action_table\": {:.12},", action_table_seconds);
    println!("    \"word_orbit_classification\": {:.12},", classified_seconds);
    println!("    \"triple_orbit_seeds\": {:.12},", seed_seconds);
    println!("    \"relation_expansion\": {:.12},", expansion_seconds);
    println!("    \"bucketed_sort\": {:.12},", sort_seconds);
    println!("    \"schedule_merge\": {:.12},", merge_seconds);
    println!("    \"total\": {:.12}", total_seconds);
    println!("  }},");
    println!("  \"validation\": {{");
    println!("    \"pclmul_matches_scalar_1000_probes\": true,");
    println!("    \"committed_shell_sorted_unique_and_reencoded\": true,");
    println!("    \"affine_word_orbit_histogram_matches_receipt\": true,");
    println!("    \"ordered_relation_count_matches_goal14\": true,");
    println!("    \"triple_orbit_sizes_sum_to_unordered_relation_count\": true,");
    println!("    \"relation_invariant_keys_sorted\": true,");
    println!("    \"schedule_records_sorted_unique\": true");
    println!("  }}");
    println!("}}");
}
