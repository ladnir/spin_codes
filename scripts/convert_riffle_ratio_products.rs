//! Convert affine ratio-product sections into a raw GF(2^64) ratio stream.

use std::env;
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::time::Instant;

const EXPECTED_SECTIONS: usize = 532;
const EXPECTED_MASS: u64 = 321_121_024;
const EXPECTED_INVERSES: u64 = 6_855_968;

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

fn read_u32(reader: &mut BufReader<File>) -> u32 {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes).expect("partial section header");
    u32::from_le_bytes(bytes)
}

fn read_u64_vector(reader: &mut BufReader<File>, count: usize) -> Vec<u64> {
    let mut bytes = vec![0u8; 8 * count];
    reader.read_exact(&mut bytes).expect("partial uint64 vector");
    bytes
        .chunks_exact(8)
        .map(|chunk| u64::from_le_bytes(chunk.try_into().unwrap()))
        .collect()
}

fn main() {
    if !std::is_x86_feature_detected!("pclmulqdq") {
        panic!("pclmulqdq is required");
    }
    let input_path = env::args().nth(1).expect(
        "usage: converter ratio-products.bin raw-ratios.bin [expected-mass]",
    );
    let output_path = env::args().nth(2).expect(
        "usage: converter ratio-products.bin raw-ratios.bin [expected-mass]",
    );
    let expected_mass = env::args()
        .nth(3)
        .map(|value| value.parse::<u64>().expect("invalid expected mass"))
        .unwrap_or(EXPECTED_MASS);
    let started = Instant::now();
    let mut reader = BufReader::with_capacity(
        16 * 1024 * 1024,
        File::open(input_path).expect("failed to open ratio products"),
    );
    let mut writer = BufWriter::with_capacity(
        16 * 1024 * 1024,
        File::create(output_path).expect("failed to create raw ratio stream"),
    );
    let mut mass = 0u64;
    let mut inverse_values = 0u64;
    let mut output = Vec::<u8>::new();
    for section in 0..EXPECTED_SECTIONS {
        let orbit_size = read_u32(&mut reader) as usize;
        let partner_count = read_u32(&mut reader) as usize;
        assert!(orbit_size == 4_064 || orbit_size == 16_256);
        let inverses = read_u64_vector(&mut reader, orbit_size);
        assert!(inverses.iter().all(|&value| value != 0));
        inverse_values += orbit_size as u64;
        output.resize(8 * orbit_size, 0);
        for _ in 0..partner_count {
            let numerators = read_u64_vector(&mut reader, orbit_size);
            for index in 0..orbit_size {
                let ratio = field_multiply(numerators[index], inverses[index]);
                assert!(ratio != 0 && ratio != 1);
                output[8 * index..8 * index + 8]
                    .copy_from_slice(&ratio.to_le_bytes());
            }
            writer.write_all(&output).expect("failed to write ratio vector");
            mass += orbit_size as u64;
        }
        if (section + 1) % 16 == 0 || section + 1 == EXPECTED_SECTIONS {
            eprintln!(
                "native ratios: {}/{} sections, {}/{} records",
                section + 1,
                EXPECTED_SECTIONS,
                mass,
                expected_mass
            );
        }
    }
    let mut trailing = [0u8; 1];
    assert_eq!(reader.read(&mut trailing).expect("failed final read"), 0);
    writer.flush().expect("failed to flush raw ratio stream");
    assert_eq!(mass, expected_mass);
    assert_eq!(inverse_values, EXPECTED_INVERSES);
    println!("{{");
    println!("  \"schema\": \"riffle-ratio-product-conversion-v1\",");
    println!("  \"sections\": {},", EXPECTED_SECTIONS);
    println!("  \"inverse_values\": {},", inverse_values);
    println!("  \"raw_ratio_mass\": {},", mass);
    println!("  \"elapsed_seconds\": {:.12},", started.elapsed().as_secs_f64());
    println!("  \"validation\": {{");
    println!("    \"pclmul_available\": true,");
    println!("    \"input_consumed_exactly\": true,");
    println!("    \"all_ratios_nonzero_and_not_one\": true");
    println!("  }}");
    println!("}}");
}
