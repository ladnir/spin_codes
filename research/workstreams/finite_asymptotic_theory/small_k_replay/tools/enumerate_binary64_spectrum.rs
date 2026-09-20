use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::Arc;
use std::thread;

const DIMENSION: usize = 32;
const LENGTH: usize = 64;
const LOW_DIMENSION: usize = 20;
const HIGH_COUNT: u64 = 1u64 << (DIMENSION - LOW_DIMENSION);
const LOW_COUNT: u64 = 1u64 << LOW_DIMENSION;

fn read_rows(path: &str) -> [u64; DIMENSION] {
    let file = File::open(path).expect("could not open generator-row file");
    let mut rows = [0u64; DIMENSION];
    let mut lines = BufReader::new(file).lines();
    for row in &mut rows {
        let line = lines.next().expect("generator-row file is truncated").unwrap();
        *row = u64::from_str_radix(line.trim(), 16).expect("invalid hexadecimal row");
    }
    assert!(lines.next().is_none(), "generator-row file has extra rows");
    rows
}

fn main() {
    let args: Vec<String> = env::args().collect();
    assert_eq!(args.len(), 3, "usage: enumerator ROWS OUT.csv");
    let rows = Arc::new(read_rows(&args[1]));
    let threads = thread::available_parallelism().map_or(1, usize::from);
    let mut handles = Vec::with_capacity(threads);
    for worker in 0..threads {
        let rows = Arc::clone(&rows);
        handles.push(thread::spawn(move || {
            let mut histogram = [0u64; LENGTH + 1];
            let mut high = worker as u64;
            while high < HIGH_COUNT {
                let mut high_word = 0u64;
                for bit in 0..(DIMENSION - LOW_DIMENSION) {
                    if ((high >> bit) & 1) != 0 {
                        high_word ^= rows[LOW_DIMENSION + bit];
                    }
                }
                let mut word = high_word;
                histogram[word.count_ones() as usize] += 1;
                for index in 1..LOW_COUNT {
                    word ^= rows[index.trailing_zeros() as usize];
                    histogram[word.count_ones() as usize] += 1;
                }
                high += threads as u64;
            }
            histogram
        }));
    }
    let mut histogram = [0u64; LENGTH + 1];
    for handle in handles {
        let local = handle.join().expect("enumeration worker failed");
        for (total, count) in histogram.iter_mut().zip(local) {
            *total += count;
        }
    }
    let mass: u64 = histogram.iter().sum();
    assert_eq!(mass, 1u64 << DIMENSION, "spectrum mass mismatch");
    let file = File::create(&args[2]).expect("could not create spectrum file");
    let mut out = BufWriter::new(file);
    writeln!(out, "weight,count").unwrap();
    for (weight, count) in histogram.iter().enumerate() {
        if *count != 0 {
            writeln!(out, "{weight},{count}").unwrap();
        }
    }
    out.flush().unwrap();
    let minimum = histogram.iter().enumerate().skip(1).find(|(_, count)| **count != 0).unwrap().0;
    println!("threads,{threads}");
    println!("mass,{mass}");
    println!("minimum_weight,{minimum}");
    println!("output,{}", args[2]);
}
