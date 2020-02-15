use std::io;
use std::io::prelude::*;
use std::fs::File;

use num_complex::Complex32;

// Reads a file of packed 32 bit floats and returns
// a Vec of its contents
fn read_file_f32(filename: &str) -> io::Result<Vec<f32>> {
    let mut f = File::open(filename)?;
    let mut buffer = Vec::new();
    f.read_to_end(&mut buffer)?;
    let mut floats = Vec::new();
    for i in (0..buffer.len()).step_by(4) {
        let mut bytes: [u8; 4] = Default::default();
        bytes.copy_from_slice(&buffer[i..i+4]);
        let float = f32::from_le_bytes(bytes);
        floats.push(float);
    }
    Ok(floats)
}

// Converts a slice of 32 bit floats to the same
// data as num_complex::Complex32
fn to_complex_32(in_slice: &[f32]) -> Vec<Complex32> {

    let mut samples = Vec::new();
    for real in in_slice.iter() {
        samples.push(Complex32::new(*real, 0.0)); // re + 0j
    }
    samples

}

fn main() {
    let data_real = read_file_f32("../data/burst_0000_raw.f32").unwrap();
    let data = to_complex_32(&data_real);
    for sample in data.iter() {
        println!("{}", sample);
    }
}
