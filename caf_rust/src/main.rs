use std::io;
use std::io::prelude::*;
use std::fs::File;


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

fn main() {
    let data_real = read_file_f32("../data/burst_0000_raw.f32").unwrap();
    for sample in data_real.iter() {
        println!("{}", sample);
    }
}
