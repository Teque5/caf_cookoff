// TODO
// Use CLAP to take in two c64 files as arguments

mod caf;
mod utils;

use caf::{CafSurface, CafRustFFTIterRayon};
use utils::read_file_c64;

fn main() {

    // Get signals 1 and 2 to compute the caf of
    let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
    let mut haystack = read_file_c64("../data/chirp_0_T+202samp_F+69.25Hz.c64").unwrap();
    haystack.resize(needle.len(), Default::default());

    // -100Hz to 100Hz, 0.5Hz step
    let mut shifts = Vec::new();
    for shift_millihz in (-100000..100000).step_by(500) {
        let shift = (shift_millihz as f64) / 1e3;
        shifts.push(shift);
    }

    // Get the CAF surface
    let surface = CafRustFFTIterRayon::caf_surface(&needle, &haystack, &shifts, 48000);
    let (freq, samp_idx) = CafRustFFTIterRayon::find_peak(surface);

    // Print the results
    println!("Frequency offset: {:.1}Hz", freq);
    println!("Time offset: {} samples ({:.3}ms)",
        samp_idx, (samp_idx as f64) / 48.0);
}
