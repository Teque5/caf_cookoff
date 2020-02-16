use caf_rust::caf::{read_file_c64, caf_surface, find_2d_peak};

fn main() {

    // Get signals 1 and 2 to compute the caf of
    let needle = read_file_c64("../data/chirp_4_raw.c64").unwrap();
    let haystack = read_file_c64("../data/chirp_4_T+70samp_F.+82.89Hz.c64").unwrap();
    let haystack = &haystack[..needle.len()];

    // -100Hz to 100Hz, 0.5Hz step
    let mut shifts = Vec::new();
    for shift_millihz in (-100000..100000).step_by(500) {
        let shift = (shift_millihz as f32) / 1e3;
        shifts.push(shift);
    }

    // Get the CAF surface
    let surface = caf_surface(&needle, &haystack, &shifts, 48000);
    let (freq_idx, samp_idx) = find_2d_peak(surface);

    // Print the results
    println!("Frequency offset: {:.1}Hz", -shifts[freq_idx]);
    println!("Time offset: {} samples ({:.3}ms)",
        4096-samp_idx, ((4096-samp_idx) as f32) / 48.0);
}
