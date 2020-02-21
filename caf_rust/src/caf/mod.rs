// TODO
// Add option for FFTW/RustFFT for direct bench comparison
// Multithreading FFTs and maybe frequency shift calculations


use std::io;
use std::io::prelude::*;
use std::f64::consts::PI;
use std::fs::File;
use std::sync::{mpsc, Arc};
use std::thread;

use num_complex::Complex64;

// mod xcor_fftw;
// use xcor_fftw::Xcor;
mod xcor_rustfft;
use xcor_rustfft::Xcor;


// Reads a file of packed 32 bit floats and returns
// a Vec of its contents as Complex64 (2x 64 bit floats)
pub fn read_file_c64(filename: &str) -> io::Result<Vec<Complex64>> {

    // Open and read a file
    let mut f = File::open(filename)?;
    let mut buffer = Vec::new();
    f.read_to_end(&mut buffer)?;
    let mut samples = Vec::new();

    // Read each real and imaginary components to a buffer
    for i in (0..buffer.len()).step_by(8) {

        // Real
        let mut real_bytes: [u8; 4] = Default::default();
        real_bytes.copy_from_slice(&buffer[i..i+4]);
        let real = f32::from_le_bytes(real_bytes);

        // Image
        let mut imag_bytes: [u8; 4] = Default::default();
        imag_bytes.copy_from_slice(&buffer[i+4..i+8]);
        let imag = f32::from_le_bytes(imag_bytes);

        // Push to the calling Vec
        samples.push(Complex64::new(real as f64, imag as f64));
    }
    Ok(samples)
}

// Read/write a slice of Complex64's to/from a file
// compatible with numpy's fromfile function
pub trait BinaryIO {
    fn write_file_binary(&self, filename: &str) -> io::Result<()>;
}
// numpy dtype=np.complex128
impl BinaryIO for Vec<Complex64> {

    fn write_file_binary(&self, filename: &str) -> io::Result<()> {

        // Open the file and output byte buffer
        let mut f = File::create(filename).unwrap();
        let mut out_buf: Vec<u8> = Vec::new();

        // Write out the real and imag components one-by-one
        for samp in self.iter() {
            for byte in f64::to_le_bytes(samp.re).iter() {
                out_buf.push(*byte);
            }
            for byte in f64::to_le_bytes(samp.im).iter() {
                out_buf.push(*byte);
            }
        }
        f.write_all(&out_buf).unwrap();
        Ok(())
    }
}

// Takes in a slice of samples at samp_rate and applies
// a frequency shift to it
fn apply_freq_shifts(samples: &[Complex64], freq_shift: f64, fs: u32)
    -> Vec<Complex64> {

    // Convert (back) to vec
    let mut samples = samples.to_vec();

    // Apply to each sample
    // x *= e^(j*2pi*fs*df*t)
    let dt = 1.0 / (fs as f64);
    let exp_common = Complex64::new(0.0, 2.0 * PI * dt * freq_shift);
    for (i, samp) in samples.iter_mut().enumerate() {
        let exp = Complex64::new(i as f64, 0.0) * exp_common;
        *samp *= Complex64::exp(&exp);
    }

    // Return our shifted samples
    samples
}

// Take in 2 signals and a range of frequency shifts to try
// and compute their CAF. Return the surface as a 2D Vec of
// cross correlation magnitudes squared (for efficiency)
#[allow(dead_code)]
pub struct CAFSurfaceRow {
    freq: f64,
    xcor_mag: Vec<f64>,
    xcor_peak_idx: usize,
    xcor_peak_val: f64,
}
pub fn caf_surface(needle: &[Complex64], haystack: &[Complex64],
    freqs_hz: &[f64], fs: u32)
    -> Vec<CAFSurfaceRow> {

    // Create our 2D surface and setup Vecs
    let mut surface = Vec::new();
    let mut needle = needle.to_vec();
    let mut haystack = haystack.to_vec();

    // Zero-pad our inputs to 2N
    needle.resize(needle.len() * 2, Default::default());
    haystack.resize(haystack.len() * 2, Default::default());

    // Setup threading channel and atomic immutable references
    let (tx, rx) = mpsc::channel();
    let needle = Arc::new(needle);
    let haystack = Arc::new(haystack);

    // Run the cross correlation against the shifted ones
    let xcor = Xcor::new(needle.len());
    for freq in freqs_hz.iter() {

        // Copy what we need to for the thread
        let tx = mpsc::Sender::clone(&tx);
        let freq = *freq; // f64 can be copied, &f64 cannot
        // Copy atomic immutable reference instead of full slice
        let needle = Arc::clone(&needle);
        let haystack = Arc::clone(&haystack);
        let mut xcor = xcor.clone(); // Also calls Arc::clone in impl

        // Spawn the thread and run
        thread::spawn(move || {

            // Generate a shifted copy and cross correlate with target
            let shifted = apply_freq_shifts(&needle, freq, fs);
            let xcor_res = xcor.run(&haystack, &shifted);

            // Take the magnitude squared of the result and find (arg)max
            let mut xcor_mag = Vec::with_capacity(xcor_res.len());
            let mut max = Default::default();
            let mut argmax = 0;
            for (i, res) in xcor_res.iter().enumerate() {
                // Use the magnitude squared (for efficiency)
                let mag_squared = res.norm_sqr();
                if mag_squared > max {
                    max = mag_squared;
                    argmax = i;
                }
                xcor_mag.push(mag_squared);
            }

            // Return our result to the main thread
            tx.send(CAFSurfaceRow {
                freq,
                xcor_mag,
                xcor_peak_idx: argmax,
                xcor_peak_val: max,
            }).unwrap();
        });
    }

    // Wait for all threads to finish and
    // Populate our results into a Vec, no longer ordered by freq
    for _ in freqs_hz.iter() {
        let row = rx.recv().unwrap();
        surface.push(row);
    }

    // Return our CAF surface
    surface
}

// Find the row with the highest correlation peak and return
// its (frequency, sample_index)
pub fn find_caf_peak(arr: Vec<CAFSurfaceRow>) -> (f64, usize) {
    let mut max: &CAFSurfaceRow = &CAFSurfaceRow {
        freq: 0.0, xcor_mag: Vec::new(),
        xcor_peak_idx: 0, xcor_peak_val: 0.0
    };
    for row in arr.iter() {
        if row.xcor_peak_val > max.xcor_peak_val {
            max = &row;
        }
    }
    (max.freq, max.xcor_peak_idx)
}


// Tests for Chirp 0-8
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chip0() {
        // Read Chirp 0 reference and modified files
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_0_T+202samp_F+69.25Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 69.25);
        assert_eq!(samp_idx, 202);
    }

    #[test]
    fn test_chirp1() {
        // Read Chirp 1 reference and modified files
        let needle = read_file_c64("../data/chirp_1_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_1_T+78samp_F+35.99Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -50Hz to 50Hz, 1Hz step
        let shifts = gen_float_shifts(-50.0, 50.0, 1.0);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 36.0);
        assert_eq!(samp_idx, 78);
    }

    #[test]
    fn test_chirp2() {
        // Read Chirp 2 reference and modified files
        let needle = read_file_c64("../data/chirp_2_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_2_T+169samp_F+32.16Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // 30Hz to 35Hz, 0.05Hz step
        let shifts = gen_float_shifts(30.0, 35.0, 0.05);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 32.15);
        assert_eq!(samp_idx, 169);
    }

    #[test]
    fn test_chip3() {
        // Read Chirp 3 reference and modified files
        let needle = read_file_c64("../data/chirp_3_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_3_T+151samp_F-76.22Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, -76.25);
        assert_eq!(samp_idx, 151);
    }

    #[test]
    fn test_chip4() {
        // Read Chirp 4 reference and modified files
        let needle = read_file_c64("../data/chirp_4_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_4_T+70samp_F+82.89Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // 80Hz to 100Hz, 0.1Hz step
        let shifts = gen_float_shifts(80.0, 100.0, 0.1);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 82.9);
        assert_eq!(samp_idx, 70);
    }

    #[test]
    fn test_chip5() {
        // Read Chirp 5 reference and modified files
        let needle = read_file_c64("../data/chirp_5_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_5_T+177samp_F-92.72Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, -92.75);
        assert_eq!(samp_idx, 177);
    }

    #[test]
    fn test_chip6() {
        // Read Chirp 6 reference and modified files
        let needle = read_file_c64("../data/chirp_6_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_6_T+15samp_F-49.69Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, -49.75);
        assert_eq!(samp_idx, 15);
    }

    #[test]
    fn test_chip7() {
        // Read Chirp 7 reference and modified files
        let needle = read_file_c64("../data/chirp_7_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_7_T+84samp_F+68.26Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 68.25);
        assert_eq!(samp_idx, 84);
    }

    #[test]
    fn test_chip8() {
        // Read Chirp 8 reference and modified files
        let needle = read_file_c64("../data/chirp_8_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_8_T+80samp_F-46.28Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, -46.25);
        assert_eq!(samp_idx, 80);
    }

    #[test]
    fn test_chip9() {
        // Read Chirp 9 reference and modified files
        let (needle, haystack) = load_files(
            "../data/chirp_9_raw.c64",
            "../data/chirp_9_T+176samp_F+61.49Hz.c64");

        // -100Hz to 100Hz, 0.5Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.5);

        // Get the CAF estimates
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = find_caf_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 61.5);
        assert_eq!(samp_idx, 176);
    }

    // Helper to load and trim files by filename
    fn load_files(needle_filename: &str, haystack_filename: &str)
        -> (Vec<Complex64>, Vec<Complex64>) {

        // Read in the two files
        let needle = read_file_c64(needle_filename).unwrap();
        let mut haystack = read_file_c64(haystack_filename).unwrap();

        // Truncate haystack if necessary
        haystack.resize(needle.len(), Default::default());

        // Return our two equal-length files
        (needle, haystack)
    }

    // Helper to easily generate a range of shifts
    // Tightest resolution is 1e-3 Hz (integer truncation)
    fn gen_float_shifts(start: f64, end: f64, step: f64)
        -> Vec<f64> {

        // Convert from float to mHz to be able to
        // iterate in for loop
        let start_milli_hz = (start * 1000.0) as i32;
        let end_milli_hz = (end * 1000.0) as i32;
        let step_milli_hz = (step * 1000.0) as usize;

        // Generate shifts as floats and return the vec
        let mut shifts = Vec::new();
        for shift_millihz in
            (start_milli_hz..end_milli_hz).step_by(step_milli_hz) {
            let shift = (shift_millihz as f64) / 1e3;
            shifts.push(shift);
        }
        shifts
    }
}
