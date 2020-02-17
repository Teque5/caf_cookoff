// TODO
// Add option for FFTW/RustFFT for direct bench comparison
// Multithreading FFTs and maybe frequency shift calculations
// impl some of these on the types directly
//      xcor, freq shift, write_file on Complex64 slice
//      2d peak on Vec<Vec<Complex64>> if possible

use std::io;
use std::io::prelude::*;
use std::f64::consts::PI;
use std::fs::File;

use num_complex::Complex64;

mod xcor_fftw;
use xcor_fftw::Xcor;
// mod xcor_rustfft;
// use xcor_rustfft::Xcor;

// Reads a file of packed 32 bit floats and returns
// a Vec of its contents
pub fn read_file_c64(filename: &str) -> io::Result<Vec<Complex64>> {
    let mut f = File::open(filename)?;
    let mut buffer = Vec::new();
    f.read_to_end(&mut buffer)?;
    let mut samples = Vec::new();
    for i in (0..buffer.len()).step_by(8) {
        let mut real_bytes: [u8; 4] = Default::default();
        real_bytes.copy_from_slice(&buffer[i..i+4]);
        let real = f32::from_le_bytes(real_bytes);
        let mut imag_bytes: [u8; 4] = Default::default();
        imag_bytes.copy_from_slice(&buffer[i+4..i+8]);
        let imag = f32::from_le_bytes(imag_bytes);
        samples.push(Complex64::new(real as f64, imag as f64));
    }
    Ok(samples)
}

// Writes a slice of Complex64s to a file compatible
// with Numpy's fromfile dtype=np.complex64
#[allow(dead_code)]
pub fn write_file_complex(filename: &str, data: &[Complex64]) -> io::Result<()> {
    let mut f = File::create(filename).unwrap();
    let mut out_buf: Vec<u8> = Vec::new();

    for bin in data.iter() {
        for byte in f64::to_le_bytes(bin.re).iter() {
            out_buf.push(*byte);
        }
        for byte in f64::to_le_bytes(bin.im).iter() {
            out_buf.push(*byte);
        }
    }
    f.write_all(&out_buf).unwrap();
    Ok(())
}

// Takes in a slice of samples at samp_rate and applies
// a frequency shift to it
fn apply_freq_shifts(samples: &[Complex64], freq_shift: f64, fs: u32)
    -> Vec<Complex64> {

    // Convert (back) to vec
    let mut samples = samples.to_vec();

    // Apply to each sample
    // x *= e^(-j*2pi*fs*df*t)
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
// and compute their CAF. Return the surface as a 2D Vec
pub fn caf_surface(needle: &[Complex64], haystack: &[Complex64],
    freqs_hz: &[f64], fs: u32)
    -> Vec<Vec<Complex64>> {

    // Create our 2D surface
    let mut surface = Vec::new();

    // Run the cross correlation against the shifted ones
    let mut xcor = Xcor::new(needle.len());
    for freq in freqs_hz.iter() {
        let shifted = apply_freq_shifts(needle, *freq, fs);
        let xcor_res = xcor.run(haystack, &shifted);
        surface.push(xcor_res);
    }

    // Return our CAF surface
    surface
}

// 2D argmax
pub fn find_2d_peak(arr: Vec<Vec<Complex64>>) -> (usize, usize) {
    let mut max: Complex64 = Default::default();
    let mut argmax = (0, 0);
    for (i, row) in arr.iter().enumerate() {
        for (j, elem) in row.iter().enumerate() {
            if elem.norm_sqr() > max.norm_sqr() {
                max = *elem;
                argmax = (i, j);
            }
        }
    }
    argmax
}


// Tests for Chirp 0-8
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chip0() {
        // Read Chirp 0 reference and modified files
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_0_T+202samp_F.+69.25Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let mut shifts = Vec::new();
        for shift_millihz in (-100000..100000).step_by(250) {
            let shift = (shift_millihz as f64) / 1e3;
            shifts.push(shift);
        }

        // Get the CAF surface
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq_idx, samp_idx) = find_2d_peak(surface);

        // Confirm correct results
        assert_eq!(shifts[freq_idx], 69.25);
        assert_eq!(samp_idx, 202);
    }
}
