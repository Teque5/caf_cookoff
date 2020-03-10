use std::f64::consts::PI;
use std::sync::{mpsc, Arc};
use std::thread;

use num_complex::Complex64;

mod xcor_fftw;
mod xcor_rustfft;


// Take in 2 signals and a range of frequency shifts to try
// and compute their CAF. Return the surface as a 2D Vec of
// cross correlation magnitudes squared (for efficiency)
#[allow(dead_code)]
pub struct CafSurfaceRow {
    freq: f64,
    xcor_mag: Vec<f64>,
    xcor_peak_idx: usize,
    xcor_peak_val: f64,
}
pub trait CafSurface {

    // Every implementation will be different
    fn caf_surface(needle: &[Complex64], haystack: &[Complex64],
        freqs_hz: &[f64], fs: u32) -> Vec<CafSurfaceRow>;

    // Find the row with the highest correlation peak and return
    // its (frequency, sample_index)
    fn find_peak(arr: Vec<CafSurfaceRow>) -> (f64, usize) {
        let mut max: &CafSurfaceRow = &CafSurfaceRow {
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

    // Takes in a slice of samples at samp_rate and applies
    // a frequency shift to it
    fn apply_freq_shift(samples: &[Complex64], freq_shift: f64, fs: u32)
        -> Vec<Complex64> {

        // Convert (back) to vec
        let mut samples = samples.to_vec();

        // Apply to each sample
        // x *= e^(j*2pi*fs*df*t)
        let dt = 1.0 / (fs as f64);
        for (i, samp) in samples.iter_mut().enumerate() {
            *samp *= Complex64::from_polar(&(1.0),
                &(2.0 * PI * freq_shift * dt * (i as f64)));
        }

        // Return our shifted samples
        samples
    }
}
pub struct CafFFTW {} // FFTW one thread
impl CafSurface for CafFFTW {

    fn caf_surface(needle: &[Complex64], haystack: &[Complex64],
        freqs_hz: &[f64], fs: u32) -> Vec<CafSurfaceRow> {

        // Create our 2D surface and setup Vecs
        let mut surface = Vec::new();
        let mut needle = needle.to_vec();
        let mut haystack = haystack.to_vec();

        // Zero-pad our inputs to 2N
        needle.resize(needle.len() * 2, Default::default());
        haystack.resize(haystack.len() * 2, Default::default());

        // Run the cross correlation against the shifted ones
        let mut xcor = xcor_fftw::Xcor::new(needle.len());
        for freq in freqs_hz.iter() {

            // Generate a shifted copy and cross correlate with target
            let shifted = Self::apply_freq_shift(&needle, *freq, fs);
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

            // Push our result to the surface
            surface.push(CafSurfaceRow {
                freq: *freq,
                xcor_mag,
                xcor_peak_idx: argmax,
                xcor_peak_val: max,
            });
        }

        // Return our CAF surface
        surface
    }
}

pub struct CafRustFFT {} // RustFFT one thread
impl CafSurface for CafRustFFT {

    fn caf_surface(needle: &[Complex64], haystack: &[Complex64],
        freqs_hz: &[f64], fs: u32) -> Vec<CafSurfaceRow> {

        // Create our 2D surface and setup Vecs
        let mut surface = Vec::new();
        let mut needle = needle.to_vec();
        let mut haystack = haystack.to_vec();

        // Zero-pad our inputs to 2N
        needle.resize(needle.len() * 2, Default::default());
        haystack.resize(haystack.len() * 2, Default::default());

        // Run the cross correlation against the shifted ones
        let mut xcor = xcor_rustfft::Xcor::new(needle.len());
        for freq in freqs_hz.iter() {

            // Generate a shifted copy and cross correlate with target
            let shifted = Self::apply_freq_shift(&needle, *freq, fs);
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

            // Push our result to the surface
            surface.push(CafSurfaceRow {
                freq: *freq,
                xcor_mag,
                xcor_peak_idx: argmax,
                xcor_peak_val: max,
            });
        }

        // Return our CAF surface
        surface
    }
}

pub struct CafRustFFTThreads {} // RustFFT using std::threads
impl CafSurface for CafRustFFTThreads {

    fn caf_surface(needle: &[Complex64], haystack: &[Complex64],
        freqs_hz: &[f64], fs: u32) -> Vec<CafSurfaceRow> {

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
        let xcor = xcor_rustfft::Xcor::new(needle.len());
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
                let shifted = Self::apply_freq_shift(&needle, freq, fs);
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
                tx.send(CafSurfaceRow {
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
}
