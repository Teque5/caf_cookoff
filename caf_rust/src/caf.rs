// TODO
// Migrate from Complex32 to Complex64 to be on parity with Go implementation
// Add option for FFTW/RustFFT for direct bench comparison
// Multithreading FFTs and maybe frequency shift calculations
// impl some of these on the types directly
//      xcor, freq shift, write_file on Complex32 slice
//      2d peak on Vec<Vec<Complex32>> if possible

use std::io;
use std::io::prelude::*;
use std::f32::consts::PI;
use std::fs::File;
use std::sync::Arc;

use fftw::array::AlignedVec;
use fftw::plan::{C2CPlan, C2CPlan32};
use fftw::types::{Sign, Flag};
use itertools::izip;
use num_complex::{Complex32, Complex};
use rustfft::{FFTplanner, num_traits::Zero, FFT};

// Reads a file of packed 32 bit floats and returns
// a Vec of its contents
pub fn read_file_c64(filename: &str) -> io::Result<Vec<Complex32>> {
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
        samples.push(Complex32::new(real, imag));
    }
    Ok(samples)
}

// Writes a slice of Complex32s to a file compatible
// with Numpy's fromfile dtype=np.complex64
#[allow(dead_code)]
pub fn write_file_c64(filename: &str, data: &[Complex32]) -> io::Result<()> {
    let mut f = File::create(filename).unwrap();
    let mut out_buf: Vec<u8> = Vec::new();

    for bin in data.iter() {
        for byte in f32::to_le_bytes(bin.re).iter() {
            out_buf.push(*byte);
        }
        for byte in f32::to_le_bytes(bin.im).iter() {
            out_buf.push(*byte);
        }
    }
    f.write_all(&out_buf).unwrap();
    Ok(())
}

// Cross correlation of 2 complex slices using FFTW
// Assumes inputs powers of 2
// Naive: ifft(fft(a) * fft(b).conj())
#[allow(dead_code)]
struct XcorFFTW {
    n: usize, // size of a, b, c
    // Aligned FFTW buffers
    a: AlignedVec<Complex32>,
    b: AlignedVec<Complex32>,
    c: AlignedVec<Complex32>,
    // Planners
    forward_planner: C2CPlan32,
    reverse_planner: C2CPlan32,
}

impl XcorFFTW {

    // Constructor
    #[allow(dead_code)]
    pub fn new(n: usize) -> Self {

        // Create planners
        let fp = C2CPlan::aligned(
            &[n], Sign::Forward, Flag::Measure).unwrap();
        let rp = C2CPlan::aligned(
            &[n], Sign::Backward, Flag::Measure).unwrap();

        // Return new struct
        XcorFFTW {
            n: n,
            a: AlignedVec::new(n),
            b: AlignedVec::new(n),
            c: AlignedVec::new(n),
            forward_planner: fp,
            reverse_planner: rp,
        }
    }

    // Run cross-correlation against any complex input buffers
    // sized N
    #[allow(dead_code)]
    pub fn run(&mut self, a: &[Complex32], b: &[Complex32]) -> Vec<Complex32> {

        // Sanity
        assert!(a.len() == self.n);
        assert!(b.len() == self.n);

        // Compute FFT(a), FFT(b)
        self.a.copy_from_slice(a);
        self.forward_planner.c2c(&mut self.a, &mut self.b).unwrap();
        self.a.copy_from_slice(b);
        self.forward_planner.c2c(&mut self.a, &mut self.c).unwrap();

        // Take complex conjugate of FFT(b) == self.c
        for bin in self.c.iter_mut() {
            *bin = bin.conj();
        }

        // Calculate FFT(a) * conj(FFT(b)) and normalize
        for (out, a, b) in izip!(self.a.iter_mut(),
                                 self.b.iter(), self.c.iter()) {

            *out = (a * b) / (self.n as f32);
        }

        // Calculate IFFT of product and return
        self.reverse_planner.c2c(&mut self.a, &mut self.b).unwrap();
        self.b.to_vec()
    }
}



// Cross correlation of 2 complex slices using RustFFT
// Assumes inputs powers of 2
// Naive: ifft(fft(a) * fft(b).conj())
#[allow(dead_code)]
struct XcorRustFFT {
    n: usize, // size of a, b, c
    // Preallocated buffers
    a: Vec<Complex32>,
    b: Vec<Complex32>,
    c: Vec<Complex32>,
    // Planners
    fft: Arc<dyn FFT<f32>>,
    ifft: Arc<dyn FFT<f32>>,
}

impl XcorRustFFT {

    // Constructor
    #[allow(dead_code)]
    pub fn new(n: usize) -> Self {

        // Create FFT plans
        let mut fp = FFTplanner::new(false);
        let fft = fp.plan_fft(n);
        let mut rp = FFTplanner::new(true);
        let ifft = rp.plan_fft(n);

        // Return new struct
        XcorRustFFT {
            n: n,
            a: vec![Default::default(); n],
            b: vec![Default::default(); n],
            c: vec![Default::default(); n],
            fft: fft,
            ifft: ifft,
        }
    }

    // Run cross-correlation against any complex input buffers
    // sized N
    #[allow(dead_code)]
    pub fn run(&mut self, a: &[Complex32], b: &[Complex32]) -> Vec<Complex32> {

        // Sanity
        assert!(a.len() == self.n);
        assert!(b.len() == self.n);

        // Compute FFT(a), FFT(b)
        self.a.copy_from_slice(a);
        self.fft.process(&mut self.a, &mut self.b);
        self.a.copy_from_slice(b);
        self.fft.process(&mut self.a, &mut self.c);

        // Take complex conjugate of FFT(b) == self.c
        for bin in self.c.iter_mut() {
            *bin = bin.conj();
        }

        // Calculate FFT(a) * conj(FFT(b)) and normalize
        for (out, a, b) in izip!(self.a.iter_mut(),
                                 self.b.iter(), self.c.iter()) {

            *out = (a * b) / (self.n as f32);
        }

        // Calculate IFFT of product and return
        self.ifft.process(&mut self.a, &mut self.b);
        self.b.to_vec()
    }
}

// Takes in a slice of samples at samp_rate and applies
// a frequency shift to it
fn apply_freq_shifts(samples: &[Complex32], freq_shift: f32, fs: u32)
    -> Vec<Complex32> {

    // Convert (back) to vec
    let mut samples = samples.to_vec();

    // Apply to each sample
    // x *= e^(-j*2pi*fs*df*t)
    let dt = 1.0 / (fs as f32);
    let exp_common = Complex32::new(0.0, -2.0 * PI * dt * freq_shift);
    for (i, samp) in samples.iter_mut().enumerate() {
        let exp = Complex32::new(i as f32, 0.0) * exp_common;
        *samp *= Complex32::exp(&exp);
    }

    // Return our shifted samples
    samples
}

// Take in 2 signals and a range of frequency shifts to try
// and compute their CAF. Return the surface as a 2D Vec
pub fn caf_surface(needle: &[Complex32], haystack: &[Complex32],
    freqs_hz: &[f32], fs: u32)
    -> Vec<Vec<Complex32>> {

    // Create our 2D surface
    let mut surface = Vec::new();

    // Run the cross correlation against the shifted ones
    // let mut xcor_fftw = XcorFFTW::new(needle.len());
    let mut xcor_rustfft = XcorRustFFT::new(needle.len());
    for freq in freqs_hz.iter() {
        let shifted = apply_freq_shifts(needle, *freq, fs);
        // let xcor_res = xcor_fftw.run(&shifted, haystack);
        let xcor_res = xcor_rustfft.run(&shifted, haystack);
        surface.push(xcor_res);
    }

    // Return our CAF surface
    surface
}

// 2D argmax
pub fn find_2d_peak(arr: Vec<Vec<Complex32>>) -> (usize, usize) {
    let mut max: Complex32 = Default::default();
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
            let shift = (shift_millihz as f32) / 1e3;
            shifts.push(shift);
        }

        // Get the CAF surface
        let surface = caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq_idx, samp_idx) = find_2d_peak(surface);

        // Confirm correct results
        assert_eq!(-shifts[freq_idx], 69.25);
        assert_eq!(4096-samp_idx, 202);
    }
}
