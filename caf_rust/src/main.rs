use std::io;
use std::io::prelude::*;
use std::f32::consts::PI;
use std::fs::File;

use fftw::array::AlignedVec;
use fftw::plan::{C2CPlan, C2CPlan32};
use fftw::types::{Sign, Flag};
use itertools::izip;
use num_complex::{Complex32, Complex};
use rustfft::{FFTplanner, num_traits::Zero};

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

// Writes a slice of Complex32s to a file compatible
// with Numpy's fromfile dtype=np.complex64
fn write_file_c32(filename: &str, data: &[Complex32]) -> io::Result<()> {
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

// Converts a slice of 32 bit floats to the same
// data as num_complex::Complex32
fn to_complex_32(in_slice: &[f32]) -> Vec<Complex32> {

    let mut samples = Vec::new();
    for real in in_slice.iter() {
        samples.push(Complex32::new(*real, 0.0)); // re + 0j
    }
    samples
}

// Cross correlation of 2 complex slices using FFTW
// Assumes inputs powers of 2
// Naive: ifft(fft(a) * fft(b).conj())
#[allow(dead_code)]
fn xcor_fftw(a: &[Complex32], b: &[Complex32]) -> Vec<Complex32> {

    // Sanity
    assert!(a.len() == b.len());

    // Allocations
    let n = a.len();
    let mut a_time = AlignedVec::new(n);
    let mut b_time = AlignedVec::new(n);
    let mut a_freq = AlignedVec::new(n);
    let mut b_freq = AlignedVec::new(n);
    let mut res_freq = AlignedVec::new(n);
    let mut res_time = AlignedVec::new(n);
    a_time.copy_from_slice(a);
    b_time.copy_from_slice(b);

    // Compute FFT(a), FFT(b)
    let mut forward_planner: C2CPlan32 = C2CPlan::aligned(
        &[n], Sign::Forward, Flag::Measure).unwrap();
    forward_planner.c2c(&mut a_time, &mut a_freq).unwrap();
    forward_planner.c2c(&mut b_time, &mut b_freq).unwrap();

    // Take complex conjugate of b
    for bin in b_freq.iter_mut() {
        *bin = bin.conj();
    }

    // Calculate a*b and normalize
    for (out, a, b) in izip!(res_freq.iter_mut(),
                             a_freq.iter(), b_freq.iter()) {
        *out = (a * b) / (n as f32);
    }

    // IFFT
    let mut reverse_planner: C2CPlan32 = C2CPlan::aligned(
        &[n], Sign::Backward, Flag::Measure).unwrap();
    reverse_planner.c2c(&mut res_freq, &mut res_time).unwrap();

    // Return IFFT output
    res_time.to_vec()
}

// Cross correlation of 2 complex slices using RustFFT
// Assumes inputs powers of 2
// Naive: ifft(fft(a) * fft(b).conj())
#[allow(dead_code)]
fn xcor_rustfft(a: &[Complex32], b: &[Complex32]) -> Vec<Complex32> {

    // Sanity
    assert!(a.len() == b.len());

    // Allocations
    let n = a.len();
    let mut a_time = a.to_vec();
    let mut b_time = b.to_vec();
    let mut a_freq: Vec<Complex32> = vec![Complex::zero(); n];
    let mut b_freq: Vec<Complex32> = vec![Complex::zero(); n];
    let mut res_freq: Vec<Complex32> = vec![Complex::zero(); n];
    let mut res_time: Vec<Complex32> = vec![Complex::zero(); n];

    // Compute FFT(a), FFT(b)
    let mut forward_planner = FFTplanner::new(false);
    let fft = forward_planner.plan_fft(n);
    fft.process(&mut a_time, &mut a_freq);
    fft.process(&mut b_time, &mut b_freq);

    // Take complex conjugate of b
    for bin in b_freq.iter_mut() {
        *bin = bin.conj();
    }

    // Calculate a*b and normalize
    for (out, a, b) in izip!(res_freq.iter_mut(),
                             a_freq.iter(), b_freq.iter()) {
        *out = (a * b) / (n as f32);
    }

    // IFFT
    let mut inverse_planner = FFTplanner::new(true);
    let ifft = inverse_planner.plan_fft(n);
    ifft.process(&mut res_freq, &mut res_time);

    // Return IFFT output
    res_time
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
fn caf_surface(needle: &[Complex32], haystack: &[Complex32],
    freqs_hz: &[f32], fs: u32)
    -> Vec<Vec<Complex32>> {

    // Create our 2D surface
    let mut surface = Vec::new();

    // Run the cross correlation against the shifted ones
    for freq in freqs_hz.iter() {
        let shifted = apply_freq_shifts(needle, *freq, fs);
        let xcor_res = xcor_rustfft(&shifted, haystack);
        // let xcor_res = xcor_fftw(&shifted, haystack);
        surface.push(xcor_res);
    }

    // Return our CAF surface
    surface
}

fn main() {

    // Get signals 1 and 2 to compute the caf of
    let needle = to_complex_32(
        &read_file_f32("../data/burst_0000_raw.f32")
        .unwrap());
    let haystack = to_complex_32(
        &read_file_f32("../data/burst_0000_t+20.007860_f+88.202617.f32")
        .unwrap());

    // -50Hz to 50Hz, 0.5Hz step
    let mut shifts = Vec::new();
    for shift_millihz in (-50000..50000).step_by(500) {
        let shift = (shift_millihz as f32) / 1e3;
        shifts.push(shift);
    }

    // Get the CAF surface
    let surface = caf_surface(&needle, &haystack, &shifts, 48000);
    let freq_idx, samp_idx = find_2d_peak(surface);

    // Write our results out to a file for Python to parse
    // write_file_c32("results.c64", &surface[0]).unwrap();
}
