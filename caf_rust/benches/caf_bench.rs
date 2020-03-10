// Benchmarks to compare strategies to compar the CAF surface
// Could definitely be consolidated into functions at some point,
// but not worth it atm


#![feature(test)]
extern crate test;


#[cfg(test)]
mod caf_benches {
    use caf_rust::caf::{CafSurface,
        CafRustFFT, CafRustFFTThreads, CafRustFFTThreadpool, CafFFTW};
    use caf_rust::utils::read_file_c64;
    use test::{black_box, Bencher};

    #[bench]
    fn bench_rustfft(b: &mut Bencher) {
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

        b.iter(|| black_box({
            // Get the CAF surface
            let surface = CafRustFFT::caf_surface(&needle, &haystack, &shifts, 48000);
            CafRustFFT::find_peak(surface);
        }));
    }

    #[bench]
    fn bench_fftw(b: &mut Bencher) {
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

        b.iter(|| black_box({
            // Get the CAF surface
            let surface = CafFFTW::caf_surface(&needle, &haystack, &shifts, 48000);
            CafFFTW::find_peak(surface);
        }));
    }

    #[bench]
    fn bench_rustfft_threads(b: &mut Bencher) {
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

        b.iter(|| black_box({
            // Get the CAF surface
            let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
            CafRustFFTThreads::find_peak(surface);
        }));
    }

    #[bench]
    fn bench_rustfft_threadpool(b: &mut Bencher) {
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

        b.iter(|| black_box({
            // Get the CAF surface
            let surface = CafRustFFTThreadpool::caf_surface(&needle, &haystack, &shifts, 48000);
            CafRustFFTThreadpool::find_peak(surface);
        }));
    }

    #[bench]
    fn bench_apply_fdoa(b: &mut Bencher) {
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let freq_hz = 77.77;
        let samp_rate = 48000;
        b.iter(|| black_box({
            // Apply frequency offset
            CafRustFFT::apply_freq_shift(&needle, freq_hz, samp_rate);
        }));
    }
}
