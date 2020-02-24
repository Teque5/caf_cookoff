#![feature(test)]
extern crate test;

#[cfg(test)]
mod caf_benches {
    use caf_rust::caf::{read_file_c64, CafSurface, CafRustFFTThreads};
    use test::{black_box, Bencher};

    #[bench]
    fn bench_rustfft_threads_chirp0(b: &mut Bencher) {
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
    fn bench_apply_fdoa(b: &mut Bencher) {
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let freq_hz = 77.77;
        let samp_rate = 48000;
        b.iter(|| black_box({
            // Apply frequency offset
            let shifted = apply_freq_shifts(&needle, freq_hz, samp_rate);
        }));
    }
}
