// Integration tests for Chirp 0-9
// Also runs chrip0 through all xcor strategies to check
// that all strategies give a right answer

extern crate caf_rust;

#[cfg(test)]
mod tests {

    use num_complex::Complex64;
    use caf_rust::caf::*;
    use caf_rust::utils::read_file_c64;

    #[test]
    fn test_rustfft_chirp0() {
        // Read Chirp 0 reference and modified files
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_0_T+202samp_F+69.25Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFT::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFT::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 69.25);
        assert_eq!(samp_idx, 202);
    }

    #[test]
    fn test_rustfft_iter_chirp0() {
        // Read Chirp 0 reference and modified files
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_0_T+202samp_F+69.25Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTIter::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTIter::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 69.25);
        assert_eq!(samp_idx, 202);
    }

    #[test]
    fn test_rustfft_rayon_chirp0() {
        // Read Chirp 0 reference and modified files
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_0_T+202samp_F+69.25Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTRayon::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTRayon::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 69.25);
        assert_eq!(samp_idx, 202);
    }

    #[test]
    fn test_rustfft_rayon_iter_chirp0() {
        // Read Chirp 0 reference and modified files
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_0_T+202samp_F+69.25Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTIterRayon::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTIterRayon::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 69.25);
        assert_eq!(samp_idx, 202);
    }

    #[test]
    fn test_fftw_chirp0() {
        // Read Chirp 0 reference and modified files
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_0_T+202samp_F+69.25Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafFFTW::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafFFTW::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 69.25);
        assert_eq!(samp_idx, 202);
    }

    #[test]
    fn test_rustfft_threads_chirp0() {
        // Read Chirp 0 reference and modified files
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_0_T+202samp_F+69.25Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 69.25);
        assert_eq!(samp_idx, 202);
    }

    #[test]
    fn test_rustfft_threadpool_chirp0() {
        // Read Chirp 0 reference and modified files
        let needle = read_file_c64("../data/chirp_0_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_0_T+202samp_F+69.25Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTThreadpool::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreadpool::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 69.25);
        assert_eq!(samp_idx, 202);
    }

    #[test]
    fn test_rustfft_threads_chirp1() {
        // Read Chirp 1 reference and modified files
        let needle = read_file_c64("../data/chirp_1_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_1_T+78samp_F+35.99Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -50Hz to 50Hz, 1Hz step
        let shifts = gen_float_shifts(-50.0, 50.0, 1.0);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 36.0);
        assert_eq!(samp_idx, 78);
    }

    #[test]
    fn test_rustfft_threads_chirp2() {
        // Read Chirp 2 reference and modified files
        let needle = read_file_c64("../data/chirp_2_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_2_T+169samp_F+32.16Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // 30Hz to 35Hz, 0.05Hz step
        let shifts = gen_float_shifts(30.0, 35.0, 0.05);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 32.15);
        assert_eq!(samp_idx, 169);
    }

    #[test]
    fn test_rustfft_threads_chirp3() {
        // Read Chirp 3 reference and modified files
        let needle = read_file_c64("../data/chirp_3_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_3_T+151samp_F-76.22Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, -76.25);
        assert_eq!(samp_idx, 151);
    }

    #[test]
    fn test_rustfft_threads_chirp4() {
        // Read Chirp 4 reference and modified files
        let needle = read_file_c64("../data/chirp_4_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_4_T+70samp_F+82.89Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // 80Hz to 100Hz, 0.1Hz step
        let shifts = gen_float_shifts(80.0, 100.0, 0.1);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 82.9);
        assert_eq!(samp_idx, 70);
    }

    #[test]
    fn test_rustfft_threads_chirp5() {
        // Read Chirp 5 reference and modified files
        let needle = read_file_c64("../data/chirp_5_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_5_T+177samp_F-92.72Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, -92.75);
        assert_eq!(samp_idx, 177);
    }

    #[test]
    fn test_rustfft_threads_chirp6() {
        // Read Chirp 6 reference and modified files
        let needle = read_file_c64("../data/chirp_6_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_6_T+15samp_F-49.69Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, -49.75);
        assert_eq!(samp_idx, 15);
    }

    #[test]
    fn test_rustfft_threads_chirp7() {
        // Read Chirp 7 reference and modified files
        let needle = read_file_c64("../data/chirp_7_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_7_T+84samp_F+68.26Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, 68.25);
        assert_eq!(samp_idx, 84);
    }

    #[test]
    fn test_rustfft_threads_chirp8() {
        // Read Chirp 8 reference and modified files
        let needle = read_file_c64("../data/chirp_8_raw.c64").unwrap();
        let haystack = read_file_c64("../data/chirp_8_T+80samp_F-46.28Hz.c64").unwrap();
        let haystack = &haystack[..needle.len()];

        // -100Hz to 100Hz, 0.25Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.25);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

        // Confirm correct results
        assert_eq!(freq, -46.25);
        assert_eq!(samp_idx, 80);
    }

    #[test]
    fn test_rustfft_threads_chirp9() {
        // Read Chirp 9 reference and modified files
        let (needle, haystack) = load_files(
            "../data/chirp_9_raw.c64",
            "../data/chirp_9_T+176samp_F+61.49Hz.c64");

        // -100Hz to 100Hz, 0.5Hz step
        let shifts = gen_float_shifts(-100.0, 100.0, 0.5);

        // Get the CAF estimates
        let surface = CafRustFFTThreads::caf_surface(&needle, &haystack, &shifts, 48000);
        let (freq, samp_idx) = CafRustFFTThreads::find_peak(surface);

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
