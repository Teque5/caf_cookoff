// Cross-correlation implementation using FFTW
// Assumes equal-length, power of 2 Complex64 slices
// in and returns their (equal length) complex
// cross-correlation
// Naive: ifft(fft(a) * fft(b).conj())

use itertools::izip;
use fftw::array::AlignedVec;
use fftw::plan::{C2CPlan, C2CPlan64};
use fftw::types::{Sign, Flag};
use num_complex::{Complex64};

#[allow(dead_code)]
pub struct Xcor {
    n: usize, // size of a, b, c
    // Aligned FFTW buffers
    a: AlignedVec<Complex64>,
    b: AlignedVec<Complex64>,
    c: AlignedVec<Complex64>,
    // Planners
    forward_planner: C2CPlan64,
    reverse_planner: C2CPlan64,
}

impl Xcor {

    // Constructor
    #[allow(dead_code)]
    pub fn new(n: usize) -> Self {

        // Create planners
        let fp = C2CPlan::aligned(
            &[n], Sign::Forward, Flag::MEASURE).unwrap();
        let rp = C2CPlan::aligned(
            &[n], Sign::Backward, Flag::MEASURE).unwrap();

        // Return new struct
        Xcor {
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
    pub fn run(&mut self, a: &[Complex64], b: &[Complex64]) -> Vec<Complex64> {

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

            *out = (a * b) / (self.n as f64);
        }

        // Calculate IFFT of product and return
        self.reverse_planner.c2c(&mut self.a, &mut self.b).unwrap();
        self.b.to_vec()
    }
}
