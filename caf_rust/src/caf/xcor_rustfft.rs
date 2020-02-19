// Cross-correlation implementation using RustFFT
// Assumes equal-length, power of 2 Complex64 slices
// in and returns their (equal length) complex
// cross-correlation
// Naive: ifft(fft(a) * fft(b).conj())

use std::sync::Arc;

use itertools::izip;
use num_complex::{Complex64};
use rustfft::{FFTplanner, FFT};

#[allow(dead_code)]
pub struct Xcor {
    n: usize, // size of a, b, c
    // Preallocated buffers
    a: Vec<Complex64>,
    b: Vec<Complex64>,
    c: Vec<Complex64>,
    // Planners
    fft: Arc<dyn FFT<f64>>,
    ifft: Arc<dyn FFT<f64>>,
}

impl Xcor {

    // Constructor
    #[allow(dead_code)]
    pub fn new(n: usize) -> Self {

        // Create FFT plans
        let mut fp = FFTplanner::new(false);
        let fft = fp.plan_fft(n);
        let mut rp = FFTplanner::new(true);
        let ifft = rp.plan_fft(n);

        // Return new struct
        Xcor {
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
    pub fn run(&mut self, a: &[Complex64], b: &[Complex64]) -> Vec<Complex64> {

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

            *out = (a * b) / (self.n as f64);
        }

        // Calculate IFFT of product and return
        self.ifft.process(&mut self.a, &mut self.b);
        self.b.to_vec()
    }
}

// Implement the very cheap copy for FFT wisdom
impl Clone for Xcor {
    fn clone(&self) -> Xcor {
        Xcor {
            n: self.n,
            a: vec![Default::default(); self.n],
            b: vec![Default::default(); self.n],
            c: vec![Default::default(); self.n],
            fft: Arc::clone(&(self.fft)),
            ifft: Arc::clone(&(self.ifft)),
        }
    }
}
