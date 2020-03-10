use std::io;
use std::io::prelude::*;
use std::fs::File;

use num_complex::Complex64;


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

