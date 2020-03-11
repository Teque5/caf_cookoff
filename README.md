# CAF Cook-Off
In our search for the language of the future, we decided to try and write a common dsp function in a few select languages. Which cross ambiguity function implementation will win? Implementations in **Rust**, **Go**, and **Python**.

## OV-1
We are trying to answer the following questions:
* Time required for basic implementation in each language?
* Time required for parralelized implementation (goroutines, threads, async)?
* Throughput performance?
* Cross-compilation?
* How simple is each implementation?

## CAF?
A cross ambiguity function (CAF) is a method of comparing complex waveforms to determine time and frequency offset.
![Signals Under Test](/docs/s0s1-time.png)
![CAF Surface](/docs/s0s1-caf.png)

## Hypothesis
Teque5 predicts that *go* and *rust* will produce the fastest implementations, but go will have the simplist parralelized version.

## Results
### Benchmarks
Time to compute a 400x8192 cross ambiguity surface using the "filterbank" CAF algorithm. I/O is all float64 and complex128 unless otherwise noted.
#### Single Thread
| lang   | backend | accel        | R9-3900X 32G | W-2135 256G | i7-8550U 16G | ARM A57 4G |
|--------|---------|:------------:|:------------:|:-----------:|:------------:|:----------:|
| rust   | fftw    |              |    109 ms    |    158 ms   |    201 ms    |      -     |
| go     | fftw*   | *c64 FFT     |    145 ms    |    182 ms   |    178 ms    |      -     |
| rust   | RustFFT |              |    177 ms    |    199 ms   |    287 ms    |      -     |
| python | scipy   | +numba       |    164 ms    |    476 ms   |    497 ms    |   2315 ms  |
| go     | go-dsp  |              |    827 ms    |    616 ms   |    795 ms    |   2386 ms  |
| python | scipy   |              |   5630 ms    |   3828 ms   |   4336 ms    |  41700 ms  |

#### Multiple Threads
| lang   | backend | accel        | R9-3900X 32G | W-3125 256G | i7-8550U 16G | ARM A57 4G |
|--------|---------|--------------|:------------:|:-----------:|:------------:|:----------:|
| rust   | RustFFT | +threadpool  |     28 ms    |     39 ms   |       -      |      -     |
| go     | fftw*   | +goroutines  |     41 ms    |     58 ms   |     82 ms    |      -     |
| rust   | RustFFT | +std::thread |     26 ms    |     58 ms   |    133 ms    |      -     |
| go     | go-dsp  | +goroutines  |     94 ms    |    106 ms   |    208 ms    |   955 ms   |
| python | scipy   | +mp +numba   |    133 ms    |    145 ms   |    161 ms    |   662 ms   |
| python | scipy   | +mp          |    599 ms    |    884 ms   |   1634 ms    |  11299 ms  |

Implementation Notes
* go fftw implementation is not saving wisdom smartly. Also data still handled as complex128, but fftw wrapper only supports complex64 so i'm casting in and out during the cross-correlation.
* `numba` uses `@numba.njit` with type hinting.
* rust was not able to crosscompile the nightly bench for `aarch64` (armv8).
* go was not able to crosscompile fftw bindings for `aarch64` (armv8).
* go without goroutines had to explicitly specify GOMAXPROCS=1 on the W-3125 (and only the W-3125).
  Failing to specify a single thread caused weird scheduling, leading to up to *3x slower performance*.
  This is particularly bizzare, as all x86 machines were running Ubuntu 18.04 with the same Go version.
* A multithreaded FFTW implementation was not attempted in Rust. Unlike RustFFT, the FFTW wrapper wasn't
  very explicit about how it handled atomic operations, if at all.

### Subjective Conclusions
|                         | python | rust  |  go   |
|-------------------------|:------:|:-----:|:-----:|
| Min Time for Viable CAF |  1 hr  | 7 hrs | 7 hrs |
| Time for Parallel Ver   | 30 min | 2 hrs | 2 hrs |
| Performance             |   ★☆☆  |  ★★★  |  ★★☆  |
| Simplicity              |   ★★★  |  ★★☆  |  ★★☆  |
| Library Avail           |   ★★★  |  ★★☆  |  ★☆☆  |
| Cross-compilation       |   ☆☆☆  |  ★★☆  |  ★★☆  |

### Observations
* Numba is **amazing** and salvages Python's reputation in 2020
* Lack of benchmarking tools in Python is quite sad 😿
* All three languages have excellent tooling
* Go and Python both have native complex types, but rust relies on a (popular) external library for support.
* Go and Rust MVPs were of comparable complexity
* Both Rust and Go threading are miles ahead of C/C++
* goroutines (green threads) end up being more "plug-and-play" than Rust's std::thread (os threads)
* Go has fftw bindings or there is a fft library in go-dsp, but the latter isn't a full implementation and the former has quite a bit of complexity. I am disappointed such a basic tool isn't better integrated. The `go-dsp` implementation only supports `complex128` types, and the `fftw` wrapper only supports `complex64`, which is a real bummer. Additionally the `math` and `math/cmplx` libraries **only** support `complex128`. WHY!
* The authors of the Rust and Go versions both had rated a 3/10 familiarity with the language before starting. Without
  some initial exposure to core Rust concepts, the MVP for Rust would have taken significantly longer.

## Run
### Requires
* python3
    * scipy
    * numpy
* Rust v1.41
* go v1.13
* GNU Radio if using grc
    * gr-sigmf

### Procedure
#### Data Generation
Install numpy, scipy for Python 3
```bash
cd utils
python3 ./generate.py
```
#### Rust
Install `rustup` from [rustup.rs](https://rustup.rs/)
```bash
rustup install nightly
cd caf_rust
cargo run
cargo test
cargo +nightly bench
```
#### Go
Install `go` from the [official downloads](https://golang.org/doc/install)
```bash
cd caf_go
go get github.com/mjibson/go-dsp/fft
go run .
go test -bench=. -benchtime=5
```
#### Python
```bash
cd caf_python
./caf.py
```
## Code Comparison
Implementations of the frequency shift function.

|             | rust   |   go   | numba  | python   |
|-------------|:------:|:------:|:------:|:--------:|
| apply_shift | 120 μs | 137 μs | 158 μs | 10300 μs |
|             |    1x  |  1.14x |  1.31x |    85x   |

#### Rust
```rust
fn apply_shift(ray: &[Complex64], freq_shift: f64, samp_rate: f64) -> Vec<Complex64> {
    // apply frequency shift
    let mut new_ray = ray.to_vec();
    let precache = Complex64::new(0.0, 2.0*PI*freq_shift/samp_rate);
    for (idx, val) in ray.iter_mut().enumerate() {
        let exp = Complex64::new(i as f64, 0.0) * exp_common;
        new_ray[idx] = val * Complex64::exp(&exp);
    }
    new_ray
}
```
#### Go
```go
func apply_shift(ray []complex128, freq_shift float64, samp_rate float64) (new_ray []complex128) {
	// apply frequency shift
	precache := complex(0, 2*math.Pi*freq_shift/samp_rate)
	new_ray = make([]complex128, len(ray))
	for idx, val := range ray {
		new_ray[idx] = val * cmplx.Exp(precache*complex(float64(idx), 0))
	}
	return
}
```
#### Python (+Numba)
```python
@numba.njit
def apply_shift(ray: np.ndarray, freq_shift: np.float64, samp_rate: np.float64) -> np.ndarray:
    '''apply frequency shift'''
    precache = 2j * np.pi * freq_shift / samp_rate
    new_ray = np.empty_like(ray)
    for idx, val in enumerate(ray):
        new_ray[idx] = val * np.exp(precache * idx)
    return new_ray
```

## References
1) S. Stein, Algorithms for ambiguity function processing,  IEEE Trans. Acoust., Speech, and Signal Processing, vol. ASSP-29, pp. 588 - 599, June 1981.
2) [Computing the Cross Ambiguity Function](http://ws.binghamton.edu/fowler/Fowler%20Personal%20Page/Publications_files/MS_Thesis_Chris_Yatrakis.pdf)
