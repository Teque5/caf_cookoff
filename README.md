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
Time to compute a 400x8192 cross ambiguity surface using the "filterbank" CAF algorithm. I/O is all float64 and complex128.
#### Single Thread
| lang   | backend | accel        | Ryzen3900X 32G | i7-8550U 16G | ARM A57 4G |
|--------|---------|:------------:|:--------------:|:------------:|:----------:|
| go     | fftw    |              |     145 ms     |    178 ms    |      -     |
| rust   | fftw    |              |     109 ms     |    201 ms    |      -     |
| rust   | RustFFT |              |     177 ms     |    287 ms    |      -     |
| python | scipy   | +numba       |     164 ms     |    497 ms    |   2315 ms  |
| go     | go-dsp  |              |     827 ms     |    795 ms    |   2386 ms  |
| python | scipy   |              |    5630 ms     |   4336 ms    |  41700 ms  |

#### Multiple Threads
| lang   | backend | accel        | Ryzen3900X 32G | i7-8550U 16G | ARM A57 4G |
|--------|---------|--------------|:--------------:|:------------:|:----------:|
| rust   | RustFFT | +std::thread |      37 ms     |    147 ms    |      -     |
| go     | fftw    | +goroutines  |      41 ms     |     82 ms    |      -     |
| go     | go-dsp  | +goroutines  |      94 ms     |    208 ms    |   955 ms   |
| python | scipy   | +mp +numba   |     133 ms     |    161 ms    |   662 ms   |
| python | scipy   | +mp          |     599 ms     |   1634 ms    |  11299 ms  |

Implementation Notes
* go fftw implementation is not saving wisdom smartly. Also data still handled as complex128, but fftw wrapper only supports complex64 so i'm casting in and out during the cross-correlation.
* `numba` uses `@numba.njit` with type hinting.
* rust was not able to crosscompile the nightly bench for `aarch64` (armv8).
* go was not able to crosscompile fftw bindings for `aarch64` (armv8).

### Subjective Conclusions
|                         | python | rust |  go  |
|-------------------------|:------:|:----:|:----:|
| Min Time for Viable CAF |  1 hr  | 7 hr | 7 hr |
| Time for Parallel Ver   | 30 min | 2 hr | 2 hr |
| Performance             |  ★☆☆ | ★★★ | ★★☆ |
| Simplicity              |  ★★★ | ★★☆ | ★★☆ |
| Library Avail           |  ★★★ | ★★☆ | ★☆☆ |
| Cross-compilation       |  ☆☆☆ | ★★★ | ★★★ |

### Observations
* Numba is **amazing** and salvages Python's reputation in 2020
* Lack of benchmarking tools in Python sad 😿
* All three languages have excellent tooling
* Go and Python both have complex types, but rust uses a struct with two floats for complex.
* Go has fftw bindings or there is a fft library in go-dsp, but the latter isn't a full implementation and the former has quite a bit of complexity. I am disappointed such a basic tool isn't better integrated. The `go-dsp` implementation only supports `complex128` types, and the `fftw` wrapper only supports `complex64`, which is a real bummer. Additionally the `math` and `math/cmplx` libraries **only** support `complex128`. WHY!

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
#### Rust
```rust
fn apply_shift(ray: &[Complex64], freq_shift: f64, samp_rate: u32) -> Vec<Complex64> {
    // apply frequency shift
    let mut ray = ray.to_vec();
    let dt = 1.0 / (samp_rate as f64);
    let exp_common = Complex64::new(0.0, 2.0 * PI * dt * freq_shift);
    for (i, samp) in ray.iter_mut().enumerate() {
        let exp = Complex64::new(i as f64, 0.0) * exp_common;
        *samp *= Complex64::exp(&exp);
    }
    ray
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
