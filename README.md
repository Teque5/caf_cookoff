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
| python | scipy   | +numba       |     231 ms     |    622 ms    |   2315 ms  |
| go     | go-dsp  |              |     827 ms     |    795 ms    |   2386 ms  |
| python | scipy   |              |    5630 ms     |    4336 ms   |  41700 ms  |

#### Multiple Threads
| lang   | backend | accel        | Ryzen3900X 32G | i7-8550U 16G | ARM A57 4G |
|--------|---------|--------------|:--------------:|:------------:|:----------:|
| rust   | RustFFT | +std::thread |      37 ms     |    147 ms    |      -     |
| go     | fftw    | +goroutines  |      41 ms     |     82 ms    |      -     |
| go     | go-dsp  | +goroutines  |      94 ms     |    208 ms    |   955 ms   |
| python | scipy   | +mp +numba   |     134 ms     |    163 ms    |   662 ms   |
| python | scipy   | +mp          |     599 ms     |    1711 ms   |  11299 ms  |

Notes
* go fftw implementation is not saving wisdom smartly. Also data still handled as complex128, but fftw wrapper only supports complex64 so i'm casting in and out during the cross-correlation.
* `numba` uses `@numba.njit` with type hinting.
* `arm64` was unable to run numba due to LLVM issues and golang was not able to compile fftw bindings on `arm64`.

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
* Go has fftw bindings or there is a fft library in go-dsp, but the latter isn't a full implementation and the former has quite a bit of complexity. I am disappointed such a basic tool isn't better integrated. The `go-dsp` implementation only supports `complex128` types, and the `fftw` wrapper only supports `complex64`, which is a real bummer. Additionally the `math` and `math/cmplx` libraries **only** support `complex128`. WHY!
* Go and Python both have complex types, but rust uses a struct with two floats.

## Run
### Requires
* python3
    * scipy
    * numpy
* Rust v1.41
* go v1.13
* GNU Radio if building signals
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
#### Golang
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

## References
1) S. Stein, Algorithms for ambiguity function processing,  IEEE Trans. Acoust., Speech, and Signal Processing, vol. ASSP-29, pp. 588 - 599, June 1981.
2) [Computing the Cross Ambiguity Function](http://ws.binghamton.edu/fowler/Fowler%20Personal%20Page/Publications_files/MS_Thesis_Chris_Yatrakis.pdf)
