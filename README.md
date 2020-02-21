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
Time to compute a 400x8192 cross ambiguity surface.
#### Single Thread
| language | method          | backend      | precision | i7-8550U 16GB | Ryzen3900X 32GB | ARM A57 4GB |
|----------|-----------------|--------------|:---------:|:-------------:|:---------------:|:-----------:|
| go       | fb              | fftw         |     c64   |     178 ms    |      145 ms     |      -      |
| rust     | fb              | fftw         |    c128   |     201 ms    |      109 ms     |      -      |
| rust     | fb              | RustFFT      |    c128   |     287 ms    |      177 ms     |      -      |
| python   | fb +numba       | scipy.signal |    c128   |     622 ms    |      422 ms     |      -      |
| go       | fb              | go-dsp       |    c128   |     795 ms    |      827 ms     |   2386 ms   |
| python   | fb              | scipy.signal |    c128   |    4336 ms    |                 |  41700 ms   |

#### Multiple Threads
| language | method          | backend      | precision | i7-8550U 16GB | Ryzen3900X 32GB | ARM A57 4GB |
|----------|-----------------|--------------|:---------:|:-------------:|:---------------:|:-----------:|
| rust     | fb +std::thread | RustFFT      |    c128   |     147 ms    |      37 ms      |      -      |
| python   | fb +mp +numba   | scipy.signal |    c128   |     181 ms    |      ? ms       |      -      |
| go       | fb +goroutines  | fftw         |     c64   |      82 ms    |      41 ms      |      -      |
| go       | fb +goroutines  | go-dsp       |    c128   |     208 ms    |      94 ms      |    955 ms   |
| python   | fb +mp          | scipy.signal |    c128   |    1711 ms    |      ? ms       |      -      |

Notes
* go fftw implementation is not saving wisdom smartly. Data still handled as complex128, but fftw wrapper only supports complex64 so i'm casting in and out during the cross-correlation.
* fb indiates the "filterbank" CAF algorithm
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
