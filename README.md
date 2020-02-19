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

## Time to Minimum Viable Product
Both *go* and *rust* versions took similar time to construct initial filterbank versions, about 7 hours. This was mostly due to our non-familiarity with these languages. By comparison we had a tidy python version in under an hour. For numba acceleration double that.

## Compute Results
Time to compute a 400x8192 cross ambiguity surface.

| language | method         | backend      | precision | i7-8550U 16GB | Ryzen 9 3900X 32GB |
|----------|----------------|--------------|:---------:|:-------------:|:------------------:|
| go       | fb +goroutines | fftw         |     c64   |      82 ms    |                    |
| go       | fb             | fftw         |     c64   |     178 ms    |                    |
| go       | fb +goroutines | go-dsp       |    c128   |     208 ms    |         94ms       |
| rust     | fb             | fftw         |    c128   |     201 ms    |        109ms       |
| rust     | fb             | RustFFT      |    c128   |     287 ms    |        177ms       |
| python   | fb +numba1     | scipy.signal |    c128   |     622 ms    |        422ms       |
| python   | fb +numba0     | scipy.signal |    c128   |     696 ms    |                    |
| go       | fb             | go-dsp       |    c128   |     795 ms    |        827ms       |
| python   | fb             | scipy.signal |    c128   |    4336 ms    |                    |

Notes
* go fftw implementation is not saving wisdom smartly. data still handled as complex128, but fftw wrapper only supports complex64 so i'm casting in and out during the cross-correlation.
* fb means caf "filterbank" implementation
* numba0 is naive wrapping of functions with `@numba.jit`
* numba1 is `@numba.njit` with type hinting

## Run
### Requires
* python3
    * scipy
    * numpy
* Rust v1.14
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

## Observations
* Go has fftw bindings or there is a fft library in go-dsp, but the latter isn't a full implementation and the former has quite a bit of complexity. I am disappointed such a basic tool isn't better integrated. The `go-dsp` implementation only supports `complex128` types, and the `fftw` wrapper only supports `complex64`, which is a real bummer. Additionally the `math` and `math/cmplx` libraries **only** support `complex128`. WHY!
* Go and Python both have complex types, but rust uses a struct with two floats.

## References
1) S. Stein, Algorithms for ambiguity function processing,  IEEE Trans. Acoust., Speech, and Signal Processing, vol. ASSP-29, pp. 588 - 599, June 1981.
2) [Computing the Cross Ambiguity Function](http://ws.binghamton.edu/fowler/Fowler%20Personal%20Page/Publications_files/MS_Thesis_Chris_Yatrakis.pdf)
