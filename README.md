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
| language | method         | backend      | precision | i7-8550U 16GB |
|----------|----------------|--------------|:---------:|:-------------:|
| rust     | fb             | fftw         |    c128   |     201 ms    |
| go       | fb +goroutines | go-dsp       |    c128   |     233 ms    |
| rust     | fb             | RustFFT      |    c128   |     287 ms    |
| python   | fb +numba1     | scipy.signal |    c128   |     622 ms    |
| python   | fb +numba0     | scipy.signal |    c128   |     696 ms    |
| go       | fb             | go-dsp       |    c128   |     877 ms    |
| python   | fb             | scipy.signal |    c128   |    4336 ms    |

* fb == caf filterbank implementation
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
go run .
go test -bench=. -benchtime=5
```
#### Python
```bash
cd caf_python
./caf.py
```

## Observations
* Go has fftw bindings or there is a fft library in go-dsp, but the latter isn't a full implementation and the former has quite a bit of complexity. I am disappointed such a basic tool isn't better integrated. The go-dsp implmentation only supports `complex128` types -> bizzare.
* Go and Python both have complex types, but rust uses a struct with two floats.

## References
1) S. Stein, Algorithms for ambiguity function processing,  IEEE Trans. Acoust., Speech, and Signal Processing, vol. ASSP-29, pp. 588 - 599, June 1981.
2) [Computing the Cross Ambiguity Function](http://ws.binghamton.edu/fowler/Fowler%20Personal%20Page/Publications_files/MS_Thesis_Chris_Yatrakis.pdf)
