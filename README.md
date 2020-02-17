# CAF Cook-Off
Which cross ambiguity function will win? Implementations in *Rust*, *Go*, and *Python*.


## OV-1
We are trying to answer the following questions:
* Time required for basic implementation in each language?
* Time required for parralelized implementation (goroutines, threads, async)?
* Throughput performance?
* Cross-compilation?
* How simple is each implementation?

## Hypothesis
Teque5 predicts that *go* and *rust* will produce the fastest implementations, but go will have the simplist parralelized version.

## Time to Minimum Viable Product
Both *go* and *rust* versions took similar time to construct initial filterbank versions, about 7 hours.

## Compute Results
TBD

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
#### Rust
```bash
cd caf_rust
cargo run
cargo +nightly bench
```
#### Golang
```bash
cd caf_go
go run .
go test -bench=.
```

## Observations
* Go has fftw bindings or there is a fft library in go-dsp, but the latter isn't a full implementation and the former has quite a bit of complexity. I am disappointed such a basic tool isn't better integrated. The go-dsp implmentation on supports `complex128` types -> bizzare.
* Go and Python both have complex types, but rust uses a struct with two floats.

## References
1) S. Stein, Algorithms for ambiguity function processing,  IEEE Trans. Acoust., Speech, and Signal Processing, vol. ASSP-29, pp. 588 - 599, June 1981.
2) [Computing the Cross Ambiguity Function](http://ws.binghamton.edu/fowler/Fowler%20Personal%20Page/Publications_files/MS_Thesis_Chris_Yatrakis.pdf)
