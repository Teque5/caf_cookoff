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

## Results
TBD

## Run (wip)
### Requires
* python3
    * scipy
    * numpy
* Rust v1.14
* go v1.13
* GNU Radio if building signals
    * gr-sigmf

### Procedure
```bash
./build_all.sh
./benchmark.sh
```

## References
[Computing the Cross Ambiguity Function](http://ws.binghamton.edu/fowler/Fowler%20Personal%20Page/Publications_files/MS_Thesis_Chris_Yatrakis.pdf)

