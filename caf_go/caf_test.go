package main

import (
	"testing"
)

func bench_amb_surf(b *testing.B) {
	data_path := "../data/"
	apple, err := load_c64(data_path + "chirp_4_raw.c64")
	if err != nil {
		panic(err)
	}
	apple_iq := c64_to_c128(apple)
	banana, err := load_c64(data_path + "chirp_4_T+70samp_F+82.89Hz.c64")
  if err != nil {
		panic(err)
	}
	banana_iq := c64_to_c128(banana)[0:4096]

	freqs_hz := arange(-100, 100, .5)

	// Bench will increase N until a sufficient # of samples have been collected
	for i := 0; i < b.N; i++ {
		amb_surf(apple_iq, banana_iq, freqs_hz, 48000)
	}

}

func bench_amb_surf_concurrent(b *testing.B) {
	data_path := "../data/"
	apple, err := load_c64(data_path + "chirp_4_raw.c64")
	if err != nil {
		panic(err)
	}
	apple_iq := c64_to_c128(apple)
	banana, err := load_c64(data_path + "chirp_4_T+70samp_F+82.89Hz.c64")
  if err != nil {
		panic(err)
	}
	banana_iq := c64_to_c128(banana)[0:4096]

	freqs_hz := arange(-100, 100, .5)

	// Bench will increase N until a sufficient # of samples have been collected
	for i := 0; i < b.N; i++ {
		amb_surf_concurrent(apple_iq, banana_iq, freqs_hz, 48000)
	}

}

func bench_apply_fdoa(b *testing.B) {
	data_path := "../data/"
	apple, err := load_c64(data_path + "chirp_4_raw.c64")
	if err != nil {
		panic(err)
	}
	apple_iq := c64_to_c128(apple)
	var freq_hz float64 = 77.77
	var samp_rate float64 = 48000
	// Bench will increase N until a sufficient # of samples have been collected
	for i := 0; i < b.N; i++ {
		apply_fdoa(apple_iq, freq_hz, samp_rate)
	}
}

// runners must start with "Benchmark"
func BenchmarkSurf(b *testing.B) { bench_amb_surf(b) }
func BenchmarkSurfConcurrent(b *testing.B) { bench_amb_surf_concurrent(b) }
func BenchmarkApplyFDOA(b *testing.B) { bench_apply_fdoa(b) }
