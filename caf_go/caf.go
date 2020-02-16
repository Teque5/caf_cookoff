package main

import (
	"fmt"
	// "io"
	"log"
	"math/cmplx"
	"os"
	// "sort"
	"math"
	"encoding/binary"
	// "log"
	// "gopkg.in/yaml.v2"
	"github.com/mjibson/go-dsp/fft"
)

func load_c64(path string) (ray []complex64, err error) {
  handle, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("cannot open handle %v", err)
	}
	defer handle.Close()
	info, err := handle.Stat()
	if err != nil {
		return nil, err
	}
	sample_count := info.Size() / 8
	fmt.Println(sample_count)
	ray = make([]complex64, sample_count)

	err = binary.Read(handle, binary.LittleEndian, &ray)
	if err != nil {
		log.Fatalf("error reading floats %v", err)
	}
	return ray, nil
}

func load_f32(path string) (ray []float32, err error) {
	handle, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("cannot open handle %v", err)
	}
	defer handle.Close()
	info, err := handle.Stat()
	if err != nil {
		return nil, err
	}
	sample_count := info.Size() / 4
	fmt.Println(sample_count)
	ray = make([]float32, sample_count)

	err = binary.Read(handle, binary.LittleEndian, &ray)
	if err != nil {
		log.Fatalf("error reading floats %v", err)
	}
	return ray, nil
}

func f32_to_c128(ray_i []float32) (ray_iq []complex128) {
	// Convert real float32 to complex64
	len_ray := len(ray_i)
	ray_iq = make([]complex128, len_ray)
	for idx := 0; idx < len_ray; idx += 1 {
		ray_iq[idx] = complex(float64(ray_i[idx]), 0)
	}
	return ray_iq
}

func c64_to_c128(ray_i []complex64) (ray_iq []complex128) {
	// Convert real float32 to complex64
	len_ray := len(ray_i)
	ray_iq = make([]complex128, len_ray)
	for idx := 0; idx < len_ray; idx += 1 {
		ray_iq[idx] = complex128(ray_i[idx])
	}
	return ray_iq
}

func xcor(apple []complex128, banana []complex128) (corr []complex128) {
	/*
	  Standard crosscorrelation implementation
	*/
	len_ray := len(apple)
  if len(banana) != len_ray {
    panic("input arrays should be same size")
  }
	apple_fft := fft.FFT(apple)
	banana_fft := fft.FFT(banana)
	fleeb := make([]complex128, len_ray)
	for idx := 0; idx < len_ray; idx += 1 {
		fleeb[idx] = apple_fft[idx] * cmplx.Conj(banana_fft[idx])
	}
	corr = fft.IFFT(fleeb)
	return corr
}

func apply_fdoa(ray []complex128, fdoa float64, samp_rate float64) []complex128 {
  len_ray := len(ray)
  pre := complex(0, -2 * math.Pi * fdoa / samp_rate)
  fmt.Println(pre)
  for idx := 0; idx <  len_ray; idx += 1 {
    ray[idx] = ray[idx] * cmplx.Exp(pre * complex(float64(idx), 0))
  }
  return ray
}

// func amb_surf(needle []complex128, haystack []complex128, freqs_hz []float64, samp_rate float64) (surf [][]complex128) {
//
// }

func main() {
	data_path := "../data/"
	apple, err := load_c64(data_path + "chirp_0_raw.c64")
	if err != nil {
		log.Fatal(err)
	}
	apple_iq := c64_to_c128(apple)
	banana, err := load_c64(data_path + "chirp_0_T+202samp_F.+69.25Hz.c64")
	if err != nil {
		log.Fatal(err)
	}
	banana_iq := c64_to_c128(banana)[0:4096]
	corr := xcor(apple_iq, banana_iq)
  trash := []complex128{1,2,3,4,5,6,7,8,9,10}
	fmt.Println(trash)
  fmt.Println(len(corr))
  z := apply_fdoa(corr, 200, 48e3)
  fmt.Println(z)
}
