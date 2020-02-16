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

func xcor(apple []complex128, banana []complex128) (corr []complex128) {
	/*
	  Standard crosscorrelation implementation
	*/
	len_ray := len(apple)
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

func main() {
	data_path := "../data/"
	apple, err := load_f32(data_path + "burst_0000_raw.f32")
	if err != nil {
		log.Fatal(err)
	}
	apple_iq := f32_to_c128(apple)
	banana, err := load_f32(data_path + "burst_0000_t+20.007860_f+88.202617.f32")
	if err != nil {
		log.Fatal(err)
	}
	banana_iq := f32_to_c128(banana)
	corr := xcor(apple_iq, banana_iq)
  trash := []complex128{1,2,3,4,5,6,7,8,9,10}
	fmt.Println(trash)
  fmt.Println(len(corr))
  z := apply_fdoa(trash, 200, 48e3)
  fmt.Println(z)


}
