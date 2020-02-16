package main

import (
	"fmt"
	// "io"
	"log"
	"math/cmplx"
	"os"
	"encoding/binary"
	"math"
	// "log"
	"github.com/mjibson/go-dsp/fft"
)

func load_c64(path string) (ray []complex64, err error) {
	// Load complex file
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
	ray = make([]complex64, sample_count)

	err = binary.Read(handle, binary.LittleEndian, &ray)
	if err != nil {
		log.Fatalf("error reading floats %v", err)
	}
	// log.Printf("loaded %f samples", sample_count)
	return ray, nil
}

func load_f32(path string) (ray []float32, err error) {
	// Load float file
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
	ray = make([]float32, sample_count)

	err = binary.Read(handle, binary.LittleEndian, &ray)
	if err != nil {
		log.Fatalf("error reading floats %v", err)
	}
	// log.Printf("loaded %f samples", sample_count)
	return ray, nil
}

func f32_to_c128(ray_i []float32) (ray_iq []complex128) {
	// Convert real float32 to complex128
	len_ray := len(ray_i)
	ray_iq = make([]complex128, len_ray)
	for idx := 0; idx < len_ray; idx += 1 {
		ray_iq[idx] = complex(float64(ray_i[idx]), 0)
	}
	return
}

func c64_to_c128(ray_i []complex64) (ray_iq []complex128) {
	// Convert real complex64 to complex128
	len_ray := len(ray_i)
	ray_iq = make([]complex128, len_ray)
	for idx := 0; idx < len_ray; idx += 1 {
		ray_iq[idx] = complex128(ray_i[idx])
	}
	return
}

func xcor(apple []complex128, banana []complex128) (corr_abs []float64) {
	// Standard crosscorrelation implementation
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
	corr := make([]complex128, len_ray)
	corr_abs = make([]float64, len_ray)
	corr = fft.IFFT(fleeb)

	for idx, val := range corr {
		corr_abs[idx] = cmplx.Abs(val)
	}
	return
}

func apply_fdoa(ray []complex128, fdoa float64, samp_rate float64) []complex128 {
	// Apply frequency shift
	len_ray := len(ray)
	precache := complex(0, -2*math.Pi*fdoa/samp_rate)
	for idx := 0; idx < len_ray; idx += 1 {
		ray[idx] = ray[idx] * cmplx.Exp(precache*complex(float64(idx), 0))
	}
	return ray
}

func amb_surf(needle []complex128, haystack []complex128, freqs_hz []float64, samp_rate float64) (surf [][]float64) {
	// Create cross ambiguity surface
	len_ray := len(needle)
	shifted := make([]complex128, len_ray)
	corr := make([]float64, len_ray)
	for _, freq_hz := range freqs_hz {
		shifted = apply_fdoa(needle, freq_hz, samp_rate)
		corr = xcor(shifted, haystack)
		surf = append(surf, corr)
	}
	return
}

func arange(start, stop, step float64) (rnge []float64) {
	N := int(math.Ceil((stop - start) / step))
	rnge = make([]float64, N, N)
	i := 0
	for x := start; x < stop; x += step {
		rnge[i] = x
		i += 1
	}
	return
}

func find_2d_peak(ray2d [][]float64) (best_idx int, best_jdx int) {
	// 2D Argmax
	var max float64 = 0
	for idx, row := range ray2d {
		for jdx, elem := range row {
			if elem > max {
				max = elem
				best_idx = idx
				best_jdx = jdx
				// fmt.Println("up", fmax, tmax)
			}
		}
	}
	return
}

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

	freqs_hz := arange(-100, 100, 0.5)
	surf := amb_surf(apple_iq, banana_iq, freqs_hz, 48000)
	fdx, tdx := find_2d_peak(surf)
	fmt.Println("out", 4096-tdx, freqs_hz[fdx])

}
