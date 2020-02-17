package main

import (
	"encoding/binary"
	"fmt"
	"github.com/mjibson/go-dsp/fft"
	"log"
	"math"
	"math/cmplx"
	"os"
	"time"
)

func dump_surf(path string, surf [][]float64) (err error) {
	// dump 2d ambiguity surface to a file
	handle, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("cannot open handle %v", err)
	}
	defer handle.Close()
	for _, row := range surf {
		err = binary.Write(handle, binary.LittleEndian, &row)
		if err != nil {
			return
		}
	}
	fmt.Println("wrote surf to file", len(surf), len(surf[0]))
	return
}

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
	// create array for zero-padding

	zero_pad := make([]complex128, len_ray)
	if len(banana) != len_ray {
		panic("input arrays should be same size")
	}
	apple_padded := append(apple, zero_pad...)
	banana_padded := append(zero_pad, banana...)
	apple_fft := fft.FFT(apple_padded)
	banana_fft := fft.FFT(banana_padded)
	fleeb := make([]complex128, 2*len_ray)
	for idx, _ := range apple_fft {
		fleeb[idx] = apple_fft[idx] * cmplx.Conj(banana_fft[idx])
	}
	// fmt.Println(banana_fft)
	// corr := make([]complex128, len_ray)
	corr_abs = make([]float64, 2*len_ray)
	corr := fft.IFFT(fleeb)

	for idx, val := range corr {
		corr_abs[idx] = cmplx.Abs(val)
	}
	return
}

func apply_fdoa(ray []complex128, fdoa float64, samp_rate float64) (new_ray []complex128) {
	// Apply frequency shift
	precache := complex(0, -2*math.Pi*fdoa/samp_rate)
	new_ray = make([]complex128, len(ray))
	for idx, val := range ray {
		new_ray[idx] = val * cmplx.Exp(precache*complex(float64(idx), 0))
	}
	return
}

func amb_surf(needle []complex128, haystack []complex128, freqs_hz []float64, samp_rate float64) ([][]float64) {
	// Create cross ambiguity surface (filterbank method)
	len_ray := len(needle)
	len_freq := len(freqs_hz)
	shifted := make([]complex128, len_ray)
	// corr := make([]float64, len_ray)
	surf := make([][]float64, len_freq, len_ray)
	for fdx, freq_hz := range freqs_hz {
		shifted = apply_fdoa(needle, freq_hz, samp_rate)
		// corr =
		surf[fdx] = xcor(shifted, haystack)
	}
	return surf
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

	start := time.Now()
	data_path := "../data/"
	apple, err := load_c64(data_path + "chirp_4_raw.c64")
	if err != nil {
		log.Fatal(err)
	}
	apple_iq := c64_to_c128(apple)
	banana, err := load_c64(data_path + "chirp_4_T+70samp_F+82.89Hz.c64")
	if err != nil {
		log.Fatal(err)
	}
	banana_iq := c64_to_c128(banana)[0:4096]

	freqs_hz := arange(-100, 100, 0.5)
	surf := amb_surf(apple_iq, banana_iq, freqs_hz, 48000)
	fdx, tdx := find_2d_peak(surf)
	fmt.Println("result:", len(apple_iq)-tdx, "samples", freqs_hz[fdx],"hz")
	t := time.Now()
	elapsed := t.Sub(start)
	fmt.Println(elapsed)
	err = dump_surf("/tmp/derp", surf)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("done")

}
