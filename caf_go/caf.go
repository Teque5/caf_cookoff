package main

import (
	"encoding/binary"
	"fmt"
	"github.com/mjibson/go-dsp/fft"
	"log"
	"math"
	"math/cmplx"
	"os"
	"sync"
	"github.com/barnex/fftw"
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
	log.Printf("wrote (%vx%v) surf to file\n", len(surf), len(surf[0]))
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
	// log.Printf("loaded %v samples from %v\n", sample_count, path)
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
	// log.Printf("loaded %v samples from %v\n", sample_count, path)
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

func c128_to_c64(ray_i []complex128) (ray_iq []complex64) {
	len_ray := len(ray_i)
	ray_iq = make([]complex64, len_ray)
	for idx, val := range(ray_i) {
		ray_iq[idx] = complex64(val)
	}
	return
}

func xcor_fftw(apple []complex128, banana []complex128) (corr_abs []float64) {
	// Standard crosscorrelation implementation
	len_ray := len(apple)
	zero_pad := make([]complex64, len_ray)
	if len(banana) != len_ray {
		panic("input arrays should be same size")
	}
	apple_padded := append(c128_to_c64(apple), zero_pad...)
	banana_padded := append(zero_pad, c128_to_c64(banana)...)

	apple_fft := make([]complex64, 2*len_ray)
  plan := fftw.PlanC2C([]int{2*len_ray}, apple_padded, apple_fft, fftw.FORWARD, fftw.ESTIMATE)
  plan.Execute()
	banana_fft := make([]complex64, 2*len_ray)
	plan = fftw.PlanC2C([]int{2*len_ray}, banana_padded, banana_fft, fftw.FORWARD, fftw.ESTIMATE)
	plan.Execute()

	fleeb := make([]complex64, 2*len_ray)
	for idx, _ := range apple_fft {
		fleeb[idx] = apple_fft[idx] * complex64(cmplx.Conj(complex128(banana_fft[idx])))
	}

	corr := make([]complex64, 2*len_ray)
	plan = fftw.PlanC2C([]int{len_ray}, fleeb, corr, fftw.BACKWARD, fftw.ESTIMATE)
	plan.Execute()

	corr_abs = make([]float64, len(corr))
	for idx, val := range corr {
		corr_abs[idx] = cmplx.Abs(complex128(val))
	}
	return
}

func xcor(apple []complex128, banana []complex128) (corr_abs []float64) {
	// Standard crosscorrelation implementation
	len_ray := len(apple)
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
	corr := fft.IFFT(fleeb)
	corr_abs = make([]float64, len(corr))
	for idx, val := range corr {
		corr_abs[idx] = cmplx.Abs(val)
	}
	return
}

func apply_fdoa(ray []complex128, fdoa float64, samp_rate float64) (new_ray []complex128) {
	// Apply frequency shift
	precache := complex(0, 2*math.Pi*fdoa/samp_rate)
	new_ray = make([]complex128, len(ray))
	for idx, val := range ray {
		new_ray[idx] = val * cmplx.Exp(precache*complex(float64(idx), 0))
	}
	return
}

type amb_row struct {
	// used for storing xcors in channels
	xcor []float64
	fdx int
}

func surf_row(needle []complex128, haystack []complex128, freq_hz float64, samp_rate float64, fdx int, wg *sync.WaitGroup, c chan amb_row) {
	defer wg.Done()
	shifted := apply_fdoa(needle, freq_hz, samp_rate)
	fleeb := new(amb_row)
	fleeb.xcor = xcor_fftw(shifted, haystack)
	fleeb.fdx = fdx
  c <- *fleeb
}

func amb_surf_concurrent(needle []complex128, haystack []complex128, freqs_hz []float64, samp_rate float64) ([][]float64) {
	// Create cross ambiguity surface (filterbank method) + concurrency
	var wg sync.WaitGroup
	len_ray := len(needle)
	len_freq := len(freqs_hz)
	c := make(chan amb_row, len_freq)
	surf := make([][]float64, len_freq, len_ray)
	for fdx, freq_hz := range freqs_hz {
		wg.Add(1)
		go surf_row(needle, haystack, freq_hz, samp_rate, fdx, &wg, c)
	}
	wg.Wait()
	close(c)
	for fleeb := range(c) {
		surf[fleeb.fdx] = fleeb.xcor
	}
	return surf
}

func amb_surf(needle []complex128, haystack []complex128, freqs_hz []float64, samp_rate float64) ([][]float64) {
	// Create cross ambiguity surface (filterbank method)
	len_ray := len(needle)
	len_freq := len(freqs_hz)
	shifted := make([]complex128, len_ray)
	surf := make([][]float64, len_freq, len_ray)
	for fdx, freq_hz := range freqs_hz {
		shifted = apply_fdoa(needle, freq_hz, samp_rate)
		surf[fdx] = xcor_fftw(shifted, haystack)
	}
	return surf
}

func arange(start, stop, step float64) (rnge []float64) {
	// implement np.arange()
	for x := start; x < stop; x += step {
		rnge = append(rnge, x)
	}
	return
}

func find_2d_peak(ray2d [][]float64) (best_idx int, best_jdx int, max float64) {
	// 2D Argmax
	for idx, row := range ray2d {
		for jdx, elem := range row {
			if elem > max {
				max = elem
				best_idx = idx
				best_jdx = jdx
			}
		}
	}
	return
}
