package main

import (
	"fmt"
	// "io"
	"log"
	"os"
	// "sort"
	// "math"
	"encoding/binary"
	// "log"
	// "gopkg.in/yaml.v2"
)

// type Sample struct {
//
// }

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

	err = binary.Read(handle, binary.BigEndian, &ray)
	if err != nil {
		log.Fatalf("error reading floats %v", err)
	}
	return ray, nil
}

func f32_to_c64(ray_i []float32) (ray_iq []complex64) {
	// Convert real float32 to complex64
	len_ray := len(ray_i)
	ray_q := make([]float32, len_ray)
	ray_iq = make([]complex64, len_ray)
	for idx := 0; idx < len_ray; idx += 1 {
		ray_iq[idx] = complex(ray_i[idx], ray_q[idx])
	}
	return ray_iq
}

// func xcor(apple []complex64, banana []complex64) (corr []complex64) {
//
// }

func main() {
	data_path := "../data/"
	apple, err := load_f32(data_path + "burst_0000_raw.f32")
	if err != nil {
		log.Fatal(err)
	}
	apple_iq := f32_to_c64(apple)
	banana, err := load_f32(data_path + "burst_0000_t+20.007860_f+88.202617.f32")
	if err != nil {
		log.Fatal(err)
	}
	banana_iq := f32_to_c64(banana)
	fmt.Println(apple_iq)
	fmt.Println(banana_iq)
}

// handle_banana, err := os.Open(data_path + "burst_0000_raw.f32")
// if err != nil {
//   log.Fatalf("cannot open apple %v", err)
// }
