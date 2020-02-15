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

// func xcor()

func main() {
	data_path := "../data/"
  apple, err := load_f32(data_path + "burst_0000_raw.f32")
  if err != nil {
    log.Fatal(err)
  }
  banana, err := load_f32(data_path + "burst_0000_t+20.007860_f+88.202617.f32")
  if err != nil {
    log.Fatal(err)
  }
  fmt.Println(apple)
  fmt.Println(banana)
}


  // handle_banana, err := os.Open(data_path + "burst_0000_raw.f32")
  // if err != nil {
  //   log.Fatalf("cannot open apple %v", err)
  // }
