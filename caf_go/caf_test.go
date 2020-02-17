package main

import (
  "testing"
)

func bench_load_c64(b *testing.B) {
  data_path := "../data/"

  for i := 0; i < b.N; i++ {
    load_c64(data_path + "chirp_8_raw.c64")
  }
}


func BenchmarkLoad(b *testing.B)  { bench_load_c64(b) }
