package main

import (
  "log"
  "time"
)

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

	freqs_hz := arange(-100, 100, .5)
	// surf := amb_surf(apple_iq, banana_iq, freqs_hz, 48000)

	// start = time.Now()
	surf := amb_surf_concurrent(apple_iq, banana_iq, freqs_hz, 48000)
	// t = time.Now()
	// elapsed = t.Sub(start)
	// log.Printf("surf calculated in %s (concurrent)\n", elapsed)

	fdx, tdx, max := find_2d_peak(surf)
	t := time.Now()
	elapsed := t.Sub(start)
	log.Printf("surf calculated in %s\n", elapsed)
	log.Println("caf result:", len(apple_iq)-tdx, "samples", freqs_hz[fdx],"hz @ amb =",max)

	err = dump_surf("/tmp/derp", surf)
	if err != nil {
		log.Fatal(err)
	}
	log.Println("done")

}
