# JONSWAP 

Calculates a JONSWAP ocean wave spectrum and optionally a realisation in time domain.

# Usage

```
 jonswap hs tp [-n num_harmonics] [-g gamma] [-s seed] [-d duration] [-t timestep] [-v] [-h]

  hs: significant wave height in meters.
  tp: spectral peak period in seconds
  -n: number of harmonics used to discretise the spectrum.
  -g: value for gamma. If ommited, DNV is used.
  -s: seed number for phase randomisation.
  -d: duration - timetrace will be shown.
  -t: time step. default is 0.1 seconds.
  -v: verbose output
  -h: hide spectrum
```

# Compile

No special dependencies required for the library and command line interface.

## Main command line interface

`make`

## Mockup multithread performance implementation

`make jonswap_multithread`

## Plots

For the plots, `gnuplot_i` is required.
Copy /symlink `gnuplot_i.c` and `gnuplot_i.h` to the current directory

`make plots` 


