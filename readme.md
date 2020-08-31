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
