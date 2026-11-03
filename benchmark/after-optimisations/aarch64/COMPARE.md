# Baseline → after the optimisation PRs (aarch64, matched cases only)

| class | kind | type | cases | FFTA speedup geomean (min–max) | FFTA/FFTW before → after | max bytes/exec before → after | plan time geomean before → after |
|:--|:--|:--|--:|--:|--:|:--|:--|
| 1d/awkward | fft | Float32 | 29 | 2.69× (1.96–3.66) | 7.45× → 2.75× | 48 MiB → 0 | 0.6 µs → 1.2 ms |
| 1d/awkward | fft | Float64 | 29 | 2.43× (1.66–3.59) | 5.63× → 2.29× | 96 MiB → 0 | 0.7 µs → 1.6 ms |
| 1d/awkward | rfft | Float32 | 29 | 2.71× (1.92–3.97) | 7.73× → 2.85× | 33 MiB → 0 | 0.4 µs → 897.1 µs |
| 1d/awkward | rfft | Float64 | 29 | 2.48× (1.75–3.96) | 6.17× → 2.43× | 66 MiB → 0 | 0.4 µs → 1.2 ms |
| 1d/pow2 | fft | Float32 | 20 | 2.26× (1.46–4.31) | 3.77× → 1.68× | 0 → 0 | 0.1 µs → 39.4 µs |
| 1d/pow2 | fft | Float64 | 47 | 2.04× (1.25–4.70) | 3.68× → 1.80× | 0 → 0 | 0.1 µs → 109.6 µs |
| 1d/pow2 | rfft | Float32 | 20 | 1.94× (1.29–2.67) | 4.47× → 2.30× | 16 MiB → 0 | 0.1 µs → 22.4 µs |
| 1d/pow2 | rfft | Float64 | 20 | 1.92× (1.29–2.79) | 3.29× → 1.74× | 32 MiB → 0 | 0.1 µs → 30.0 µs |
| 1d/prime | fft | Float32 | 19 | 3.30× (1.78–8.57) | 7.54× → 2.27× | 96 MiB → 0 | 0.1 µs → 127.1 µs |
| 1d/prime | fft | Float64 | 19 | 3.19× (1.85–6.57) | 5.81× → 1.80× | 192 MiB → 0 | 0.1 µs → 165.5 µs |
| 1d/prime | rfft | Float32 | 19 | 3.23× (2.03–4.74) | 6.57× → 2.01× | 116 MiB → 0 | 0.1 µs → 129.7 µs |
| 1d/prime | rfft | Float64 | 19 | 2.90× (1.47–4.75) | 6.25× → 2.12× | 232 MiB → 0 | 0.1 µs → 171.0 µs |
| 1d/smooth | fft | Float32 | 24 | 2.37× (1.60–3.09) | 10.18× → 4.28× | 0 → 0 | 0.8 µs → 73.8 µs |
| 1d/smooth | fft | Float64 | 24 | 2.65× (1.64–3.97) | 7.33× → 2.79× | 0 → 0 | 0.9 µs → 97.7 µs |
| 1d/smooth | rfft | Float32 | 24 | 2.27× (1.53–3.13) | 10.47× → 4.52× | 16 MiB → 0 | 0.8 µs → 55.8 µs |
| 1d/smooth | rfft | Float64 | 24 | 2.52× (1.49–3.70) | 8.66× → 3.49× | 32 MiB → 0 | 0.8 µs → 73.5 µs |
| 2d | fft | Float32 | 9 | 2.41× (1.53–3.75) | 6.06× → 2.48× | 32 KiB → 0 | 0.2 µs → 2.4 µs |
| 2d | fft | Float64 | 18 | 3.47× (1.45–10.26) | 6.48× → 1.85× | 9 MiB → 8 KiB | 0.2 µs → 9.2 µs |
| 2d | rfft | Float32 | 9 | 4.45× (3.21–5.67) | 13.23× → 2.92× | 80 MiB → 0 | 0.2 µs → 2.2 µs |
| 2d | rfft | Float64 | 14 | 5.01× (3.30–7.60) | 11.83× → 2.33× | 160 MiB → 0 | 0.2 µs → 6.9 µs |
| 3d | fft | Float32 | 5 | 2.69× (2.39–3.36) | 11.63× → 4.31× | 2 KiB → 0 | 0.3 µs → 1.3 µs |
| 3d | fft | Float64 | 5 | 2.62× (2.11–3.52) | 7.56× → 2.96× | 4 KiB → 0 | 0.3 µs → 1.6 µs |
| batched_dim1 | fft | Float32 | 6 | 2.04× (1.66–2.54) | 5.07× → 2.43× | 0 → 0 | 0.1 µs → 11.7 µs |
| batched_dim1 | fft | Float64 | 10 | 3.75× (1.42–11.87) | 6.17× → 1.57× | 0 → 4 KiB | 0.1 µs → 27.8 µs |
| batched_dim1 | rfft | Float32 | 6 | 2.56× (2.15–3.17) | 6.39× → 2.40× | 32 MiB → 0 | 0.1 µs → 6.5 µs |
| batched_dim1 | rfft | Float64 | 9 | 2.58× (2.12–3.26) | 4.90× → 1.69× | 65 MiB → 0 | 0.1 µs → 7.4 µs |
| batched_dim2 | fft | Float32 | 6 | 1.61× (1.15–2.19) | 2.90× → 1.66× | 0 → 0 | 0.1 µs → 11.8 µs |
| batched_dim2 | fft | Float64 | 6 | 1.53× (1.13–2.19) | 2.24× → 1.42× | 0 → 0 | 0.1 µs → 15.8 µs |
| batched_dim2 | rfft | Float32 | 6 | 2.10× (1.45–3.06) | 3.54× → 1.64× | 32 MiB → 0 | 0.1 µs → 6.7 µs |
| batched_dim2 | rfft | Float64 | 6 | 1.93× (1.25–3.17) | 3.05× → 1.51× | 65 MiB → 0 | 0.1 µs → 9.1 µs |

Largest slowdowns / speedups (planned execution, FFTA before → after; FFTW for reference):
- 1.13× — fft Float64 64×65536 dims=(2,): 324.29 ms → 288.25 ms (FFTW 285.84 ms)
- 1.15× — fft Float32 64×65536 dims=(2,): 279.49 ms → 242.85 ms (FFTW 209.93 ms)
- 1.24× — fft Float64 64×16384 dims=(2,): 62.99 ms → 50.61 ms (FFTW 39.57 ms)
- 1.25× — fft Float64 1048576 dims=(1,): 58.40 ms → 46.63 ms (FFTW 6.90 ms)
- 9.21× — fft Float64 2048×2048 dims=(1, 2): 334.95 ms → 36.36 ms (FFTW 26.05 ms)
- 9.43× — fft Float64 1024×64 dims=(1,): 1.37 ms → 145.4 µs (FFTW 76.1 µs)
- 10.26× — fft Float64 512×512 dims=(1, 2): 16.61 ms → 1.62 ms (FFTW 933.0 µs)
- 10.77× — fft Float64 4096×64 dims=(1,): 6.83 ms → 634.0 µs (FFTW 385.2 µs)
- 11.51× — fft Float64 16384×64 dims=(1,): 37.12 ms → 3.22 ms (FFTW 2.16 ms)
- 11.87× — fft Float64 65536×64 dims=(1,): 187.55 ms → 15.80 ms (FFTW 13.83 ms)

510 matched cases; geometric-mean speedup 2.57×; 0 cases slower by >5%.
