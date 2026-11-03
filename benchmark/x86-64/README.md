# x86-64 results for the comprehensive benchmark suite

Companion to the aarch64 run in `benchmark/REPORT.md`. Same `suite.jl` and
`report.jl` (bench branch @ 73b190c for every sweep, so before/after pairs are
methodologically identical), same knobs (`--seconds 0.5`, `--maxlog2 22`, `-t 8`),
same FFTA baseline (0.3.1, `src/` byte-identical to `main` @ 7aeb327).

## Machine

Intel Core Ultra 7 165H (Meteor Lake-H), 6 P-cores (×2 SMT) + 8 E-cores + 2 LP-E,
**AVX2 + FMA, no AVX-512** (LLVM target `alderlake`), L3 24 MiB, 31 GiB RAM,
Linux under WSL2, Julia 1.12.6, FFTW 3.3.11, provider `fftw` (not MKL).

Because this part has no AVX-512, FFTW holds 2 `ComplexF64` per vector register
rather than 4, so the FFTA/FFTW ratios here are a **lower bound** on the x86-64
gap, not a typical value for server parts.

## Headline

FFTA/FFTW execution ratio, single-threaded, FFTW `ESTIMATE`, baseline → the full
optimisation stack (lower is better):

| type | pow2 | smooth | prime | awkward | 2D | 3D | batched d=1 | batched d=2 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| ComplexF64 fft | 3.09 → **1.02** | 13.14 → 3.95 | 7.58 → **1.31** | 5.98 → 1.65 | 5.77 → 1.59 | 6.98 → 2.65 | 3.18 → **1.15** | 1.89 → 1.31 |
| ComplexF32 fft | 4.68 → 1.75 | 16.80 → 6.59 | 10.06 → 1.89 | 8.66 → 2.54 | 5.22 → 1.77 | 11.61 → 3.64 | 6.47 → 2.35 | 2.22 → 1.28 |
| Float64 rfft | 3.83 → 1.67 | 13.67 → 4.32 | 6.12 → **1.32** | 6.09 → 1.67 | 14.86 → **2.10** | ✗ → 3.65 | 6.86 → 1.96 | 2.41 → **1.10** |
| Float32 rfft | 5.54 → 2.56 | 14.75 → 6.21 | 6.21 → 1.44 | 8.37 → 2.38 | 14.13 → 2.76 | ✗ → 4.00 | 7.86 → 2.77 | 2.62 → 1.19 |

Stack vs baseline: **3.24× geomean over 510 matched cases**, best 18.82×,
worst 0.91×, one row below 0.95 (inside this machine's noise band; see below).

Against FFTW `MEASURE` rather than `ESTIMATE`, pow2 `ComplexF64` is 1.5–1.8× at
n = 2^14–2^19 and 2.7–3.2× at n ≥ 2^20 — quote the `ESTIMATE` figure only with
that qualifier attached.

## What the second architecture changed

1. **The Bluestein constants did not transfer.** The 3-smooth pad cost factor
   measures 2.01–3.11 here against 1.28–2.26 on aarch64; the DFT/Bluestein
   crossover is n = 23 here against 47 there. Both constants were retuned
   upstream (`BLUESTEIN_SMOOTH_FACTOR` 1.9 → 2.1, `DEFAULT_BLUESTEIN_CUTOFF`
   47 → 29). See `calib_bluestein_x86.md`, `calib_cutoff_x86.md`,
   `calib_padchoice_x86.md`.
2. **The x86 penalty is a radix-5/7 penalty, not a SIMD penalty.** Per class,
   x86/aarch64 is 1.21× for pow2 but 1.60× for smooth composites, and within the
   smooth class it tracks factor content monotonically (2^15·3 → 1.21×;
   2^3·3^2·5^3·7^2 → 3.77×).
3. **A regression invisible on aarch64.** Wide `64×N dims=2` real transforms
   regressed 1.2–1.3×; localised to the `mapslices` → strided-views rewrite,
   whose copy had been an unlabelled copy-in optimisation for strided pencils.
   Fixed upstream; the rows are now 1.31–1.61× *faster* than baseline.
   See `branches/probe_*.md`.
4. **The codelet ceiling is lower on x86.** Up to n = 64 the codelet/FFTW ratios
   match aarch64 (0.51–1.54×), but n = 256 is 3.62× here against 1.73× there —
   AVX2's 16 vector registers against NEON's 32. A fixed codelet set should stop
   at 64. See `codelet_gen2_x86.md`, `codelet_rec_x86.md`.

## Measurement discipline

This is a laptop part under WSL2 with hybrid P/E cores that the guest cannot
identify or pin (all 22 vCPUs report `cpu_capacity` 1024, no `cpufreq`).
`thread_scaling_x86.md` quantifies it: a single-worker control re-measured inside
every trial round holds 0.2–1.7% spread, while multi-threaded runs spread 5–37%.

**Consequently, a cross-run delta below roughly 1.3× on this machine is not
evidence.** Six apparent regressions across these runs failed to reproduce when
re-measured back-to-back in a single session, including one that looked like a
30% regression and is in fact an 11% win. Every claim above rests on same-session
measurement; the scripts are included so the checks can be repeated.

## Files

| file | contents |
|:--|:--|
| `REPORT.md` | full rendered report, 540 cases (520 + the 20-row prime band) |
| `results_suite.json` | raw baseline results, band rows merged |
| `plots/` | SVG ratio and timing plots |
| `branches/*.json` | per-branch sweeps (integration pre- and post-fix, A, B, C, D, E) |
| `branches/probe_*.md` | back-to-back five-branch regression attribution |
| `calib_*_x86.md` | Bluestein pad factor, DFT/Bluestein crossover, pad choice |
| `codelet_*_x86.md` | §5.2 codelet tables reproduced on x86-64 |
| `thread_scaling_x86.md` | thread-placement probe with the 1-worker control |
| `check_B_small_x86.md` | back-to-back check of the small-`rfft` sizes |
| `*.jl` | every script used, so each measurement is reproducible |
| `PROVENANCE.md` | which commit each run measured |
