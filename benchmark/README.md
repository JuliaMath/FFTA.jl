# FFTA.jl Performance Benchmarks

This directory contains two benchmark suites comparing FFTA.jl against FFTW.jl:

1. **`suite.jl` + `report.jl`** — a comprehensive sweep (1D/2D/3D, batched
   `dims`, `Float32`/`Float64`, complex and real, plan vs. execution time,
   allocations, FFTW threading) that renders a markdown report with tables
   and SVG ratio plots. See [Comprehensive suite](#comprehensive-suite) below.
2. **`run_benchmarks.jl`** — the original 1D `fft`/`rfft` sweep that produces
   the interactive HTML report used by the documentation and CI.

## Comprehensive suite

```bash
cd benchmark
julia --project=. -t 8 suite.jl              # ~45 min; writes results_suite.json
julia --project=. report.jl results_suite.json   # writes REPORT.md + plots/*.svg
```

Options for `suite.jl`:

| flag | meaning |
|:--|:--|
| `--quick` | small sizes (≤ 2^16) and short time budget, for smoke testing |
| `--maxlog2 K` | largest 1D size 2^K (default 22) |
| `--seconds S` | time budget per measurement (default 0.5) |
| `--only 1d,nd,batched,threads` | run a subset of the sections |
| `--kinds fft,rfft` | run only complex (`fft`) or only real (`rfft`) cases |
| `--sizes 23,29,31` | restrict the 1D sweep to these sizes (to add rows to an existing run) |
| `--out FILE` | output JSON (default `results_suite.json`) |

`-t N` controls the thread count used for the FFTW-threading section (FFTA
is single-threaded). `report.jl` accepts several JSON files and merges them.

`report.jl` regenerates the whole `REPORT.md`; sections added by hand (cross
references to other runs) have to be re-added after re-rendering.

To compare two runs (e.g. before and after a change to FFTA):

```bash
julia --project=. compare.jl before.json after.json      # prints a markdown table
```

The suite loads FFTA and FFTW into the **same process** and reaches FFTA's
`plan_*` methods with `invoke`, because FFTW's methods on `StridedArray` are
more specific and would otherwise be selected. This also allows every FFTA
result to be checked against FFTW (the `rel. err` column of the report). The
timing loop is a small custom one (minimum over samples within a time budget)
rather than `@benchmark`, because BenchmarkTools' per-call-site compilation
dominates a sweep of several hundred cases.

What is measured per case:

* **exec** — planned execution (`mul!(y, p, x)`; FFTA real plans only support
  `p * x`, which includes output allocation)
* **plan** — plan construction (FFTW with `FFTW.ESTIMATE`; a `FFTW.MEASURE`
  column is included for 1D power-of-two `ComplexF64`)
* **cold** — plan + execute, i.e. the one-shot `fft(x)` path
* **alloc** — bytes allocated by one planned execution
* **rel. err** — ‖y_FFTA − y_FFTW‖ / ‖y_FFTW‖

The committed [`REPORT.md`](REPORT.md) was produced on the machine described
in its *Environment* section; re-run the suite to obtain numbers for your own
hardware.

## Original 1D suite (`run_benchmarks.jl`)

## Structure

```
benchmark/
├── suite.jl                   # Comprehensive FFTA-vs-FFTW sweep (see above)
├── report.jl                  # Renders REPORT.md + plots/ from suite.jl output
├── REPORT.md, plots/          # Committed results of the comprehensive suite
├── run_benchmarks.jl          # Main script to run the original 1D benchmarks
├── generate_html_report.jl    # Script to generate HTML report with Plotly.js
├── Project.toml               # Dependencies (FFTW, BenchmarkTools, JSON, ...)
├── ffta_env/                  # Isolated environment for FFTA
│   ├── bench_ffta.jl         # FFTA benchmark script
│   └── Project.toml          # FFTA dependencies
├── fftw_env/                  # Isolated environment for FFTW
│   ├── bench_fftw.jl         # FFTW benchmark script
│   └── Project.toml          # FFTW dependencies
└── README.md                  # This file
```

## Why Separate Environments?

FFTW.jl takes precedence over other FFT implementations when loaded in the same environment. To ensure fair and accurate benchmarks, we run FFTA and FFTW benchmarks in completely separate Julia processes with their own isolated environments.

## Usage

### Quick Start

Run all benchmarks and generate plots:

```bash
cd benchmark
julia run_benchmarks.jl
```

This will:
1. Run FFTA benchmarks in an isolated environment
2. Run FFTW benchmarks in an isolated environment
3. Generate an interactive HTML report with embedded Plotly.js charts

### Individual Benchmarks

Run FFTA benchmark only:
```bash
cd benchmark/ffta_env
julia --project=. bench_ffta.jl
```

Run FFTW benchmark only:
```bash
cd benchmark/fftw_env
julia --project=. bench_fftw.jl
```

Generate HTML report from existing JSON results:
```bash
cd benchmark
julia --project=. generate_html_report.jl
```

### Building Documentation with Benchmarks

To build the documentation with benchmark results included:

```bash
# 1. Run benchmarks
cd benchmark
julia run_benchmarks.jl
cd ..

# 2. Build documentation
julia --project=docs docs/make.jl
```

The `docs/make.jl` script will automatically detect and copy benchmark results from `benchmark/` to `docs/src/assets/benchmarks/` before building. The documentation will include the interactive benchmark report.

## Output

The benchmark suite generates:

1. **results_ffta.json**: Raw benchmark data for FFTA.jl (complex FFT)
2. **results_fftw.json**: Raw benchmark data for FFTW.jl (complex FFT)
3. **results_ffta_rfft.json**: Raw benchmark data for FFTA.jl (real FFT)
4. **results_fftw_rfft.json**: Raw benchmark data for FFTW.jl (real FFT)
5. **benchmark_report.html**: Self-contained interactive HTML report with:
   - Embedded JSON data
   - Client-side Plotly.js charts (no external files needed)
   - Combined Runtime/N vs N plot for all categories (both complex and real FFT)
   - Absolute runtime plot for all categories (both complex and real FFT)
   - Individual plots for each category (odd/even powers of 2, powers of 3, composite, primes)
   - Detailed results tables with speedup comparisons
   - Separate sections for complex FFT and real FFT results

## Metrics

For each array size, we measure:
- **median_time**: Median execution time
- **runtime_per_element**: Runtime divided by array length (shows scaling efficiency)
- **mean_time**: Mean execution time
- **min_time**: Minimum execution time
- **max_time**: Maximum execution time

## Array Sizes Tested

The benchmarks test various array sizes categorized by their mathematical structure to understand FFT performance characteristics:

### Categories

1. **Odd Powers of 2**: 2¹, 2³, 2⁵, 2⁷, 2⁹, 2¹¹, 2¹³, 2¹⁵, 2¹⁷, 2¹⁹
   - Sizes: 2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288
   - Tests radix-2 FFT with odd exponents

2. **Even Powers of 2**: 2², 2⁴, 2⁶, 2⁸, 2¹⁰, 2¹², 2¹⁴, 2¹⁶, 2¹⁸, 2²⁰
   - Sizes: 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576
   - Tests radix-2 FFT with even exponents (often doubly-even cases)

3. **Powers of 3**: 3¹, 3², 3³, 3⁴, 3⁵, 3⁶, 3⁷, 3⁸, 3⁹
   - Sizes: 3, 9, 27, 81, 243, 729, 2187, 6561, 19683
   - Tests radix-3 FFT algorithms

4. **Composite**: 3, 12, 60, 300, 2100, 23100
   - Cumulative products of 3, 4, 5, 5, 7, 11
   - Tests mixed-radix FFT factorization with increasing complexity

5. **Prime Numbers**: 20 logarithmically-spaced primes up to 20,000
   - Tests FFT performance on prime-sized arrays with logarithmic spacing
   - Prime sizes require specialized FFT algorithms (e.g., Bluestein's algorithm)
   - Logarithmic spacing ensures coverage from small to large primes

All tests are run for both:
- **Complex FFT**: Complex double-precision input arrays (`ComplexF64`)
- **Real FFT**: Real double-precision input arrays (`Float64`)

The real FFT (rfft) is optimized for real-valued input and exploits conjugate symmetry, typically achieving ~2x speedup over complex FFT for real data.

## Interpreting Results

### Runtime/N vs N Plot

This plot shows how efficiently each implementation scales with array size:
- **Ideal FFT**: Should show O(log N) growth (since FFT is O(N log N), Runtime/N is O(log N))
- **Flat line**: Indicates optimal scaling
- **Upward trend**: Indicates scaling overhead (cache effects, algorithm inefficiencies)

### Absolute Runtime Plot

This plot shows the raw execution time for each array size:
- Lower is better
- Should show approximately O(N log N) growth
- Useful for comparing absolute performance at specific sizes

## Dependencies

The benchmark suite requires:
- Julia 1.x
- FFTA.jl (the package being benchmarked)
- FFTW.jl (for comparison)
- BenchmarkTools.jl (for accurate timing)
- JSON.jl (for storing results)
- Primes.jl (for generating prime-sized arrays)
- Dates.jl (standard library, for timestamps)

All dependencies are automatically installed when running the benchmarks.

**Note:** Plots are generated client-side using Plotly.js (loaded from CDN in the HTML). No Julia plotting packages are required.

## Continuous Integration

The benchmark suite integrates with GitHub Actions via the `.github/workflows/benchmarks.yml` workflow:

- **Automatic Runs**: Benchmarks run automatically on:
  - Pull requests that modify source code or benchmarks
  - Pushes to the main branch
  - Manual workflow dispatch

- **Artifacts**: Each CI run uploads:
  - Interactive HTML report (`benchmark_report.html`) with embedded Plotly.js charts
  - Raw benchmark results (`.json` files)
  - Benchmark logs for debugging
  - Artifacts are retained for 30 days

To view benchmark results from a CI run:
1. Go to the Actions tab in the repository
2. Click on a "Benchmarks" workflow run
3. Download the `benchmark-results` artifact
4. Open `benchmark_report.html` in a browser
