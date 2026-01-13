# FFTA.jl Performance Benchmarks

This directory contains a comprehensive benchmark suite to compare the performance of FFTA.jl against FFTW.jl.

## Structure

```
benchmark/
├── run_benchmarks.jl          # Main script to run all benchmarks
├── generate_html_report.jl    # Script to generate HTML report with Plotly.js
├── Project.toml               # Dependencies (JSON, Dates)
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

1. **results_ffta.json**: Raw benchmark data for FFTA.jl
2. **results_fftw.json**: Raw benchmark data for FFTW.jl
3. **benchmark_report.html**: Self-contained interactive HTML report with:
   - Embedded JSON data
   - Client-side Plotly.js charts (no external files needed)
   - Combined Runtime/N vs N plot for all categories
   - Absolute runtime plot for all categories
   - Individual plots for each category (odd/even powers of 2, powers of 3, composite, primes)
   - Detailed results tables with speedup comparisons

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

All tests use complex double-precision arrays (`ComplexF64`)

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
