# FFTA.jl Performance Benchmarks

This directory contains a comprehensive benchmark suite to compare the performance of FFTA.jl against FFTW.jl.

## Structure

```
benchmark/
├── run_benchmarks.jl          # Main script to run all benchmarks
├── plot_results.jl            # Script to generate comparison plots
├── generate_html_report.jl    # Script to generate HTML report
├── Project.toml               # Dependencies for plotting
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
3. Generate comparison plots (combined and per-category)
4. Generate an interactive HTML report with embedded plots

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

Generate plots from existing results:
```bash
cd benchmark
julia --project=. plot_results.jl
```

Generate HTML report from existing results:
```bash
cd benchmark
julia --project=. generate_html_report.jl
```

## Output

The benchmark suite generates:

1. **results_ffta.json**: Raw benchmark data for FFTA.jl
2. **results_fftw.json**: Raw benchmark data for FFTW.jl
3. **benchmark_report.html**: Interactive HTML report with all plots and detailed tables
4. **performance_comparison_all.png**: Combined plot of Runtime/N vs N for all categories
5. **absolute_runtime_comparison.png**: Plot of absolute runtime vs N for all categories
6. **performance_*.png**: Individual plots for each category (odd powers of 2, even powers of 2, powers of 3, composite)

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

1. **Odd Powers of 2**: 2¹, 2³, 2⁵, 2⁷, 2⁹, 2¹¹, 2¹³, 2¹⁵
   - Sizes: 2, 8, 32, 128, 512, 2048, 8192, 32768
   - Tests radix-2 FFT with odd exponents

2. **Even Powers of 2**: 2², 2⁴, 2⁶, 2⁸, 2¹⁰, 2¹², 2¹⁴
   - Sizes: 4, 16, 64, 256, 1024, 4096, 16384
   - Tests radix-2 FFT with even exponents (often doubly-even cases)

3. **Powers of 3**: 3¹, 3², 3³, 3⁴, 3⁵, 3⁶, 3⁷, 3⁸, 3⁹
   - Sizes: 3, 9, 27, 81, 243, 729, 2187, 6561, 19683
   - Tests radix-3 FFT algorithms

4. **Composite Number**: 2×3×4×5×7 = 840
   - Tests mixed-radix FFT factorization

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
- Plots.jl (for visualization)

All dependencies are automatically installed when running the benchmarks.

## Continuous Integration

The benchmark suite integrates with GitHub Actions via the `.github/workflows/benchmarks.yml` workflow:

- **Automatic Runs**: Benchmarks run automatically on:
  - Pull requests that modify source code or benchmarks
  - Pushes to the main branch
  - Manual workflow dispatch

- **Artifacts**: Each CI run uploads:
  - HTML report (`benchmark_report.html`)
  - All plots (`.png` files)
  - Raw results (`.json` files)
  - Artifacts are retained for 30 days

- **PR Comments**: For pull requests, the workflow automatically posts a summary table with key results

To view benchmark results from a CI run:
1. Go to the Actions tab in the repository
2. Click on a "Benchmarks" workflow run
3. Download the `benchmark-results` artifact
4. Open `benchmark_report.html` in a browser
