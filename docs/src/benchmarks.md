# Performance Benchmarks

This page contains performance benchmarks comparing FFTA.jl against FFTW.jl.

## Interactive Benchmark Report

```@raw html
<div style="margin: 20px 0; padding: 15px; background-color: #e8f4f8; border-left: 4px solid #3498db; border-radius: 4px;">
    <p><strong>Note:</strong> The interactive benchmark report below is generated automatically by the CI pipeline when benchmarks are run.</p>
    <p>If you don't see the report, it means benchmarks haven't been run yet for this version of the documentation.</p>
</div>
```

```@raw html
<div id="benchmark-container">
    <iframe src="benchmarks/benchmark_report.html" style="width: 100%; height: 2000px; border: 1px solid #ddd; border-radius: 4px;" onload="this.style.height=(this.contentWindow.document.body.scrollHeight+50)+'px';">
        <p>Your browser does not support iframes. <a href="benchmarks/benchmark_report.html">Click here to view the benchmark report</a>.</p>
    </iframe>
</div>

<script>
// Check if benchmark report exists
fetch('benchmarks/benchmark_report.html')
    .then(response => {
        if (!response.ok) {
            document.getElementById('benchmark-container').innerHTML =
                '<div style="padding: 20px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; margin: 20px 0;">' +
                '<h3>⚠️ Benchmark Report Not Available</h3>' +
                '<p>The benchmark report has not been generated yet. Benchmarks are automatically run on:</p>' +
                '<ul>' +
                '<li>Pull requests that modify source code</li>' +
                '<li>Pushes to the main branch</li>' +
                '<li>Manual workflow dispatch</li>' +
                '</ul>' +
                '<p>To generate benchmarks locally, run:</p>' +
                '<pre><code>cd benchmark\njulia run_benchmarks.jl</code></pre>' +
                '</div>';
        }
    })
    .catch(error => {
        document.getElementById('benchmark-container').innerHTML =
            '<div style="padding: 20px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; margin: 20px 0;">' +
            '<h3>⚠️ Benchmark Report Not Available</h3>' +
            '<p>The benchmark report could not be loaded.</p>' +
            '</div>';
    });
</script>
```

## Running Benchmarks Locally

To run the benchmarks on your local machine:

```bash
cd benchmark
julia run_benchmarks.jl
```

This will:
1. Run FFTA benchmarks in an isolated environment
2. Run FFTW benchmarks in an isolated environment
3. Generate an interactive HTML report at `benchmark/benchmark_report.html`

For more details, see the [benchmark README](https://github.com/dannys4/FFTA.jl/tree/main/benchmark).

## Benchmark Methodology

The benchmark suite compares FFTA.jl (a pure Julia FFT implementation) against FFTW.jl (Julia bindings to the FFTW C library).

### Array Size Categories

Benchmarks are organized into categories based on array size structure:

- **Odd Powers of 2**: 2¹, 2³, 2⁵, ..., 2¹⁵ (2, 8, 32, 128, 512, 2048, 8192, 32768)
- **Even Powers of 2**: 2², 2⁴, 2⁶, ..., 2¹⁴ (4, 16, 64, 256, 1024, 4096, 16384)
- **Powers of 3**: 3¹, 3², 3³, ..., 3⁹ (3, 9, 27, 81, 243, 729, 2187, 6561, 19683)
- **Composite**: 3, 12, 60, 300, 2100, 23100 (cumulative products of 3, 4, 5, 5, 7, 11)
- **Prime Numbers**: 20 logarithmically-spaced primes up to 20,000

### Metrics

For each array size, we measure:
- **Median time**: Median execution time across 100 samples
- **Runtime/N**: Runtime divided by array length (shows scaling efficiency)
- **Mean/Min/Max time**: Statistical measures of performance

### Isolation

Each package is benchmarked in a completely separate Julia process to ensure:
- FFTW doesn't take precedence over FFTA when both are loaded
- Fair and accurate performance comparison
- No cross-contamination between implementations
