#!/usr/bin/env julia

"""
Main script to run FFTA and FFTW benchmarks in separate environments
This ensures that FFTW doesn't take precedence over FFTA
"""

println("=" ^ 70)
println("FFT Performance Benchmark Suite")
println("Comparing FFTA.jl vs FFTW.jl")
println("=" ^ 70)
println()

# Get the benchmark directory
benchmark_dir = @__DIR__

# Run FFTA benchmark in separate process
println("Step 1/4: Running FFTA.jl benchmark...")
println("-" ^ 70)
ffta_script = joinpath(benchmark_dir, "ffta_env", "bench_ffta.jl")
ffta_project = joinpath(benchmark_dir, "ffta_env")

ffta_cmd = `julia --project=$ffta_project $ffta_script`

result = run(ffta_cmd)
if !success(result)
    error("FFTA benchmark failed with exit code ", result.exitcode)
end
println()

# Run FFTW benchmark in separate process
println("Step 2/4: Running FFTW.jl benchmark...")
println("-" ^ 70)
fftw_script = joinpath(benchmark_dir, "fftw_env", "bench_fftw.jl")
fftw_project = joinpath(benchmark_dir, "fftw_env")

fftw_cmd = `julia --project=$fftw_project $fftw_script`

result = run(fftw_cmd)
if !success(result)
    error("FFTW benchmark failed with exit code ", result.exitcode)
end
println()

# Generate plots
println("Step 3/4: Generating comparison plots...")
println("-" ^ 70)
plot_script = joinpath(benchmark_dir, "plot_results.jl")
plot_project = benchmark_dir

plot_cmd = `julia --project=$plot_project $plot_script`

result = run(plot_cmd)
if !success(result)
    error("Plotting failed with exit code ", result.exitcode)
end
println()

# Generate HTML report
println("Step 4/4: Generating HTML report...")
println("-" ^ 70)
html_script = joinpath(benchmark_dir, "generate_html_report.jl")

html_cmd = `julia --project=$plot_project $html_script`

result = run(html_cmd)
if !success(result)
    error("HTML generation failed with exit code ", result.exitcode)
end
println()

println("=" ^ 70)
println("Benchmark suite completed successfully!")
println("=" ^ 70)
println()
println("Results:")
println("  - FFTA results: ", joinpath(benchmark_dir, "results_ffta.json"))
println("  - FFTW results: ", joinpath(benchmark_dir, "results_fftw.json"))
println("  - HTML report: ", joinpath(benchmark_dir, "benchmark_report.html"))
println("  - Combined plot: ", joinpath(benchmark_dir, "performance_comparison_all.png"))
println("  - Absolute runtime plot: ", joinpath(benchmark_dir, "absolute_runtime_comparison.png"))
println("  - Category plots: performance_*.png")
println()
