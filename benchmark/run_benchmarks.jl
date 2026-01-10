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
println("Step 1/3: Running FFTA.jl benchmark...")
println("-" ^ 70)
ffta_script = joinpath(benchmark_dir, "ffta_env", "bench_ffta.jl")
ffta_project = joinpath(benchmark_dir, "ffta_env")

ffta_cmd = `julia --project=$ffta_project -e "
    import Pkg
    Pkg.instantiate()
    include(\"$ffta_script\")
"`

success(run(ffta_cmd)) || error("FFTA benchmark failed!")
println()

# Run FFTW benchmark in separate process
println("Step 2/3: Running FFTW.jl benchmark...")
println("-" ^ 70)
fftw_script = joinpath(benchmark_dir, "fftw_env", "bench_fftw.jl")
fftw_project = joinpath(benchmark_dir, "fftw_env")

fftw_cmd = `julia --project=$fftw_project -e "
    import Pkg
    Pkg.instantiate()
    include(\"$fftw_script\")
"`

success(run(fftw_cmd)) || error("FFTW benchmark failed!")
println()

# Generate plots
println("Step 3/3: Generating comparison plots...")
println("-" ^ 70)
plot_script = joinpath(benchmark_dir, "plot_results.jl")
plot_project = benchmark_dir

plot_cmd = `julia --project=$plot_project -e "
    import Pkg
    Pkg.instantiate()
    include(\"$plot_script\")
"`

success(run(plot_cmd)) || error("Plotting failed!")
println()

println("=" ^ 70)
println("Benchmark suite completed successfully!")
println("=" ^ 70)
println()
println("Results:")
println("  - FFTA results: $(joinpath(benchmark_dir, \"results_ffta.json\"))")
println("  - FFTW results: $(joinpath(benchmark_dir, \"results_fftw.json\"))")
println("  - Performance plot: $(joinpath(benchmark_dir, \"performance_comparison.png\"))")
println("  - Absolute runtime plot: $(joinpath(benchmark_dir, \"absolute_runtime_comparison.png\"))")
println()
