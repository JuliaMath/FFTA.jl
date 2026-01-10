#!/usr/bin/env julia

"""
Benchmark script for FFTW.jl
This script runs in an isolated environment with only FFTW.jl loaded.
"""

# Ensure packages are installed
import Pkg
Pkg.instantiate()

using FFTW
using BenchmarkTools
using JSON

# Load shared benchmark definitions
include(joinpath(@__DIR__, "..", "common_defs.jl"))

# Use shared definitions
const SIZE_CATEGORIES = create_size_categories()
const ALL_SIZES = sort(get_all_sizes())

# Number of samples for each benchmark
const SAMPLES = 100
const EVALS = 10

function benchmark_fftw()
    results = Dict{String, Any}()
    results["package"] = "FFTW"
    results["data"] = []
    results["categories"] = SIZE_CATEGORIES

    println("Benchmarking FFTW.jl...")
    println("=" ^ 50)

    for n in ALL_SIZES
        category = SIZE_CATEGORIES[n]
        println("Testing array size: $n (category: $category)")

        # Benchmark complex FFT
        x = randn(ComplexF64, n)
        trial = @benchmark fft($x) samples=SAMPLES evals=EVALS

        median_time = median(trial).time * 1e-9  # Convert to seconds
        runtime_per_element = median_time / n

        push!(results["data"], Dict(
            "size" => n,
            "category" => category,
            "median_time" => median_time,
            "runtime_per_element" => runtime_per_element,
            "mean_time" => mean(trial).time * 1e-9,
            "min_time" => minimum(trial).time * 1e-9,
            "max_time" => maximum(trial).time * 1e-9
        ))

        println("  Median time: $(median_time * 1e6) μs")
        println("  Time per element: $(runtime_per_element * 1e9) ns")
    end

    # Save results to JSON
    output_file = joinpath(@__DIR__, "..", "results_fftw.json")
    open(output_file, "w") do io
        JSON.print(io, results, 2)
    end

    println("\nResults saved to: $output_file")
    return results
end

# Run benchmark
benchmark_fftw()
