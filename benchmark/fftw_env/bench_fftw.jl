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

# Load shared benchmark definitions from parent project
# (common_defs.jl uses Primes which is in the parent benchmark project)
const SIZE_CATEGORIES, ALL_SIZES = let
    parent_project = joinpath(@__DIR__, "..")
    Pkg.activate(parent_project)
    Pkg.instantiate()
    include(joinpath(parent_project, "common_defs.jl"))
    # Restore this project
    Pkg.activate(@__DIR__)
    (create_size_categories(), sort(get_all_sizes()))
end

# Number of samples for each benchmark
const SAMPLES = 100
const EVALS = 10

function benchmark_fftw_complex()
    results = Dict{String, Any}()
    results["package"] = "FFTW"
    results["fft_type"] = "complex"
    results["data"] = []
    results["categories"] = SIZE_CATEGORIES

    println("Benchmarking FFTW.jl (Complex FFT)...")
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

function benchmark_fftw_real()
    results = Dict{String, Any}()
    results["package"] = "FFTW"
    results["fft_type"] = "real"
    results["data"] = []
    results["categories"] = SIZE_CATEGORIES

    println("\nBenchmarking FFTW.jl (Real FFT)...")
    println("=" ^ 50)

    for n in ALL_SIZES
        category = SIZE_CATEGORIES[n]
        println("Testing array size: $n (category: $category)")

        # Benchmark real FFT
        x = randn(Float64, n)
        trial = @benchmark rfft($x) samples=SAMPLES evals=EVALS

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
    output_file = joinpath(@__DIR__, "..", "results_fftw_rfft.json")
    open(output_file, "w") do io
        JSON.print(io, results, 2)
    end

    println("\nResults saved to: $output_file")
    return results
end

# Run benchmarks
benchmark_fftw_complex()
benchmark_fftw_real()
