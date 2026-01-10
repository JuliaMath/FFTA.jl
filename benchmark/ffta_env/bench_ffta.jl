#!/usr/bin/env julia

"""
Benchmark script for FFTA.jl
This script runs in an isolated environment with only FFTA.jl loaded.
"""

# Add FFTA from the parent directory
import Pkg
ffta_path = joinpath(@__DIR__, "..", "..")
Pkg.develop(path=ffta_path)
Pkg.instantiate()

using FFTA
using BenchmarkTools
using JSON

# Array sizes categorized by structure
const ODD_POWERS_OF_2 = [2^p for p in 1:2:15]  # 2, 8, 32, 128, 512, 2048, 8192, 32768
const EVEN_POWERS_OF_2 = [2^p for p in 2:2:14]  # 4, 16, 64, 256, 1024, 4096, 16384
const POWERS_OF_3 = [3^p for p in 1:9]  # 3, 9, 27, 81, 243, 729, 2187, 6561, 19683
const COMPOSITE = [2*3*4*5*7]  # 840

# Combine all sizes and categorize them
const SIZE_CATEGORIES = Dict{Int, String}()
for n in ODD_POWERS_OF_2
    SIZE_CATEGORIES[n] = "odd_power_of_2"
end
for n in EVEN_POWERS_OF_2
    SIZE_CATEGORIES[n] = "even_power_of_2"
end
for n in POWERS_OF_3
    SIZE_CATEGORIES[n] = "power_of_3"
end
for n in COMPOSITE
    SIZE_CATEGORIES[n] = "composite"
end

const ALL_SIZES = sort(collect(keys(SIZE_CATEGORIES)))

# Number of samples for each benchmark
const SAMPLES = 100
const EVALS = 10

function benchmark_ffta()
    results = Dict{String, Any}()
    results["package"] = "FFTA"
    results["data"] = []
    results["categories"] = SIZE_CATEGORIES

    println("Benchmarking FFTA.jl...")
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
    output_file = joinpath(@__DIR__, "..", "results_ffta.json")
    open(output_file, "w") do io
        JSON.print(io, results, 2)
    end

    println("\nResults saved to: $output_file")
    return results
end

# Run benchmark
benchmark_ffta()
