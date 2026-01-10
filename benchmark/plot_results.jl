#!/usr/bin/env julia

"""
Plot benchmark results for FFTA.jl vs FFTW.jl
Plots runtime/N vs N to show scaling behavior
"""

using JSON
using Plots

function plot_benchmark_results()
    # Load results
    ffta_file = joinpath(@__DIR__, "results_ffta.json")
    fftw_file = joinpath(@__DIR__, "results_fftw.json")

    if !isfile(ffta_file) || !isfile(fftw_file)
        error("Benchmark results not found. Run benchmarks first!")
    end

    ffta_results = JSON.parsefile(ffta_file)
    fftw_results = JSON.parsefile(fftw_file)

    # Extract data
    ffta_data = ffta_results["data"]
    fftw_data = fftw_results["data"]

    ffta_sizes = [d["size"] for d in ffta_data]
    ffta_runtime_per_n = [d["runtime_per_element"] for d in ffta_data]

    fftw_sizes = [d["size"] for d in fftw_data]
    fftw_runtime_per_n = [d["runtime_per_element"] for d in fftw_data]

    # Create plot: Runtime/N vs N
    plt = plot(
        ffta_sizes,
        ffta_runtime_per_n .* 1e9,  # Convert to nanoseconds
        label="FFTA.jl",
        marker=:circle,
        markersize=6,
        linewidth=2,
        xscale=:log10,
        yscale=:log10,
        xlabel="Array Length (N)",
        ylabel="Runtime / N (nanoseconds)",
        title="FFT Performance Comparison: Runtime/N vs N",
        legend=:best,
        grid=true,
        size=(800, 600),
        dpi=150
    )

    plot!(
        plt,
        fftw_sizes,
        fftw_runtime_per_n .* 1e9,
        label="FFTW.jl",
        marker=:square,
        markersize=6,
        linewidth=2
    )

    # Add theoretical O(log N) reference line for comparison
    # FFT complexity is O(N log N), so runtime/N should be O(log N)
    ref_sizes = ffta_sizes
    ref_baseline = minimum(fftw_runtime_per_n .* 1e9)
    ref_line = ref_baseline .* log2.(ref_sizes) ./ log2(ref_sizes[1])
    plot!(
        plt,
        ref_sizes,
        ref_line,
        label="O(log N) reference",
        linestyle=:dash,
        linewidth=1.5,
        color=:gray,
        alpha=0.7
    )

    # Save plot
    output_file = joinpath(@__DIR__, "performance_comparison.png")
    savefig(plt, output_file)
    println("Plot saved to: $output_file")

    # Also create a plot of absolute runtime
    plt2 = plot(
        ffta_sizes,
        [d["median_time"] * 1e6 for d in ffta_data],
        label="FFTA.jl",
        marker=:circle,
        markersize=6,
        linewidth=2,
        xscale=:log10,
        yscale=:log10,
        xlabel="Array Length (N)",
        ylabel="Runtime (microseconds)",
        title="FFT Performance Comparison: Absolute Runtime",
        legend=:best,
        grid=true,
        size=(800, 600),
        dpi=150
    )

    plot!(
        plt2,
        fftw_sizes,
        [d["median_time"] * 1e6 for d in fftw_data],
        label="FFTW.jl",
        marker=:square,
        markersize=6,
        linewidth=2
    )

    output_file2 = joinpath(@__DIR__, "absolute_runtime_comparison.png")
    savefig(plt2, output_file2)
    println("Plot saved to: $output_file2")

    # Print summary statistics
    println("\n" * "=" ^ 60)
    println("Performance Summary")
    println("=" ^ 60)

    for (i, n) in enumerate(ffta_sizes)
        ffta_time = ffta_data[i]["median_time"] * 1e6
        fftw_time = fftw_data[i]["median_time"] * 1e6
        speedup = fftw_time / ffta_time

        println("N = $n:")
        println("  FFTA: $(round(ffta_time, digits=3)) μs")
        println("  FFTW: $(round(fftw_time, digits=3)) μs")
        if speedup > 1
            println("  FFTA is $(round(speedup, digits=2))x faster")
        else
            println("  FFTW is $(round(1/speedup, digits=2))x faster")
        end
        println()
    end

    return plt, plt2
end

# Run plotting
plot_benchmark_results()
