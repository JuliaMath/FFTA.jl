#!/usr/bin/env julia

"""
Plot benchmark results for FFTA.jl vs FFTW.jl
Plots runtime/N vs N to show scaling behavior, with separate plots for each category
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

    # Category information
    categories = ["odd_power_of_2", "even_power_of_2", "power_of_3", "composite"]
    category_labels = Dict(
        "odd_power_of_2" => "Odd Powers of 2 (2¹, 2³, 2⁵, ...)",
        "even_power_of_2" => "Even Powers of 2 (2², 2⁴, 2⁶, ...)",
        "power_of_3" => "Powers of 3 (3¹, 3², 3³, ...)",
        "composite" => "Composite (2×3×4×5×7)"
    )

    # Category colors
    category_colors = Dict(
        "odd_power_of_2" => :blue,
        "even_power_of_2" => :red,
        "power_of_3" => :green,
        "composite" => :purple
    )

    # Create combined plot: Runtime/N vs N (all categories)
    plt_combined = plot(
        xlabel="Array Length (N)",
        ylabel="Runtime / N (nanoseconds)",
        title="FFT Performance: Runtime/N vs N (All Categories)",
        xscale=:log10,
        yscale=:log10,
        legend=:outertopright,
        grid=true,
        size=(1200, 700),
        dpi=150
    )

    # Plot each category
    for category in categories
        # Filter data for this category
        ffta_cat = [d for d in ffta_data if get(d, "category", "") == category]
        fftw_cat = [d for d in fftw_data if get(d, "category", "") == category]

        if !isempty(ffta_cat)
            ffta_sizes = [d["size"] for d in ffta_cat]
            ffta_runtime_per_n = [d["runtime_per_element"] for d in ffta_cat]

            plot!(
                plt_combined,
                ffta_sizes,
                ffta_runtime_per_n .* 1e9,
                label="FFTA: $(category_labels[category])",
                marker=:circle,
                markersize=5,
                linewidth=2,
                color=category_colors[category],
                linestyle=:solid
            )
        end

        if !isempty(fftw_cat)
            fftw_sizes = [d["size"] for d in fftw_cat]
            fftw_runtime_per_n = [d["runtime_per_element"] for d in fftw_cat]

            plot!(
                plt_combined,
                fftw_sizes,
                fftw_runtime_per_n .* 1e9,
                label="FFTW: $(category_labels[category])",
                marker=:square,
                markersize=5,
                linewidth=2,
                color=category_colors[category],
                linestyle=:dash
            )
        end
    end

    # Save combined plot
    output_file = joinpath(@__DIR__, "performance_comparison_all.png")
    savefig(plt_combined, output_file)
    println("Combined plot saved to: $output_file")

    # Create separate plots for each category
    category_plots = []
    for category in categories
        ffta_cat = [d for d in ffta_data if get(d, "category", "") == category]
        fftw_cat = [d for d in fftw_data if get(d, "category", "") == category]

        if isempty(ffta_cat) && isempty(fftw_cat)
            continue
        end

        plt_cat = plot(
            xlabel="Array Length (N)",
            ylabel="Runtime / N (nanoseconds)",
            title="$(category_labels[category])",
            xscale=:log10,
            yscale=:log10,
            legend=:best,
            grid=true,
            size=(600, 500),
            dpi=150
        )

        if !isempty(ffta_cat)
            ffta_sizes = [d["size"] for d in ffta_cat]
            ffta_runtime_per_n = [d["runtime_per_element"] for d in ffta_cat]

            plot!(
                plt_cat,
                ffta_sizes,
                ffta_runtime_per_n .* 1e9,
                label="FFTA.jl",
                marker=:circle,
                markersize=6,
                linewidth=2,
                color=category_colors[category]
            )
        end

        if !isempty(fftw_cat)
            fftw_sizes = [d["size"] for d in fftw_cat]
            fftw_runtime_per_n = [d["runtime_per_element"] for d in fftw_cat]

            plot!(
                plt_cat,
                fftw_sizes,
                fftw_runtime_per_n .* 1e9,
                label="FFTW.jl",
                marker=:square,
                markersize=6,
                linewidth=2,
                color=:gray
            )
        end

        push!(category_plots, (category, plt_cat))

        # Save individual category plot
        cat_output = joinpath(@__DIR__, "performance_$(category).png")
        savefig(plt_cat, cat_output)
        println("Category plot saved to: $cat_output")
    end

    # Create absolute runtime comparison
    plt_absolute = plot(
        xlabel="Array Length (N)",
        ylabel="Runtime (microseconds)",
        title="FFT Absolute Runtime Comparison",
        xscale=:log10,
        yscale=:log10,
        legend=:outertopright,
        grid=true,
        size=(1200, 700),
        dpi=150
    )

    for category in categories
        ffta_cat = [d for d in ffta_data if get(d, "category", "") == category]
        fftw_cat = [d for d in fftw_data if get(d, "category", "") == category]

        if !isempty(ffta_cat)
            ffta_sizes = [d["size"] for d in ffta_cat]
            ffta_times = [d["median_time"] * 1e6 for d in ffta_cat]

            plot!(
                plt_absolute,
                ffta_sizes,
                ffta_times,
                label="FFTA: $(category_labels[category])",
                marker=:circle,
                markersize=5,
                linewidth=2,
                color=category_colors[category],
                linestyle=:solid
            )
        end

        if !isempty(fftw_cat)
            fftw_sizes = [d["size"] for d in fftw_cat]
            fftw_times = [d["median_time"] * 1e6 for d in fftw_cat]

            plot!(
                plt_absolute,
                fftw_sizes,
                fftw_times,
                label="FFTW: $(category_labels[category])",
                marker=:square,
                markersize=5,
                linewidth=2,
                color=category_colors[category],
                linestyle=:dash
            )
        end
    end

    output_file2 = joinpath(@__DIR__, "absolute_runtime_comparison.png")
    savefig(plt_absolute, output_file2)
    println("Absolute runtime plot saved to: $output_file2")

    # Print summary statistics
    println("\n" * "=" ^ 70)
    println("Performance Summary by Category")
    println("=" ^ 70)

    for category in categories
        ffta_cat = [d for d in ffta_data if get(d, "category", "") == category]
        fftw_cat = [d for d in fftw_data if get(d, "category", "") == category]

        if isempty(ffta_cat) && isempty(fftw_cat)
            continue
        end

        println("\n$(category_labels[category]):")
        println("-" ^ 70)

        all_sizes = sort(unique(vcat([d["size"] for d in ffta_cat], [d["size"] for d in fftw_cat])))

        for n in all_sizes
            ffta_entry = findfirst(d -> d["size"] == n, ffta_cat)
            fftw_entry = findfirst(d -> d["size"] == n, fftw_cat)

            if !isnothing(ffta_entry) && !isnothing(fftw_entry)
                ffta_time = ffta_cat[ffta_entry]["median_time"] * 1e6
                fftw_time = fftw_cat[fftw_entry]["median_time"] * 1e6
                speedup = fftw_time / ffta_time

                println("  N = $n:")
                println("    FFTA: $(round(ffta_time, digits=3)) μs")
                println("    FFTW: $(round(fftw_time, digits=3)) μs")
                if speedup > 1
                    println("    FFTA is $(round(speedup, digits=2))x faster")
                else
                    println("    FFTW is $(round(1/speedup, digits=2))x faster")
                end
            end
        end
    end
    println("\n" * "=" ^ 70)

    return plt_combined, plt_absolute, category_plots
end

# Run plotting
plot_benchmark_results()
