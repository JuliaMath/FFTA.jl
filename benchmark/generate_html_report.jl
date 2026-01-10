#!/usr/bin/env julia

"""
Generate HTML report with embedded plots for FFTA vs FFTW benchmarks
"""

using JSON
using Base64

function generate_html_report()
    # Load results
    ffta_file = joinpath(@__DIR__, "results_ffta.json")
    fftw_file = joinpath(@__DIR__, "results_fftw.json")

    if !isfile(ffta_file) || !isfile(fftw_file)
        error("Benchmark results not found. Run benchmarks first!")
    end

    ffta_results = JSON.parsefile(ffta_file)
    fftw_results = JSON.parsefile(fftw_file)

    # Find all plot images
    plot_dir = @__DIR__
    plot_files = [
        ("performance_comparison_all.png", "Combined Performance: Runtime/N vs N"),
        ("absolute_runtime_comparison.png", "Absolute Runtime Comparison"),
        ("performance_odd_power_of_2.png", "Odd Powers of 2"),
        ("performance_even_power_of_2.png", "Even Powers of 2"),
        ("performance_power_of_3.png", "Powers of 3"),
        ("performance_composite.png", "Composite Number (840)")
    ]

    # Start HTML document
    html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>FFTA.jl vs FFTW.jl Performance Benchmark</title>
        <style>
            body {
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
                max-width: 1400px;
                margin: 0 auto;
                padding: 20px;
                background-color: #f5f5f5;
                color: #333;
            }
            h1 {
                color: #2c3e50;
                border-bottom: 3px solid #3498db;
                padding-bottom: 10px;
            }
            h2 {
                color: #34495e;
                margin-top: 30px;
                border-bottom: 2px solid #95a5a6;
                padding-bottom: 5px;
            }
            h3 {
                color: #7f8c8d;
            }
            .summary {
                background-color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin: 20px 0;
            }
            .plot-container {
                background-color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin: 20px 0;
                text-align: center;
            }
            .plot-container img {
                max-width: 100%;
                height: auto;
                border-radius: 4px;
            }
            .category-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(600px, 1fr));
                gap: 20px;
                margin: 20px 0;
            }
            table {
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                background-color: white;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            th, td {
                padding: 12px;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }
            th {
                background-color: #3498db;
                color: white;
                font-weight: 600;
            }
            tr:hover {
                background-color: #f5f5f5;
            }
            .faster {
                color: #27ae60;
                font-weight: 600;
            }
            .slower {
                color: #e74c3c;
                font-weight: 600;
            }
            .metadata {
                background-color: #ecf0f1;
                padding: 15px;
                border-radius: 4px;
                margin: 10px 0;
                font-family: monospace;
                font-size: 0.9em;
            }
            .info-box {
                background-color: #e8f4f8;
                border-left: 4px solid #3498db;
                padding: 15px;
                margin: 20px 0;
                border-radius: 4px;
            }
        </style>
    </head>
    <body>
        <h1>FFTA.jl vs FFTW.jl Performance Benchmark Report</h1>

        <div class="info-box">
            <strong>Note:</strong> This benchmark compares FFTA.jl (a pure Julia FFT implementation) against
            FFTW.jl (Julia bindings to the FFTW C library). Each package was benchmarked in a separate
            Julia process to ensure isolation and prevent interference.
        </div>

        <div class="summary">
            <h2>Benchmark Configuration</h2>
            <div class="metadata">
                <strong>Benchmarking Tool:</strong> BenchmarkTools.jl<br>
                <strong>Samples per size:</strong> 100<br>
                <strong>Evaluations per sample:</strong> 10<br>
                <strong>Data type:</strong> ComplexF64<br>
                <strong>Total array sizes tested:</strong> $(length(ffta_results["data"]))<br>
            </div>

            <h3>Array Size Categories</h3>
            <ul>
                <li><strong>Odd Powers of 2:</strong> 2¹, 2³, 2⁵, 2⁷, 2⁹, 2¹¹, 2¹³, 2¹⁵ (2, 8, 32, 128, 512, 2048, 8192, 32768)</li>
                <li><strong>Even Powers of 2:</strong> 2², 2⁴, 2⁶, 2⁸, 2¹⁰, 2¹², 2¹⁴ (4, 16, 64, 256, 1024, 4096, 16384)</li>
                <li><strong>Powers of 3:</strong> 3¹, 3², 3³, 3⁴, 3⁵, 3⁶, 3⁷, 3⁸, 3⁹ (3, 9, 27, 81, 243, 729, 2187, 6561, 19683)</li>
                <li><strong>Composite:</strong> 2×3×4×5×7 = 840</li>
            </ul>
        </div>
    """

    # Add plots
    html *= "\n<h2>Performance Visualizations</h2>\n"

    for (filename, title) in plot_files
        plot_path = joinpath(plot_dir, filename)
        if isfile(plot_path)
            # Read and encode image as base64
            img_data = read(plot_path)
            img_base64 = base64encode(img_data)

            html *= """
            <div class="plot-container">
                <h3>$title</h3>
                <img src="data:image/png;base64,$img_base64" alt="$title">
            </div>
            """
        end
    end

    # Add detailed results table
    html *= """
    <h2>Detailed Results</h2>
    """

    # Category information
    categories = ["odd_power_of_2", "even_power_of_2", "power_of_3", "composite"]
    category_labels = Dict(
        "odd_power_of_2" => "Odd Powers of 2",
        "even_power_of_2" => "Even Powers of 2",
        "power_of_3" => "Powers of 3",
        "composite" => "Composite (840)"
    )

    ffta_data = ffta_results["data"]
    fftw_data = fftw_results["data"]

    for category in categories
        ffta_cat = [d for d in ffta_data if get(d, "category", "") == category]
        fftw_cat = [d for d in fftw_data if get(d, "category", "") == category]

        if isempty(ffta_cat) && isempty(fftw_cat)
            continue
        end

        html *= """
        <h3>$(category_labels[category])</h3>
        <table>
            <thead>
                <tr>
                    <th>Array Size (N)</th>
                    <th>FFTA Time (μs)</th>
                    <th>FFTW Time (μs)</th>
                    <th>FFTA Runtime/N (ns)</th>
                    <th>FFTW Runtime/N (ns)</th>
                    <th>Speedup</th>
                </tr>
            </thead>
            <tbody>
        """

        all_sizes = sort(unique(vcat([d["size"] for d in ffta_cat], [d["size"] for d in fftw_cat])))

        for n in all_sizes
            ffta_entry = findfirst(d -> d["size"] == n, ffta_cat)
            fftw_entry = findfirst(d -> d["size"] == n, fftw_cat)

            if !isnothing(ffta_entry) && !isnothing(fftw_entry)
                ffta_time = ffta_cat[ffta_entry]["median_time"] * 1e6
                fftw_time = fftw_cat[fftw_entry]["median_time"] * 1e6
                ffta_per_n = ffta_cat[ffta_entry]["runtime_per_element"] * 1e9
                fftw_per_n = fftw_cat[fftw_entry]["runtime_per_element"] * 1e9
                speedup = fftw_time / ffta_time

                speedup_class = speedup > 1 ? "faster" : "slower"
                speedup_text = speedup > 1 ? "$(round(speedup, digits=2))x (FFTA faster)" : "$(round(1/speedup, digits=2))x (FFTW faster)"

                html *= """
                <tr>
                    <td>$n</td>
                    <td>$(round(ffta_time, digits=3))</td>
                    <td>$(round(fftw_time, digits=3))</td>
                    <td>$(round(ffta_per_n, digits=3))</td>
                    <td>$(round(fftw_per_n, digits=3))</td>
                    <td class="$speedup_class">$speedup_text</td>
                </tr>
                """
            end
        end

        html *= """
            </tbody>
        </table>
        """
    end

    # Close HTML
    html *= """
        <div class="summary" style="margin-top: 40px;">
            <h2>Interpretation</h2>
            <p>
                The plots show <strong>Runtime/N vs N</strong>, which normalizes the runtime by the array size.
                For an optimal FFT implementation with O(N log N) complexity, this should scale as O(log N).
            </p>
            <p>
                <strong>Key Observations:</strong>
            </p>
            <ul>
                <li>Powers of 2 typically show the best performance for FFT algorithms due to the Cooley-Tukey radix-2 algorithm</li>
                <li>Powers of 3 may use different factorization strategies</li>
                <li>Composite numbers test the FFT implementation's ability to handle general factorizations</li>
                <li>FFTA.jl is a pure Julia implementation, while FFTW.jl wraps optimized C code</li>
            </ul>
        </div>

        <footer style="margin-top: 40px; padding: 20px; text-align: center; color: #7f8c8d; border-top: 1px solid #ddd;">
            <p>Generated on $(now())</p>
            <p>FFTA.jl Performance Benchmark Suite</p>
        </footer>
    </body>
    </html>
    """

    # Write HTML file
    output_file = joinpath(@__DIR__, "benchmark_report.html")
    write(output_file, html)
    println("HTML report saved to: $output_file")

    return output_file
end

# Generate report
using Dates
generate_html_report()
