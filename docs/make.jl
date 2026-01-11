using Documenter
using FFTA

# Check for and copy local benchmark results before building docs
benchmark_src = joinpath(@__DIR__, "..", "benchmark", "benchmark_report.html")
benchmark_dest_dir = joinpath(@__DIR__, "src", "assets", "benchmarks")
benchmark_dest = joinpath(benchmark_dest_dir, "benchmark_report.html")

if isfile(benchmark_src)
    @info "Found local benchmark results, copying to docs/src/assets/benchmarks/"
    mkpath(benchmark_dest_dir)
    cp(benchmark_src, benchmark_dest, force=true)

    # Also copy JSON files if they exist
    for json_file in ["results_ffta.json", "results_fftw.json"]
        src = joinpath(@__DIR__, "..", "benchmark", json_file)
        if isfile(src)
            cp(src, joinpath(benchmark_dest_dir, json_file), force=true)
        end
    end
    @info "Benchmark results copied successfully"
else
    @info "No local benchmark results found - documentation will show placeholder"
end

makedocs(
    sitename = "FFTA",
    format = Documenter.HTML(),
    pages = [
        "Development Tools" => "dev.md",
        "Benchmarks" => "benchmarks.md"
    ],
    modules = [FFTA]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JuliaMath/FFTA.jl.git",
    push_preview = true
)
