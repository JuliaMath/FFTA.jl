using Documenter
using FFTA

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
