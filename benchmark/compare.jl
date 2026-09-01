#!/usr/bin/env julia
# Compare FFTA execution times between two `suite.jl` result files.
#
#     julia --project=. compare.jl before.json after.json [--all]
#
# Prints a markdown table of FFTA exec time before/after, the speedup, and the
# FFTA/FFTW ratio after (from the `after` file). Cases are matched on
# (kind, T, size, dims, fftw_threads, fftw_flags). By default only cases where
# the FFTA time changed by more than 5% are listed; `--all` lists everything.
import Pkg; Pkg.instantiate()
using JSON, Printf, Statistics

files = filter(a -> endswith(a, ".json"), ARGS)
length(files) == 2 || error("usage: compare.jl before.json after.json [--all]")
showall = "--all" in ARGS
key(e) = (e["kind"], e["T"], join(e["size"], "×"), join(e["dims"], ","), e["fftw_threads"], e["fftw_flags"])
before = Dict(key(e) => e for e in JSON.parsefile(files[1])["results"])
after  = JSON.parsefile(files[2])["results"]

fmt_t(t) = t < 1e-6 ? @sprintf("%.0f ns", t * 1e9) : t < 1e-3 ? @sprintf("%.2f µs", t * 1e6) :
           t < 1 ? @sprintf("%.3f ms", t * 1e3) : @sprintf("%.3f s", t)
fmt_b(b) = b == 0 ? "0" : b < 1024 ? "$(b) B" : b < 1024^2 ? @sprintf("%.1f KiB", b / 1024) : @sprintf("%.2f MiB", b / 1024^2)

println("| case | FFTA before | FFTA after | speedup | alloc before → after | FFTA/FFTW after |")
println("|:--|--:|--:|--:|:--|--:|")
speedups = Float64[]
for e in after
    k = key(e); haskey(before, k) || continue
    b = before[k]
    (haskey(b, "ffta_exec_min") && haskey(e, "ffta_exec_min")) || continue
    sp = b["ffta_exec_min"] / e["ffta_exec_min"]; push!(speedups, sp)
    (showall || abs(sp - 1) > 0.05) || continue
    @printf("| %s %s %s dims=%s%s | %s | %s | **%.2f×** | %s → %s | %.2f× |\n",
            e["kind"], e["T"], k[3], k[4], e["fftw_threads"] > 1 ? " ($(e["fftw_threads"]) thr)" : "",
            fmt_t(b["ffta_exec_min"]), fmt_t(e["ffta_exec_min"]), sp,
            fmt_b(b["ffta_exec_alloc"]), fmt_b(e["ffta_exec_alloc"]), e["ffta_exec_min"] / e["fftw_exec_min"])
end
isempty(speedups) || @printf("\n%d matched cases; geometric-mean speedup %.2f× (min %.2f×, max %.2f×)\n",
                             length(speedups), exp(mean(log.(speedups))), minimum(speedups), maximum(speedups))
