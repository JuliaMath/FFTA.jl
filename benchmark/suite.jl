#!/usr/bin/env julia
#=
Comprehensive FFTA-vs-FFTW benchmark suite.

Both packages are loaded into the *same* process. Because FFTW.jl defines
`plan_*` methods on `StridedArray{ComplexF64}` etc. that are more specific
than FFTA's `AbstractArray` methods, calling `plan_fft(x)` would silently give
an FFTW plan. We therefore reach FFTA's methods explicitly with `invoke`
(see `ffta_plan_*` below). This also lets us cross-check every FFTA result
against FFTW in the same process.

Usage (from the `benchmark/` directory):

    julia --project=. -t 8 suite.jl [--quick] [--out results_suite.json]
                                    [--only 1d,nd,batched,threads] [--kinds fft,rfft]
                                    [--sizes 23,29,31]
                                    [--maxlog2 22] [--seconds 0.5]

Then `julia --project=. report.jl results_suite.json` renders REPORT.md.
=#

import Pkg
Pkg.develop(path = joinpath(@__DIR__, ".."))
Pkg.instantiate()

using AbstractFFTs
using FFTA
using FFTW
using BenchmarkTools
using JSON
using Primes
using Statistics
using LinearAlgebra
using Dates

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
function getopt(flag, default)
    i = findfirst(==(flag), ARGS)
    i === nothing && return default
    return ARGS[i + 1]
end
const QUICK    = "--quick" in ARGS
const OUTFILE  = getopt("--out", joinpath(@__DIR__, "results_suite.json"))
const ONLY     = split(getopt("--only", "1d,nd,batched,threads"), ",")
const KINDS    = Symbol.(split(getopt("--kinds", "fft,rfft"), ","))
# restrict the 1D sweep to these sizes (comma separated), e.g. to add rows to an existing run
const SIZES    = let v = getopt("--sizes", "")
    isempty(v) ? Int[] : parse.(Int, split(v, ","))
end
const MAXLOG2  = parse(Int, getopt("--maxlog2", QUICK ? "16" : "22"))
const SECONDS  = parse(Float64, getopt("--seconds", QUICK ? "0.2" : "0.5"))
const NTHREADS = Threads.nthreads()

# BenchmarkTools is used for `@allocated`-style checks and as a cross-check; the
# timing loop below avoids its per-callsite compilation cost (~0.5 s each),
# which dominates a sweep of several hundred cases.

# ---------------------------------------------------------------------------
# FFTA plan constructors that bypass FFTW's more specific methods
# ---------------------------------------------------------------------------
# `num_threads` is FFTA's keyword for the number of workers a plan may use
# (versions without it ignore unknown keywords).
ffta_plan_fft(x::AbstractArray{T,N}, dims; num_threads = 1) where {T<:Complex,N} =
    invoke(AbstractFFTs.plan_fft, Tuple{AbstractArray{T,N}, Any}, x, dims; num_threads)
ffta_plan_bfft(x::AbstractArray{T,N}, dims; num_threads = 1) where {T<:Complex,N} =
    invoke(AbstractFFTs.plan_bfft, Tuple{AbstractArray{T,N}, Any}, x, dims; num_threads)
# (older FFTA annotates `region::RegionTypes`; the signature tuple must match the
# method actually defined, so try the unannotated one first)
function ffta_plan_rfft(x::AbstractArray{T,N}, dims; num_threads = 1) where {T<:Real,N}
    try
        invoke(AbstractFFTs.plan_rfft, Tuple{AbstractArray{T,N}, Any}, x, dims; num_threads)
    catch
        invoke(AbstractFFTs.plan_rfft, Tuple{AbstractArray{T,N}, FFTA.RegionTypes}, x, dims; num_threads)
    end
end
function ffta_plan_brfft(x::AbstractArray{T,N}, len, dims; num_threads = 1) where {T<:Complex,N}
    try
        invoke(AbstractFFTs.plan_brfft, Tuple{AbstractArray{T,N}, Integer, Any}, x, len, dims; num_threads)
    catch
        invoke(AbstractFFTs.plan_brfft, Tuple{AbstractArray{T,N}, Int, FFTA.RegionTypes}, x, len, dims; num_threads)
    end
end

# With both packages loaded, `plan_rfft(::Vector{Float64}, ::Int)` is *ambiguous*
# (FFTW's `StridedArray` method vs FFTA's `region::RegionTypes` method), so FFTW's
# real plans have to be reached explicitly as well.
fftw_plan_fft(x::StridedArray{T,N}, dims; kw...) where {T<:Complex,N} =
    invoke(AbstractFFTs.plan_fft, Tuple{StridedArray{T,N}, Any}, x, dims; kw...)
fftw_plan_rfft(x::StridedArray{T,N}, dims; kw...) where {T<:Real,N} =
    invoke(AbstractFFTs.plan_rfft, Tuple{StridedArray{T,N}, Any}, x, dims; kw...)

# Sanity: these must really be FFTA plans
let x = randn(ComplexF64, 8), r = randn(Float64, 8)
    @assert ffta_plan_fft(x, 1)  isa FFTA.FFTAPlan
    @assert ffta_plan_rfft(r, 1) isa FFTA.FFTAPlan
    @assert fftw_plan_fft(x, 1)  isa FFTW.FFTWPlan
    @assert fftw_plan_rfft(r, 1) isa FFTW.FFTWPlan
end

# ---------------------------------------------------------------------------
# Size classes
# ---------------------------------------------------------------------------
function size_classes(maxlog2)
    nmax = 1 << maxlog2
    pow2   = [1 << k for k in 3:maxlog2]
    # 2^a 3^b 5^c 7^d, not powers of two, spread log-uniformly
    smooth_all = Int[]
    for a in 0:maxlog2, b in 0:14, c in 0:9, d in 0:8
        n = 2^a * 3^b * 5^c * 7^d
        (n <= nmax && n >= 8 && !ispow2(n)) && push!(smooth_all, n)
    end
    sort!(unique!(smooth_all))
    smooth = Int[]
    for t in exp.(range(log(12), log(nmax), length = 18))
        push!(smooth, smooth_all[argmin(abs.(log.(smooth_all) .- log(t)))])
    end
    # explicitly include some "classic" smooth sizes
    for n in (12, 60, 120, 360, 720, 1000, 1000_000)
        n <= nmax && push!(smooth, n)
    end
    sort!(unique!(smooth))
    # primes: below (DFT path) and above (Bluestein path) the cutoff of 73
    prime = Int[]
    for t in exp.(range(log(7), log(nmax), length = 16))
        push!(prime, nextprime(round(Int, t)))
    end
    # straddle the DFT/Bluestein crossover: the O(n²) leaf is used below the
    # cutoff (73 originally, 29 since the twiddle-table work), so sample the
    # band 23–46 as well as 61–79
    for n in (23, 29, 31, 37, 43, 61, 71, 73, 79)
        push!(prime, n)
    end
    filter!(<=(nmax), prime); sort!(unique!(prime))
    # awkward: (prime > cutoff) × small factor
    awk = Int[]
    for p in (101, 1009, 4099, 16411, 65537, 262147), f in (2, 3, 4, 6, 16)
        n = p * f
        n <= nmax && push!(awk, n)
    end
    sort!(unique!(awk))
    return (pow2 = pow2, smooth = smooth, prime = prime, awkward = awk)
end

const CLASSES = size_classes(MAXLOG2)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
const RESULTS = Dict{String,Any}[]

"""
    timeit(f; seconds) -> (min, median, samples)

Minimal benchmark loop: warm up, pick an eval count so each sample lasts at
least ~20 µs (timer granularity), then collect samples for `seconds`. Returns
per-call times in seconds.
"""
function timeit(f; seconds = SECONDS)
    f(); f()                                # compile + warm caches
    t0 = time_ns(); f(); t1 = time_ns()
    evals = max(1, ceil(Int, 20_000 / max(t1 - t0, 1)))
    times = Float64[]
    deadline = time_ns() + seconds * 1e9
    while time_ns() < deadline || length(times) < 5
        ts = time_ns()
        for _ in 1:evals
            f()
        end
        push!(times, (time_ns() - ts) / evals)
        length(times) >= 100_000 && break
    end
    return (min = minimum(times) * 1e-9, median = median(times) * 1e-9, samples = length(times))
end
bench_min(f)   = timeit(f)
bench_plan(mk) = timeit(mk)

relerr(a, b) = norm(a .- b) / max(norm(b), floatmin(real(eltype(b))))

# The output file is rewritten after every case so that a long run can be
# followed (and salvaged); `complete` is only set by the final write.
function record!(entry; complete::Bool = false)
    entry === nothing || push!(RESULTS, entry)
    open(OUTFILE, "w") do io
        JSON.print(io, Dict("meta" => META, "complete" => complete, "results" => RESULTS), 1)
    end
end

function class_of(n)
    for (k, v) in pairs(CLASSES)
        n in v && return String(k)
    end
    return "other"
end

sizestr(sz) = join(sz, "x")

# One planned transform benchmark (both packages) for a given problem.
# kind ∈ :fft, :rfft  ;  dims: Int or tuple ;  T: real element type
function bench_case!(kind::Symbol, T::Type, sz::Tuple, dims; shape::String,
                     fftw_threads::Int = 1, measure::Bool = false, extra = Dict())
    kind in KINDS || return nothing
    FFTW.set_num_threads(fftw_threads)
    x = kind === :fft ? randn(Complex{T}, sz) : randn(T, sz)
    n_along = prod(size(x, d) for d in dims)
    entry = Dict{String,Any}(
        "kind" => String(kind), "T" => string(T), "size" => collect(sz),
        "dims" => collect(dims), "shape" => shape, "n" => n_along,
        "class" => length(sz) == 1 ? class_of(sz[1]) : "nd",
        "fftw_threads" => fftw_threads, "fftw_flags" => measure ? "MEASURE" : "ESTIMATE",
    )
    merge!(entry, extra)
    flags = measure ? FFTW.MEASURE : FFTW.ESTIMATE

    # ---- FFTW
    mk_w = kind === :fft ? (() -> fftw_plan_fft(x, dims; flags)) :
                           (() -> fftw_plan_rfft(x, dims; flags))
    pw = mk_w()
    yw = pw * x
    if !measure
        tp = bench_plan(mk_w); entry["fftw_plan_min"] = tp.min; entry["fftw_plan_median"] = tp.median
    end
    te = bench_min(() -> mul!(yw, pw, x)); entry["fftw_exec_min"] = te.min; entry["fftw_exec_median"] = te.median
    entry["fftw_exec_alloc"] = (@allocated mul!(yw, pw, x))
    entry["fftw_oneshot_min"] = bench_min(() -> mk_w() * x).min

    # ---- FFTA
    try
        mk_a = kind === :fft ? (() -> ffta_plan_fft(x, dims; num_threads = fftw_threads)) :
                               (() -> ffta_plan_rfft(x, dims; num_threads = fftw_threads))
        pa = mk_a()
        ya = pa * x
        entry["ffta_relerr"] = relerr(ya, yw)
        tp = bench_plan(mk_a); entry["ffta_plan_min"] = tp.min; entry["ffta_plan_median"] = tp.median
        # (`applicable`/`hasmethod` are fooled by LinearAlgebra's generic
        # 3-argument `mul!` fallback, so try the call)
        has_mul = try
            mul!(ya, pa, x); true
        catch err
            err isa MethodError || rethrow()
            false
        end
        if has_mul
            te = bench_min(() -> mul!(ya, pa, x))
            entry["ffta_exec_alloc"] = (@allocated mul!(ya, pa, x))
            entry["ffta_exec_api"] = "mul!"
        else
            # older FFTA real plans only implement `*` (no mul!)
            te = bench_min(() -> pa * x)
            entry["ffta_exec_alloc"] = (@allocated pa * x)
            entry["ffta_exec_api"] = "*"
        end
        entry["ffta_exec_min"] = te.min; entry["ffta_exec_median"] = te.median
        entry["ffta_oneshot_min"] = bench_min(() -> mk_a() * x).min
    catch err
        entry["ffta_error"] = sprint(showerror, err)[1:min(end, 200)]
        @warn "FFTA failed" kind T sz dims err = entry["ffta_error"]
    end
    record!(entry)
    r = haskey(entry, "ffta_exec_min") ? round(entry["ffta_exec_min"] / entry["fftw_exec_min"]; digits = 2) : "n/a"
    println(rpad("$kind $(T) $(sizestr(sz)) dims=$dims [$(shape)] thr=$fftw_threads", 60),
            " fftw=", round(entry["fftw_exec_min"] * 1e6; digits = 2), "us  ratio=", r)
    flush(stdout)   # progress lines are block-buffered when stdout is a file
    return entry
end

# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
cpuinfo() = try
    s = read("/proc/cpuinfo", String)
    m = match(r"model name\s*:\s*(.*)", s)
    if m === nothing   # aarch64 /proc/cpuinfo has no model name; ask lscpu
        m = match(r"Model name:\s*(.*)", read(`lscpu`, String))
    end
    m === nothing ? Sys.CPU_NAME : strip(String(m.captures[1])) * " ($(Sys.CPU_NAME))"
catch
    Sys.CPU_NAME
end
const META = Dict(
    "date" => string(now(UTC)), "julia" => string(VERSION), "arch" => String(Sys.ARCH),
    "os" => string(Sys.KERNEL), "cpu" => cpuinfo(), "cpu_name" => Sys.CPU_NAME,
    "ncores" => Sys.CPU_THREADS, "julia_threads" => NTHREADS,
    "fftw_version" => string(FFTW.version), "fftw_provider" => FFTW.get_provider(),
    "ffta_version" => string(Pkg.dependencies()[Base.UUID("b86e33f2-c0db-4aa1-a6e0-ab43e668529e")].version),
    "ffta_sha" => try strip(read(`git -C $(joinpath(@__DIR__, "..")) rev-parse --short HEAD`, String)) catch; "?" end,
    "benchmark_seconds" => SECONDS, "quick" => QUICK, "maxlog2" => MAXLOG2,
    "sizes" => Dict(String(k) => v for (k, v) in pairs(CLASSES)),
)
println("== FFTA vs FFTW benchmark suite")
for (k, v) in META
    k == "sizes" && continue
    println("   $k = $v")
end
for (k, v) in pairs(CLASSES)
    println("   sizes.$k = $v")
end

# ---------------------------------------------------------------------------
# 1. 1D sweep: all classes × {ComplexF64, ComplexF32} fft and {Float64, Float32} rfft
# ---------------------------------------------------------------------------
if "1d" in ONLY
    println("\n== 1D sweep")
    for cls in (:pow2, :smooth, :prime, :awkward), n in CLASSES[cls]
        isempty(SIZES) || n in SIZES || continue
        for T in (Float64, Float32)
            bench_case!(:fft,  T, (n,), 1; shape = "1d")
            bench_case!(:rfft, T, (n,), 1; shape = "1d")
        end
    end
    # FFTW with MEASURE planning, ComplexF64 pow2 only (what a tuned consumer would see)
    println("\n== 1D pow2 ComplexF64 with FFTW.MEASURE")
    for n in CLASSES.pow2
        isempty(SIZES) || n in SIZES || continue
        bench_case!(:fft, Float64, (n,), 1; shape = "1d", measure = true)
    end
end

# ---------------------------------------------------------------------------
# 2. Multidimensional: 2D n×n, 2D non-square, 3D n³
# ---------------------------------------------------------------------------
if "nd" in ONLY
    println("\n== 2D / 3D")
    nd2 = filter(n -> n * n <= (1 << MAXLOG2), [8, 16, 32, 64, 128, 256, 512, 1024, 2048])
    for n in nd2, T in (Float64, Float32)
        bench_case!(:fft,  T, (n, n), (1, 2); shape = "2d")
        bench_case!(:rfft, T, (n, n), (1, 2); shape = "2d")
    end
    for sz in ((1000, 1000), (720, 480), (1009, 64), (64, 1009), (127, 257))
        prod(sz) <= (1 << MAXLOG2) || continue
        bench_case!(:fft,  Float64, sz, (1, 2); shape = "2d")
        bench_case!(:rfft, Float64, sz, (1, 2); shape = "2d")
    end
    nd3 = filter(n -> n^3 <= (1 << MAXLOG2), [8, 16, 32, 64, 128])
    for n in nd3, T in (Float64, Float32)
        bench_case!(:fft,  T, (n, n, n), (1, 2, 3); shape = "3d")
        bench_case!(:rfft, T, (n, n, n), (1, 2, 3); shape = "3d")
    end
end

# ---------------------------------------------------------------------------
# 3. Batched: transform along one dim of a matrix (dims=1 contiguous, dims=2 strided)
# ---------------------------------------------------------------------------
if "batched" in ONLY
    println("\n== batched (dims keyword)")
    for n in filter(n -> n * 64 <= (1 << MAXLOG2), [64, 256, 1024, 4096, 16384, 65536]),
        T in (Float64, Float32)
        bench_case!(:fft,  T, (n, 64), 1; shape = "batched_dim1")
        bench_case!(:fft,  T, (64, n), 2; shape = "batched_dim2")
        bench_case!(:rfft, T, (n, 64), 1; shape = "batched_dim1")
        bench_case!(:rfft, T, (64, n), 2; shape = "batched_dim2")
    end
    # a Welch/periodogram-like shape: many short windows
    for (n, m) in ((256, 4096), (1024, 1024), (4096, 256))
        n * m <= (1 << MAXLOG2) || continue
        bench_case!(:rfft, Float64, (n, m), 1; shape = "batched_dim1")
    end
end

# ---------------------------------------------------------------------------
# 4. Threading: FFTW multi-threaded vs single (FFTA has no threading today)
# ---------------------------------------------------------------------------
if "threads" in ONLY && NTHREADS > 1
    println("\n== FFTW threading (julia -t $NTHREADS)")
    for n in filter(n -> n >= 1 << 16, CLASSES.pow2)
        bench_case!(:fft, Float64, (n,), 1; shape = "1d", fftw_threads = NTHREADS)
    end
    for n in filter(n -> n * n <= (1 << MAXLOG2) && n >= 256, [256, 512, 1024, 2048])
        bench_case!(:fft, Float64, (n, n), (1, 2); shape = "2d", fftw_threads = NTHREADS)
    end
    for n in filter(n -> n * 64 <= (1 << MAXLOG2) && n >= 1024, [1024, 4096, 16384, 65536])
        bench_case!(:fft, Float64, (n, 64), 1; shape = "batched_dim1", fftw_threads = NTHREADS)
    end
    FFTW.set_num_threads(1)
end

record!(nothing; complete = true)
println("\nWrote $(length(RESULTS)) results to $OUTFILE")
