#!/usr/bin/env julia
# Render `REPORT.md` (+ SVG plots in `plots/`) from one or more `suite.jl` JSON
# result files.
#
#     julia --project=. report.jl results_suite.json [more.json ...] [--out REPORT.md]
#
# Multiple files are concatenated (e.g. a single-threaded run and a `-t 8` run).
# No plotting package is needed: plots are written as plain SVG.
import Pkg; Pkg.instantiate()
using JSON, Statistics, Printf, Dates

function getopt(flag, default)
    i = findfirst(==(flag), ARGS)
    i === nothing ? default : ARGS[i + 1]
end
const OUT = getopt("--out", joinpath(@__DIR__, "REPORT.md"))
const PLOTDIR = joinpath(dirname(OUT), "plots")
mkpath(PLOTDIR)
const INPUTS = filter(a -> endswith(a, ".json"), ARGS)
isempty(INPUTS) && error("give at least one results JSON")

metas = Dict{String,Any}[]
results = Dict{String,Any}[]
for f in INPUTS
    d = JSON.parsefile(f)
    push!(metas, d["meta"]); append!(results, d["results"])
end
const META = first(metas)

# ---------------------------------------------------------------------------
# formatting helpers
# ---------------------------------------------------------------------------
fmt_t(t) = t === nothing ? "—" : t < 1e-6 ? @sprintf("%.0f ns", t * 1e9) :
           t < 1e-3 ? @sprintf("%.2f µs", t * 1e6) : t < 1 ? @sprintf("%.3f ms", t * 1e3) :
           @sprintf("%.3f s", t)
fmt_r(r) = r === nothing ? "—" : @sprintf("%.2f×", r)
fmt_b(b) = b === nothing ? "—" : b == 0 ? "0" : b < 1024 ? "$(b) B" : b < 1024^2 ?
           @sprintf("%.1f KiB", b / 1024) : @sprintf("%.2f MiB", b / 1024^2)
fmt_e(e) = e === nothing ? "—" : @sprintf("%.1e", e)
sizestr(e) = join(e["size"], "×")
ratio(e) = haskey(e, "ffta_exec_min") ? e["ffta_exec_min"] / e["fftw_exec_min"] : nothing
cold_ratio(e) = haskey(e, "ffta_oneshot_min") ? e["ffta_oneshot_min"] / e["fftw_oneshot_min"] : nothing
plan_ratio(e) = haskey(e, "ffta_plan_min") && haskey(e, "fftw_plan_min") ? e["ffta_plan_min"] / e["fftw_plan_min"] : nothing
label(e) = (e["kind"] == "fft" ? "Complex" * replace(e["T"], "Float" => "F") : e["T"]) * " " * e["kind"]
gmean(v) = isempty(v) ? nothing : exp(mean(log.(v)))

sel(; kw...) = filter(results) do e
    all(e[string(k)] == v for (k, v) in kw)
end
sel1(; kw...) = filter(e -> e["fftw_threads"] == 1 && e["fftw_flags"] == "ESTIMATE", sel(; kw...))

# ---------------------------------------------------------------------------
# SVG log-log plot
# ---------------------------------------------------------------------------
const COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#17becf", "#7f7f7f"]
function svg_loglog(path; series, title, xlabel, ylabel, hline = nothing, width = 860, height = 480)
    ml, mr, mt, mb = 70, 20, 40, 55
    pw, ph = width - ml - mr, height - mt - mb
    xs = vcat((s.x for s in series)...); ys = vcat((s.y for s in series)...)
    xs = filter(>(0), xs); ys = filter(>(0), ys)
    (isempty(xs) || isempty(ys)) && return
    lx0, lx1 = log10(minimum(xs)), log10(maximum(xs))
    ly0, ly1 = log10(minimum(ys)), log10(maximum(ys))
    hline !== nothing && (ly0 = min(ly0, log10(hline)); ly1 = max(ly1, log10(hline)))
    ly0 = floor(ly0 * 2) / 2 ; ly1 = ceil(ly1 * 2) / 2
    lx0 = floor(lx0); lx1 = ceil(lx1)
    X(v) = ml + (log10(v) - lx0) / (lx1 - lx0) * pw
    Y(v) = mt + ph - (log10(v) - ly0) / (ly1 - ly0) * ph
    io = IOBuffer()
    print(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height" font-family="Helvetica, Arial, sans-serif" font-size="12">
<rect width="100%" height="100%" fill="white"/>
<text x="$(width ÷ 2)" y="22" text-anchor="middle" font-size="15" font-weight="bold">$title</text>
""")
    for d in ceil(Int, lx0):floor(Int, lx1)
        x = X(10.0^d)
        print(io, """<line x1="$x" y1="$mt" x2="$x" y2="$(mt + ph)" stroke="#ddd"/><text x="$x" y="$(mt + ph + 16)" text-anchor="middle">1e$d</text>\n""")
    end
    for d in ceil(Int, ly0 * 2):floor(Int, ly1 * 2)
        v = 10.0^(d / 2); y = Y(v)
        lab = isinteger(d / 2) ? "1e$(d ÷ 2)" : @sprintf("%.2g", v)
        print(io, """<line x1="$ml" y1="$y" x2="$(ml + pw)" y2="$y" stroke="#eee"/><text x="$(ml - 6)" y="$(y + 4)" text-anchor="end">$lab</text>\n""")
    end
    if hline !== nothing
        y = Y(hline)
        print(io, """<line x1="$ml" y1="$y" x2="$(ml + pw)" y2="$y" stroke="#000" stroke-dasharray="6,4"/>\n""")
    end
    print(io, """<rect x="$ml" y="$mt" width="$pw" height="$ph" fill="none" stroke="#333"/>
<text x="$(ml + pw ÷ 2)" y="$(height - 12)" text-anchor="middle">$xlabel</text>
<text transform="translate(16,$(mt + ph ÷ 2)) rotate(-90)" text-anchor="middle">$ylabel</text>
""")
    for (i, s) in enumerate(series)
        c = COLORS[mod1(i, length(COLORS))]
        pts = [(X(x), Y(y)) for (x, y) in zip(s.x, s.y) if x > 0 && y > 0]
        isempty(pts) && continue
        print(io, """<polyline fill="none" stroke="$c" stroke-width="1.5" points="$(join(("$(round(x; digits=1)),$(round(y; digits=1))" for (x, y) in pts), ' '))"/>\n""")
        for (x, y) in pts
            print(io, """<circle cx="$(round(x; digits=1))" cy="$(round(y; digits=1))" r="3" fill="$c"/>""")
        end
        ly = mt + 10 + 16 * (i - 1)
        print(io, """\n<rect x="$(ml + 10)" y="$(ly - 5)" width="12" height="10" fill="$c"/><text x="$(ml + 26)" y="$(ly + 4)">$(s.name)</text>\n""")
    end
    print(io, "</svg>\n")
    write(path, take!(io))
end

# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------
function table_1d(io, rows; cold = true)
    println(io, "| n | class | FFTW exec | FFTA exec | **exec ratio** | FFTW plan | FFTA plan | cold ratio | FFTA alloc/exec | rel. err |")
    println(io, "|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for e in rows
        if haskey(e, "ffta_error")
            println(io, "| $(e["n"]) | $(e["class"]) | $(fmt_t(e["fftw_exec_min"])) | ✗ `$(e["ffta_error"][1:min(end,60)])` | — | $(fmt_t(get(e,"fftw_plan_min",nothing))) | — | — | — | — |")
            continue
        end
        println(io, "| $(e["n"]) | $(e["class"]) | $(fmt_t(e["fftw_exec_min"])) | $(fmt_t(e["ffta_exec_min"])) | **$(fmt_r(ratio(e)))** | $(fmt_t(get(e,"fftw_plan_min",nothing))) | $(fmt_t(e["ffta_plan_min"])) | $(fmt_r(cold_ratio(e))) | $(fmt_b(e["ffta_exec_alloc"])) | $(fmt_e(e["ffta_relerr"])) |")
    end
end
function table_nd(io, rows)
    println(io, "| size | dims | FFTW exec | FFTA exec | **exec ratio** | FFTW plan | FFTA plan | cold ratio | FFTW alloc | FFTA alloc | rel. err |")
    println(io, "|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    for e in rows
        d = join(e["dims"], ",")
        if haskey(e, "ffta_error")
            println(io, "| $(sizestr(e)) | $d | $(fmt_t(e["fftw_exec_min"])) | ✗ `$(e["ffta_error"][1:min(end,60)])` | — | $(fmt_t(get(e,"fftw_plan_min",nothing))) | — | — | $(fmt_b(e["fftw_exec_alloc"])) | — | — |")
            continue
        end
        println(io, "| $(sizestr(e)) | $d | $(fmt_t(e["fftw_exec_min"])) | $(fmt_t(e["ffta_exec_min"])) | **$(fmt_r(ratio(e)))** | $(fmt_t(get(e,"fftw_plan_min",nothing))) | $(fmt_t(e["ffta_plan_min"])) | $(fmt_r(cold_ratio(e))) | $(fmt_b(e["fftw_exec_alloc"])) | $(fmt_b(e["ffta_exec_alloc"])) | $(fmt_e(e["ffta_relerr"])) |")
    end
end

# ---------------------------------------------------------------------------
# write report
# ---------------------------------------------------------------------------
open(OUT, "w") do io
    println(io, "# FFTA.jl vs FFTW.jl benchmark report\n")
    println(io, "Generated $(Dates.format(now(UTC), "yyyy-mm-dd HH:MM")) UTC by `benchmark/report.jl` from `", join(basename.(INPUTS), "`, `"), "`.\n")
    println(io, "## Environment\n")
    println(io, "| | |\n|:--|:--|")
    for (k, lab) in (("cpu", "CPU"), ("arch", "Architecture"), ("ncores", "Cores"), ("os", "OS"),
                     ("julia", "Julia"), ("julia_threads", "Julia threads"), ("ffta_version", "FFTA version"),
                     ("ffta_sha", "FFTA commit"), ("fftw_version", "FFTW version"), ("fftw_provider", "FFTW provider"),
                     ("benchmark_seconds", "Time budget per measurement (s)"), ("date", "Run date (UTC)"))
        println(io, "| $lab | $(get(META, k, "?")) |")
    end
    println(io, """

## Methodology

* Both packages are loaded in one Julia process; FFTA plans are created via
  `invoke` on FFTA's `AbstractFFTs.plan_*` methods (FFTW's more specific
  `StridedArray` methods would otherwise be selected). Every FFTA result is
  compared with FFTW's (`rel. err` = ‖y_FFTA − y_FFTW‖ / ‖y_FFTW‖).
* **exec** = execution of a pre-built plan (`mul!(y, p, x)` where supported;
  if a plan only implements `p * x`, that is timed instead and includes the
  output allocation — see the *FFTA API* column of the allocation table). Times are the **minimum** over samples collected during the
  time budget (batched evals for sub-20 µs calls).
* **plan** = time to construct a plan (FFTW with `FFTW.ESTIMATE`, the
  FFTW.jl default, unless stated otherwise). **cold ratio** = (FFTA plan+exec)
  / (FFTW plan+exec), i.e. the one-shot `fft(x)` cost ratio.
* **exec ratio** = FFTA exec / FFTW exec — values > 1 mean FFTA is slower.
* **alloc/exec** = bytes allocated by one planned execution (`@allocated`).
* FFTW and FFTA are single-threaded (FFTA: one worker) except in the
  *Threading* section, where both are planned with the same thread count
  (FFTA versions with a `num_threads` keyword thread across the pencils of
  multidimensional and batched transforms).
* Size classes: `pow2` = 2^k; `smooth` = 2^a·3^b·5^c·7^d (not a power of 2);
  `prime`; `awkward` = (prime > 73) × small factor. Primes < 73 use FFTA's
  O(n²) DFT leaf, primes ≥ 73 use Bluestein.
""")

    # ---- summary
    println(io, "## Summary: geometric-mean exec ratio (FFTA / FFTW), single-threaded\n")
    println(io, "| type | pow2 | smooth | prime | awkward | 2D | 3D | batched dim=1 | batched dim=2 |")
    println(io, "|:--|--:|--:|--:|--:|--:|--:|--:|--:|")
    for kind in ("fft", "rfft"), T in ("Float64", "Float32")
        cells = String[]
        for cls in ("pow2", "smooth", "prime", "awkward")
            rs = filter(!isnothing, ratio.(sel1(kind = kind, T = T, shape = "1d", class = cls)))
            g = gmean(rs)
            push!(cells, g === nothing ? "—" : @sprintf("%.2f× (%.1f–%.1f)", g, minimum(rs), maximum(rs)))
        end
        for sh in ("2d", "3d", "batched_dim1", "batched_dim2")
            es = sel1(kind = kind, T = T, shape = sh)
            rs = filter(!isnothing, ratio.(es))
            g = gmean(rs)
            nerr = count(e -> haskey(e, "ffta_error"), es)
            push!(cells, g === nothing ? (nerr > 0 ? "✗ unsupported" : "—") :
                  @sprintf("%.2f× (%.1f–%.1f)%s", g, minimum(rs), maximum(rs), nerr > 0 ? " ($nerr ✗)" : ""))
        end
        println(io, "| $(label(Dict("kind"=>kind,"T"=>T))) | ", join(cells, " | "), " |")
    end
    println(io, "\nCells are `geomean (min–max)`. ✗ = FFTA threw (unsupported).\n")

    # ---- plots
    println(io, "## Ratio plots (single-threaded, planned execution)\n")
    for (kind, T) in (("fft", "Float64"), ("rfft", "Float64"), ("fft", "Float32"), ("rfft", "Float32"))
        series = []
        for cls in ("pow2", "smooth", "prime", "awkward")
            es = sort(filter(e -> ratio(e) !== nothing, sel1(kind = kind, T = T, shape = "1d", class = cls)); by = e -> e["n"])
            isempty(es) && continue
            push!(series, (name = cls, x = [e["n"] for e in es], y = ratio.(es)))
        end
        fn = "ratio_1d_$(kind)_$(T).svg"
        svg_loglog(joinpath(PLOTDIR, fn); series, title = "1D $(label(Dict("kind"=>kind,"T"=>T))): FFTA exec time / FFTW exec time",
                   xlabel = "n", ylabel = "ratio (FFTA / FFTW)", hline = 1.0)
        println(io, "![]($(joinpath("plots", fn)))\n")
    end
    # absolute time per element for ComplexF64
    let series = []
        for (pkg, key) in (("FFTW", "fftw_exec_min"), ("FFTA", "ffta_exec_min")), cls in ("pow2", "smooth", "prime")
            es = sort(filter(e -> haskey(e, key), sel1(kind = "fft", T = "Float64", shape = "1d", class = cls)); by = e -> e["n"])
            isempty(es) && continue
            push!(series, (name = "$pkg $cls", x = [e["n"] for e in es], y = [e[key] / (e["n"] * log2(e["n"])) * 1e9 for e in es]))
        end
        fn = "time_per_nlogn_1d_fft_Float64.svg"
        svg_loglog(joinpath(PLOTDIR, fn); series, title = "1D ComplexF64 fft: exec time / (n log2 n)", xlabel = "n", ylabel = "ns per n·log2(n)")
        println(io, "![]($(joinpath("plots", fn)))\n")
    end
    let series = []
        for (pkg, key) in (("FFTW", "fftw_plan_min"), ("FFTA", "ffta_plan_min")), cls in ("pow2", "smooth", "prime")
            es = sort(filter(e -> haskey(e, key), sel1(kind = "fft", T = "Float64", shape = "1d", class = cls)); by = e -> e["n"])
            isempty(es) && continue
            push!(series, (name = "$pkg $cls", x = [e["n"] for e in es], y = [e[key] for e in es]))
        end
        fn = "plan_time_1d_fft_Float64.svg"
        svg_loglog(joinpath(PLOTDIR, fn); series, title = "1D ComplexF64: plan creation time (FFTW.ESTIMATE)", xlabel = "n", ylabel = "seconds")
        println(io, "![]($(joinpath("plots", fn)))\n")
    end

    # ---- 1D tables
    println(io, "## 1D transforms\n")
    for kind in ("fft", "rfft"), T in ("Float64", "Float32")
        rows = sort(sel1(kind = kind, T = T, shape = "1d"); by = e -> e["n"])
        isempty(rows) && continue
        open_ = T == "Float64"
        println(io, "<details$(open_ ? " open" : "")><summary><b>1D $(label(rows[1]))</b> ($(length(rows)) sizes)</summary>\n")
        table_1d(io, rows)
        println(io, "\n</details>\n")
    end
    let rows = sort(filter(e -> e["fftw_flags"] == "MEASURE", sel(kind = "fft", T = "Float64", shape = "1d")); by = e -> e["n"])
        if !isempty(rows)
            println(io, "<details><summary><b>1D ComplexF64 fft vs FFTW planned with FFTW.MEASURE</b> (pow2 only)</summary>\n")
            println(io, "| n | FFTW exec (MEASURE) | FFTA exec | ratio |\n|--:|--:|--:|--:|")
            for e in rows
                println(io, "| $(e["n"]) | $(fmt_t(e["fftw_exec_min"])) | $(fmt_t(get(e, "ffta_exec_min", nothing))) | $(fmt_r(ratio(e))) |")
            end
            println(io, "\n</details>\n")
        end
    end

    # ---- ND
    println(io, "## 2D and 3D transforms (all dims)\n")
    for sh in ("2d", "3d"), kind in ("fft", "rfft"), T in ("Float64", "Float32")
        rows = sort(sel1(kind = kind, T = T, shape = sh); by = e -> prod(e["size"]))
        isempty(rows) && continue
        println(io, "<details$(T == "Float64" ? " open" : "")><summary><b>$(uppercase(sh)) $(label(rows[1]))</b></summary>\n")
        table_nd(io, rows); println(io, "\n</details>\n")
    end

    # ---- batched
    println(io, "## Batched transforms along one dimension of a matrix (`dims` keyword)\n")
    println(io, "`dims=1` transforms contiguous columns; `dims=2` transforms strided rows.\n")
    for kind in ("fft", "rfft"), T in ("Float64", "Float32")
        rows = sort(filter(e -> startswith(e["shape"], "batched"), sel1(kind = kind, T = T)); by = e -> (e["dims"][1], prod(e["size"])))
        isempty(rows) && continue
        println(io, "<details$(T == "Float64" ? " open" : "")><summary><b>Batched $(label(rows[1]))</b></summary>\n")
        table_nd(io, rows); println(io, "\n</details>\n")
    end

    # ---- threading
    let rows = filter(e -> e["fftw_threads"] > 1, results)
        if !isempty(rows)
            nt = rows[1]["fftw_threads"]
            println(io, "## Threading: FFTW vs FFTA with $nt threads\n")
            println(io, "| size | dims | FFTW 1 thread | FFTW $nt threads | FFTW speedup | FFTA 1 thread | FFTA $nt threads | FFTA speedup | FFTA / FFTW ($nt thr) |\n|:--|:--|--:|--:|--:|--:|--:|--:|--:|")
            for e in sort(rows; by = e -> (e["shape"], prod(e["size"])))
                base = filter(b -> b["kind"] == e["kind"] && b["T"] == e["T"] && b["size"] == e["size"] && b["dims"] == e["dims"] &&
                                   b["fftw_threads"] == 1 && b["fftw_flags"] == "ESTIMATE", results)
                t1 = isempty(base) ? nothing : base[1]["fftw_exec_min"]
                a1 = isempty(base) ? nothing : get(base[1], "ffta_exec_min", nothing)
                tn = e["fftw_exec_min"]; ta = get(e, "ffta_exec_min", nothing)
                println(io, "| $(sizestr(e)) | $(join(e["dims"], ",")) | $(fmt_t(t1)) | $(fmt_t(tn)) | $(t1 === nothing ? "—" : fmt_r(t1 / tn)) | $(fmt_t(a1)) | $(fmt_t(ta)) | $(a1 === nothing || ta === nothing ? "—" : fmt_r(a1 / ta)) | $(ta === nothing ? "—" : fmt_r(ta / tn)) |")
            end
            println(io)
        end
    end

    # ---- allocations
    println(io, "## Allocations per planned execution\n")
    println(io, "| case | FFTA API | zero-alloc? | typical bytes |\n|:--|:--|:--|--:|")
    for (name, f) in (("1D complex fft, mul!", e -> e["kind"] == "fft" && e["shape"] == "1d" && e["class"] != "prime" && e["class"] != "awkward"),
                      ("1D complex fft, prime/awkward (Bluestein)", e -> e["kind"] == "fft" && e["shape"] == "1d" && (e["class"] == "prime" || e["class"] == "awkward") && e["n"] >= 73),
                      ("1D rfft", e -> e["kind"] == "rfft" && e["shape"] == "1d"),
                      ("2D/3D complex fft, mul!", e -> e["kind"] == "fft" && e["shape"] in ("2d", "3d")),
                      ("2D rfft", e -> e["kind"] == "rfft" && e["shape"] == "2d"),
                      ("batched complex fft, mul!", e -> e["kind"] == "fft" && startswith(e["shape"], "batched")),
                      ("batched rfft", e -> e["kind"] == "rfft" && startswith(e["shape"], "batched")))
        es = filter(e -> f(e) && haskey(e, "ffta_exec_alloc") && e["fftw_threads"] == 1, results)
        isempty(es) && continue
        al = [e["ffta_exec_alloc"] for e in es]
        println(io, "| $name | `$(es[1]["ffta_exec_api"])` | $(all(==(0), al) ? "yes" : "no") | $(fmt_b(round(Int, median(al)))) (max $(fmt_b(maximum(al)))) |")
    end
    println(io)
    println(io, "FFTW planned execution allocates 0 bytes in all cases", any(e -> e["fftw_exec_alloc"] != 0, results) ? " except: " * join(("$(sizestr(e)) $(e["kind"])" for e in filter(e -> e["fftw_exec_alloc"] != 0, results)), ", ") : ".", "\n")
end
println("wrote $OUT and plots in $PLOTDIR")
