# Verify the strided-pencil copy-in fix: the dims=2 regression must go away,
# the dims=1 win must survive, allocations must stay zero, and results must
# still match FFTW.
using FFTA, FFTW, AbstractFFTs, LinearAlgebra, Printf
FFTW.set_num_threads(1)
function rplan(x, d)
    try
        invoke(AbstractFFTs.plan_rfft, Tuple{AbstractArray{eltype(x),2},Any}, x, d; num_threads = 1)
    catch
        invoke(AbstractFFTs.plan_rfft, Tuple{AbstractArray{eltype(x),2},FFTA.RegionTypes}, x, d; num_threads = 1)
    end
end
function tmin(f; reps=15)
    f(); minimum(@elapsed(f()) for _ in 1:reps)
end
println("| case | FFTA (ms) | alloc/exec | rel.err vs FFTW |")
println("|:--|--:|--:|--:|")
for (rows, cols) in ((64, 16384), (64, 65536)), d in (2, 1)
    for T in (Float64, Float32)
        x = randn(T, rows, cols)
        p = rplan(x, d)
        y = p * x
        usemul = try; mul!(y, p, x); true; catch; false; end
        f = usemul ? (() -> mul!(y, p, x)) : (() -> p * x)
        t = tmin(f)
        a = @allocated f()
        # FFTA and FFTW's plan_rfft methods are mutually ambiguous (the known
        # coexistence bug), so reach FFTW's StridedArray method explicitly.
        ref = invoke(AbstractFFTs.plan_rfft, Tuple{StridedArray{T,2},Any}, x, d) * x
        err = norm(y .- ref) / max(norm(ref), eps())
        @printf("| rfft %s %dx%d dims=%d | %.2f | %s | %.1e |\n", T, rows, cols, d, t*1e3,
                a == 0 ? "0" : @sprintf("%.1f KiB", a/1024), err)
    end
end
