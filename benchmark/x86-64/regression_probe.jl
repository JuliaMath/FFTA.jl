# Localise the x86-only 64xN dims=2 regression seen in integration/all.
# Run in each branch's env; compare the same rows across branches.
# dims=1 at the same sizes is the control: it did not regress.
using FFTA, FFTW, AbstractFFTs, LinearAlgebra, Printf

FFTW.set_num_threads(1)
ffta_plan(x, d) = invoke(AbstractFFTs.plan_fft, Tuple{AbstractArray{eltype(x),2},Any}, x, d; num_threads = 1)
function ffta_rplan(x, d)
    try
        invoke(AbstractFFTs.plan_rfft, Tuple{AbstractArray{eltype(x),2},Any}, x, d; num_threads = 1)
    catch
        invoke(AbstractFFTs.plan_rfft, Tuple{AbstractArray{eltype(x),2},FFTA.RegionTypes}, x, d; num_threads = 1)
    end
end
function tmin(f; reps=8)
    f(); minimum(@elapsed(f()) for _ in 1:reps)
end

println("| case | FFTA (ms) |")
println("|:--|--:|")
for (rows, cols) in ((64, 16384), (64, 65536)), d in (2, 1)
    for T in (ComplexF64, ComplexF32)
        x = randn(T, rows, cols)
        p = ffta_plan(x, d); y = p * x
        @printf("| fft %s %dx%d dims=%d | %.2f |\n", T, rows, cols, d, tmin(() -> mul!(y, p, x))*1e3)
    end
    for T in (Float64, Float32)
        x = randn(T, rows, cols)
        p = ffta_rplan(x, d)
        t = try
            y = p * x
            try
                mul!(y, p, x); tmin(() -> mul!(y, p, x))
            catch
                tmin(() -> p * x)      # older branches: no mul! for real plans
            end
        catch
            NaN
        end
        @printf("| rfft %s %dx%d dims=%d | %.2f |\n", T, rows, cols, d, t*1e3)
    end
end
