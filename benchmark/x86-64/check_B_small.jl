# Back-to-back baseline-vs-B check for the small 1D rfft sizes where the
# cross-run comparison suggested B is slower. These are microsecond-scale, so
# cross-run drift dominates; same-session measurement is the only useful test.
using FFTA, AbstractFFTs, LinearAlgebra, Printf
function rp(x)
    try
        invoke(AbstractFFTs.plan_rfft, Tuple{AbstractArray{Float64,1},Any}, x, 1)
    catch
        invoke(AbstractFFTs.plan_rfft, Tuple{AbstractArray{Float64,1},FFTA.RegionTypes}, x, 1)
    end
end
function tmin(n; reps=2000)
    x = randn(Float64, n)
    p = rp(x)
    y = p * x
    usemul = try; mul!(y, p, x); true; catch; false; end
    f = usemul ? (() -> mul!(y, p, x)) : (() -> p * x)
    f(); (minimum(@elapsed(f()) for _ in 1:reps), usemul)
end
for n in (256, 1024, 2048, 4096)
    t, m = tmin(n)
    @printf("| %d | %.3f | %s |\n", n, t*1e6, m ? "mul!" : "*")
end
