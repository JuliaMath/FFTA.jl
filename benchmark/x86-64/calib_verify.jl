using FFTA, AbstractFFTs, LinearAlgebra, Printf
function tmin(n; reps=300)
    x = randn(ComplexF64, n)
    p = invoke(AbstractFFTs.plan_fft, Tuple{AbstractArray{ComplexF64,1},Any}, x, 1)
    y = p * x; mul!(y, p, x)
    minimum(@elapsed(mul!(y, p, x)) for _ in 1:reps)
end
haspad = isdefined(FFTA, :bluestein_pad_length)
for n in (4099, 8443, 65537)
    ts = [tmin(n) for _ in 1:3]
    pad = haspad ? FFTA.bluestein_pad_length(n) : nextpow(2, 2n - 1)
    @printf("n=%6d pad=%7d  trials(us): %s   best=%.2f\n", n, pad,
            join((@sprintf("%.1f", t*1e6) for t in ts), ", "), minimum(ts)*1e6)
end
