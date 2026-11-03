# Does the Bluestein pad chooser make the right call on x86-64?
# Times the actual prime transforms the aarch64 fit was validated on.
using FFTA, AbstractFFTs, LinearAlgebra, Printf
function tmin(n; reps=50)
    x = randn(ComplexF64, n)
    p = invoke(AbstractFFTs.plan_fft, Tuple{AbstractArray{ComplexF64,1},Any}, x, 1)
    y = p * x; mul!(y, p, x)
    minimum(@elapsed(mul!(y, p, x)) for _ in 1:reps)
end
haspad = isdefined(FFTA, :bluestein_pad_length)
println(haspad ? "branch: HAS bluestein_pad_length (E)" : "branch: no bluestein_pad_length (A / baseline)")
for n in (73, 79, 4099, 8443, 65537)
    pad = haspad ? FFTA.bluestein_pad_length(n) : nextpow(2, 2n - 1)
    @printf("n=%6d  pad=%7d  t=%9.2f us\n", n, pad, tmin(n)*1e6)
end
