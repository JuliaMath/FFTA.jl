# x86-64 replication of the aarch64 fit behind the Bluestein 3-smooth pad
# cost factor (1.9) and the 2048 floor in bluestein_pad_length.
# Method, matched to the aarch64 run: ComplexF64, planned execution,
# t = min over 200 samples of mul!(y,p,x); c(n) = t / (n*log2 n);
# factor = c(3-smooth) / c(pow2 in the same group).
using FFTA, AbstractFFTs, LinearAlgebra, Printf

function tmin(n; reps=200)
    x = randn(ComplexF64, n)
    p = invoke(AbstractFFTs.plan_fft, Tuple{AbstractArray{ComplexF64,1},Any}, x, 1)
    y = p * x
    mul!(y, p, x)                       # warm
    minimum(@elapsed(mul!(y, p, x)) for _ in 1:reps)
end

# (pow2 reference, candidate 3-smooth / other pad lengths) — same groups as aarch64
groups = [ (256,    [270, 288, 300, 320]),
           (2048,   [2025, 2160, 2250, 2400]),
           (16384,  [8640, 8748, 9000, 9216, 10240]),
           (262144, [131220, 131250, 135000, 138240]) ]

factorize3(n) = (a=0; b=0; m=n; while m%2==0; m÷=2; a+=1; end; while m%3==0; m÷=3; b+=1; end; (a,b,m))

println("## Bluestein pad-length cost calibration — x86-64\n")
println("| group | n | factorization | t (µs) | c = t/(n·log2 n) [ps] | factor vs pow2 |")
println("|:--|--:|:--|--:|--:|--:|")
results = Tuple{Int,Int,Float64}[]
for (p2, cands) in groups
    t2 = tmin(p2); c2 = t2 / (p2 * log2(p2))
    @printf("| %d | %d | 2^%d | %.2f | %.1f | 1.00× (ref) |\n",
            p2, p2, round(Int, log2(p2)), t2*1e6, c2*1e12)
    for n in cands
        t = tmin(n); c = t / (n * log2(n))
        (a,b,r) = factorize3(n)
        fs = r == 1 ? "2^$a·3^$b" : "2^$a·3^$b·$r"
        @printf("| %d | %d | %s | %.2f | %.1f | **%.2f×** |\n", p2, n, fs, t*1e6, c*1e12, c/c2)
        push!(results, (p2, n, c/c2))
    end
end

println("\n### Factors for the lengths the chooser actually considers (3-smooth, i.e. 2^a·3^b only)\n")
pure = [(g,n,f) for (g,n,f) in results if factorize3(n)[3] == 1]
for (g,n,f) in pure
    (a,b,_) = factorize3(n)
    @printf("* %6d (2^%d·3^%d) vs %6d : **%.2f×**\n", n, a, b, g, f)
end
if !isempty(pure)
    fs = [f for (_,_,f) in pure]
    @printf("\nrange %.2f–%.2f, geomean %.2f×\n", minimum(fs), maximum(fs), exp(sum(log.(fs))/length(fs)))
end
