# x86-64 replication of the DFT-vs-Bluestein crossover measurement that set
# BLUESTEIN_CUTOFF = 47 on aarch64.
# BLUESTEIN_CUTOFF=1000 forces the O(n^2) DFT leaf; =2 forces Bluestein.
using FFTA, AbstractFFTs, LinearAlgebra, Printf

function tmin(n, cutoff; reps=200)
    x = randn(ComplexF64, n)
    p = invoke(AbstractFFTs.plan_fft, Tuple{AbstractArray{ComplexF64,1},Any}, x, 1;
               BLUESTEIN_CUTOFF = cutoff)
    y = p * x
    mul!(y, p, x)
    minimum(@elapsed(mul!(y, p, x)) for _ in 1:reps)
end

primes = [13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
println("## DFT vs Bluestein crossover — x86-64 (ComplexF64, planned exec, min of 200)\n")
println("| n | DFT leaf (µs) | Bluestein (µs) | Bluestein/DFT | cheaper |")
println("|--:|--:|--:|--:|:--|")
cross = Ref{Union{Nothing,Int}}(nothing)
for n in primes
    td = tmin(n, 1000); tb = tmin(n, 2)
    win = tb < td ? "**Bluestein**" : "DFT"
    if cross[] === nothing && tb < td; cross[] = n; end
    @printf("| %d | %.2f | %.2f | %.2f× | %s |\n", n, td*1e6, tb*1e6, tb/td, win)
end
println()
println(cross[] === nothing ? "DFT wins at every n tested — cutoff should be above $(last(primes))" :
        "**x86-64 crossover: Bluestein first wins at n = $(cross[])** (aarch64 value: 47)")
