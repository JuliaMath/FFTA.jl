using FFTA, Test, Random, LinearAlgebra

# Regression test for issue #118. Singleton's trigonometric recurrence
# (see src/singleton_twiddle.jl) brings the twiddle error from O(N)
# down to roughly O(log N · √log N) + an O(√N) twiddle term. Bounds
# below are set with ~2-3× headroom over what the current
# implementation produces (max over 5 seeds on aarch64 NEON); they
# still fail comfortably against the pre-fix naive `w *= step`
# recurrence, which ballooned past ~4000 ULP at N = 16384.

# (N, max eps ratio) across the power-of-2 ladder. Covers both even
# powers (= powers of 4, recursion bottoms at N = 4) and odd powers
# (recursion bottoms at N = 2), which hit different base cases in
# `fft_pow2_radix4!`.
const POWERS_OF_2 = (
    (1 << 4,   3.0),   # 16    = 4^2
    (1 << 5,   3.0),   # 32
    (1 << 6,   3.0),   # 64    = 4^3
    (1 << 7,   3.0),   # 128
    (1 << 8,   4.0),   # 256   = 4^4
    (1 << 9,   4.0),   # 512
    (1 << 10,  8.0),   # 1024  = 4^5
    (1 << 11,  8.0),   # 2048
    (1 << 12, 10.0),   # 4096  = 4^6
    (1 << 13, 12.0),   # 8192
    (1 << 14, 14.0),   # 16384 = 4^7
    (1 << 15, 20.0),   # 32768
    (1 << 16, 22.0),   # 65536 = 4^8
    (1 << 17, 28.0),   # 131072
    (1 << 18, 28.0),   # 262144 = 4^9
)

const POWERS_OF_3 = (
    (3^1,    3.0),
    (3^2,    3.0),
    (3^3,    3.0),
    (3^4,    4.0),
    (3^5,    5.0),
    (3^6,    7.0),
    (3^7,   10.0),
    (3^8,   13.0),
    (3^9,   18.0),
)

function _worst_relerr(N::Int)
    worst = 0.0
    for seed in 1:5
        rng = Xoshiro(seed)
        x64 = randn(rng, ComplexF64, N)
        x32 = ComplexF32.(x64)
        y32 = fft(x32)
        y_ref = ComplexF32.(fft(x64))
        relerr = norm(y32 .- y_ref) / norm(y_ref)
        worst = max(worst, relerr / eps(Float32))
    end
    return worst
end

@testset "twiddle accuracy (issue #118)" begin
    @testset "powers of 2" begin
        @testset "N = $N" for (N, bound) in POWERS_OF_2
            @test _worst_relerr(N) ≤ bound
        end
    end
    @testset "powers of 3" begin
        @testset "N = $N" for (N, bound) in POWERS_OF_3
            @test _worst_relerr(N) ≤ bound
        end
    end
end
