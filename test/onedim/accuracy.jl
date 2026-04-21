using FFTA, Test, Random, LinearAlgebra

# Regression test for issue #118. Singleton's trigonometric recurrence
# (see src/singleton_twiddle.jl) brings the twiddle error from O(N)
# down to roughly O(log N · √log N) without the per-step trig cost
# that direct `cispi` evaluation would impose. These bounds are set
# with ~3× headroom over what the current implementation produces;
# they still fail against the pre-fix naive `w *= step` recurrence.

Random.seed!(42)

const ACCURACY_CASES = (
    (256,    5.0),
    (1024,  15.0),
    (4096,  20.0),
    (16384, 30.0),
    (81,     5.0),
    (243,    5.0),
    (729,   10.0),
)

@testset "twiddle accuracy (issue #118)" begin
    @testset "N = $N" for (N, bound) in ACCURACY_CASES
        worst = 0.0
        for seed in 1:5
            Random.seed!(seed)
            x64 = randn(ComplexF64, N)
            x32 = ComplexF32.(x64)
            y32 = fft(x32)
            y_ref = ComplexF32.(fft(x64))
            relerr = norm(y32 .- y_ref) / norm(y_ref)
            worst = max(worst, relerr / eps(Float32))
        end
        @test worst ≤ bound
    end
end
