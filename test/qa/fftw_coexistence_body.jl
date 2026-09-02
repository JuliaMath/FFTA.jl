using FFTA, FFTW, Test

# When FFTW.jl is loaded alongside FFTA, AbstractFFTs' design is that FFTW's
# `StridedArray` methods take over. Make sure FFTA's methods do not make the
# generic entry points ambiguous instead (this used to turn `rfft(x)` into a
# `MethodError` when `region` was annotated with `RegionTypes`).
@testset "coexistence with FFTW.jl" begin
    x = randn(16)
    y = randn(ComplexF64, 16)
    @test rfft(x) ≈ FFTW.rfft(x)
    @test irfft(rfft(x), 16) ≈ x
    @test fft(y) ≈ FFTW.fft(y)
    @test plan_rfft(x) isa FFTW.FFTWPlan
    @test plan_brfft(rfft(x), 16) isa FFTW.FFTWPlan
    @test plan_fft(y) isa FFTW.FFTWPlan
    @test plan_bfft(y) isa FFTW.FFTWPlan
    X = randn(8, 6)
    @test rfft(X, 2) ≈ FFTW.rfft(X, 2)
    @test irfft(rfft(X, 1:2), 8, 1:2) ≈ X
    # FFTA's own methods are still reachable for non-strided arrays
    @test rfft(view(x, [1:16;])) ≈ FFTW.rfft(x)
    @test plan_rfft(view(x, [1:16;])) isa FFTA.FFTAPlan
    @test isempty(filter(Test.detect_ambiguities(FFTA, FFTW)) do (m1, m2)
        m1.name in (:plan_fft, :plan_bfft, :plan_rfft, :plan_brfft)
    end)
end
