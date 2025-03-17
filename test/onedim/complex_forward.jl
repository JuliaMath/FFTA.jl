using FFTA, Test
test_nums = [8, 11, 15, 16, 27, 100]
@testset verbose = true " forward" begin
    for N in test_nums
        x = ones(ComplexF64, N)
        y = fft(x)
        y_ref = 0*y
        y_ref[1] = N
        @test y ≈ y_ref atol=1e-12
    end
end

@testset verbose = true "against naive implementation. Size: $n" for n in 1:64
    x = complex.(randn(n), randn(n))
    @test naive_1d_fourier_transform(x, FFTA.FFT_FORWARD) ≈ fft(x)
end

@testset "error messages" begin
    @test_throws DimensionMismatch fft(zeros(0))
end
