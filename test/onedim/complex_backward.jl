using FFTA, Test

@testset "backward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(ComplexF64, N)
    y = bfft(x)
    y_ref = 0*y
    y_ref[1] = N
    @test y ≈ y_ref atol=1e-12
end

@testset "More backward tests. Size: $n" for n in 1:64
    x = complex.(randn(n), randn(n))

    @testset "against naive implementation" begin
        @test naive_1d_fourier_transform(x, FFTA.FFT_BACKWARD) ≈ bfft(x)
    end

    @testset "allocation" begin
        @test (@allocations bfft(x)) <= 44
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch bfft(zeros(0))
end
