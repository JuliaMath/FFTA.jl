using FFTA, Test
test_nums = [8, 11, 15, 16, 27, 100]
@testset "backward" begin
    for N in test_nums
        x = ones(ComplexF64, N, N)
        y = bfft(x)
        y_ref = 0*y
        y_ref[1] = length(x)
        @test y ≈ y_ref
    end
end

@testset verbose = true "against naive implementation" for n in 1:64
    @testset "size: ($m, $n)" for m in n:(n + 1)
        X = complex.(randn(m, n), randn(m, n))
        Y = similar(X)
        @test naive_2d_fourier_transform(X, FFTA.FFT_BACKWARD) ≈ bfft(X)
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch bfft(zeros(ComplexF64, 0, 0))
end
