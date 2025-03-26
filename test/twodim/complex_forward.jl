using FFTA, Test

@testset " forward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(ComplexF64, N, N)
    y = fft(x)
    y_ref = 0*y
    y_ref[1] = length(x)
    @test y ≈ y_ref
end

@testset "More forward tests" for n in 1:64
    @testset "size: ($m, $n)" for m in n:(n + 1)
        X = complex.(randn(m, n), randn(m, n))

        @testset "against naive implementation" begin
            @test naive_2d_fourier_transform(X, FFTA.FFT_FORWARD) ≈ fft(X)
        end

        @testset "allocations" begin
            @test (@allocations fft(X)) <= 111
        end
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch fft(zeros(0, 0))
end
