using FFTA, Test

@testset verbose = true " forward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(ComplexF64, N)
    y = fft(x)
    y_ref = 0*y
    y_ref[1] = N
    @test y ≈ y_ref atol=1e-12
    @test y == fft(reshape(x,1,1,N),3)[1,1,:]
    @test y == fft(reshape(x,N,1), 1)[:,1]
end

@testset "1D plan, 1D array. Size: $n" for n in 1:64
    x = complex.(randn(n), randn(n))

    @testset "against naive implementation" begin
        @test naive_1d_fourier_transform(x, FFTA.FFT_FORWARD) ≈ fft(x)
    end

    @testset "allocation regression" begin
        @test (@test_allocations fft(x)) <= 47
    end
end

@testset "1D plan, ND array. Size: $n" for n in 1:64
    x = complex.(randn(n, n + 1, n + 2), randn(n, n + 1, n + 2))

    @testset "against 1D array with mapslices, r=$r" for r in 1:3
        @test fft(x, r) == mapslices(fft, x; dims = r)
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch fft(zeros(0))
end
