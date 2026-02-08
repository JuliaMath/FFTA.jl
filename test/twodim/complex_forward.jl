using FFTA, Test

@testset " forward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(ComplexF64, N, N)
    y = fft(x)
    y_ref = 0*y
    y_ref[1] = length(x)
    @test y ≈ y_ref
    x = randn(N,N)
    @test fft(x) ≈ fft(reshape(x,1,N,N), [2,3])[1,:,:]
    @test fft(x) ≈ fft(reshape(x,1,N,N,1), [2,3])[1,:,:,1]
    @test fft(x) ≈ fft(reshape(x,1,1,N,N,1), [3,4])[1,1,:,:,1]
end

@testset "2D plan, 2D array. Size: $n" for n in 1:64
    @testset "size: ($m, $n)" for m in n:(n + 1)
        X = complex.(randn(m, n), randn(m, n))

        @testset "against naive implementation" begin
            @test naive_2d_fourier_transform(X, FFTA.FFT_FORWARD) ≈ fft(X)
        end

        @testset "allocations" begin
            @test (@test_allocations fft(X)) <= 116
        end
    end
end

@testset "2D plan, ND array. Size: $n" for n in 1:64
    x = randn(ComplexF64, n, n + 1, n + 2)

    @testset "against 1D array with mapslices, r=$r" for r in [[1,2], [1,3], [2,3]]
        @test fft(x, r) == mapslices(fft, x; dims = r)
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch fft(zeros(0, 0))
end
