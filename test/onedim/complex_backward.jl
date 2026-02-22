using FFTA, Test

# add N=103 to test Bluestein (for large prime lengths)
@testset "backward. N=$N" for N in [8, 11, 15, 16, 27, 100, 103]
    x = ones(ComplexF64, N)
    y = bfft(x)
    y_ref = zero(y)
    y_ref[1] = N
    # atol = N == 103 ? 0.0 : 1e-12   # Bluestein has larger error, can use rtol
    @test isapprox(y, y_ref; atol=1e-12)
end

@testset "1D plan, 1D array. Size: $n" for n in 1:64
    x = randn(ComplexF64, n)
    # Assuming that fft works since it is tested independently
    y = fft(x)

    @testset "round tripping with ifft" begin
        @test ifft(y) ≈ x
    end

    @testset "allocation regression" begin
        @test (@test_allocations bfft(y)) <= 47
    end
end

@testset "1D plan, ND array. Size: $n" for n in 1:64
    x = randn(ComplexF64, n, n + 1, n + 2)

    @testset "against 1D array with mapslices, r=$r" for r in 1:3
        @test bfft(x, r) == mapslices(bfft, x; dims = r)
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch bfft(complex.(zeros(0)))
end
