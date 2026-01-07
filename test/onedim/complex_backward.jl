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
    # Assuming that fft works since it is tested independently
    y = fft(x)

    @testset "round tripping with ifft" begin
        @test ifft(y) ≈ x
    end

    @testset "allocation regression" begin
        @test (@test_allocations bfft(y)) <= 47
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch bfft(complex.(zeros(0)))
end
