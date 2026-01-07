using FFTA, Test, LinearAlgebra

@testset "backward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(Complex{Float64}, N)
    y = brfft(x, 2*(N-1))
    y_ref = 0*y
    y_ref[1] = 2*(N-1)
    if !isapprox(y_ref, y, atol=1e-12)
        println(norm(y_ref - y))
    end
    @test y_ref ≈ y atol=1e-12
    @test y isa Vector{<:Real}
end

@testset "More backward tests. Size: $n" for n in 1:64
    x = randn(n)
    # Assuming that rfft works since it is tested separately
    y = rfft(x)

    @testset "round tripping with irfft" begin
        @test irfft(y, n) ≈ x
    end

    @testset "allocation regression" begin
        @test (@test_allocations brfft(y, n)) <= 55
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch brfft(zeros(ComplexF64, 0), 0)
end
