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

@testset "1D plan, 1D array. Size: $n" for n in 1:64
    x = randn(n)
    # Assuming that rfft works since it is tested separately
    y = rfft(x)

    @testset "round tripping with irfft" begin
        @test irfft(y, n) ≈ x
    end

    @testset "allocation regression" begin
        brfft(y, n) # warmup
        @test (@test_allocations brfft(y, n)) <= 55
    end
end

@testset "1D plan, ND array. Size: $n" for n in 1:64
    x = randn(n, n + 1, n + 2)

    @testset "round tripping with irfft, r=$r" for r in 1:3
        @test irfft(rfft(x, r), size(x,r), r) ≈ x
    end

    @testset "against 1D array with mapslices, r=$r" for r in 1:3
        y = rfft(x, r)
        @test brfft(y, size(x, r), r) == mapslices(t -> brfft(t, size(x, r)), y; dims = r)
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch brfft(zeros(ComplexF64, 0), 0)
end
