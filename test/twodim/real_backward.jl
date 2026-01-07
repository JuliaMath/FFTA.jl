using FFTA, Test

@testset "backward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(Float64, N, N)
    y = brfft(x, 2(N-1))
    y_ref = 0*y
    y_ref[1] = N*(2(N-1))
    @test y_ref ≈ y atol=1e-10
end

@testset "More backward tests" for n in 1:64
    @testset "size: ($m, $n)" for m in n:(n + 1)
        X = randn(m, n)
        # Assuming that rfft works since it is tested independently
        Y = rfft(X)

        @testset "round trip" begin
            @test X ≈ irfft(Y, m)
        end

        @testset "allocations" begin
            @test (@test_allocations bfft(X)) <= 12050
        end
    end
end

@testset "allocations" begin
    X = randn(256, 256)
    Y = rfft(X)
    brfft(Y, 256) # compile
    @test (@test_allocations brfft(Y, 256)) <= 68
end
