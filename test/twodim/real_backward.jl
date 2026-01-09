using FFTA, Test

@testset "backward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(Complex{Float64}, N, N)
    y = brfft(x, 2(N-1))
    y_ref = 0*y
    y_ref[1] = N*(2(N-1))
    @test y_ref ≈ y atol=1e-10
end

@testset "2D plan, 2D array. Size: $n" for n in 1:64
    @testset "size: ($m, $n)" for m in n:(n + 1)
        X = randn(m, n)
        # Assuming that rfft works since it is tested independently
        Y = rfft(X)

        @testset "round trip with irfft" begin
            @test X ≈ irfft(Y, m)
        end

        @testset "allocations" begin
            @test (@test_allocations brfft(Y, m)) <= 12050
        end
    end
end

@testset "2D plan, ND array. Size: $n" for n in 1:64
    x = randn(n, n + 1, n + 2)

    @testset "round trip with irfft, r=$r" for r in [(1,2), (1,3), (2,3)]
        @test x ≈ irfft(rfft(x, r), size(x, r[1]), r)
    end

    @testset "against 2D array with mapslices, r=$r" for r in [(1,2), (1,3), (2,3)]
        y = rfft(x, r)
        @test brfft(y, size(x, r[1]), r) == mapslices(t -> brfft(t, size(x, r[1])), y; dims = [r...])
    end
end

@testset "allocations" begin
    X = randn(256, 256)
    Y = rfft(X)
    brfft(Y, 256) # compile
    @test (@test_allocations brfft(Y, 256)) <= 82
end
