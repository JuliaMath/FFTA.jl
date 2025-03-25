using FFTA, Test

@testset "backward. N=$N" for N in [8]
    x = ones(Float64, N, N)
    y = brfft(x, 2(N-1))
    y_ref = 0*y
    y_ref[1] = N*(2(N-1))
    @test y_ref ≈ y atol=1e-12
end

@testset "allocations" begin
    X = randn(256, 256)
    rfft(X) # compile
    @test (@allocations rfft(X)) <= 51
end
