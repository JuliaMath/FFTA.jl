using FFTA, Test

@testset " forward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(Float64, N, N)
    y = rfft(x)
    y_ref = 0*y
    y_ref[1] = length(x)
    @test y ≈ y_ref
end

@testset "allocations" begin
    X = randn(256, 256)
    Y = rfft(X)
    brfft(Y, 256) # compile
    @test (@test_allocations brfft(Y, 256)) <= 54
end
