using FFTA, Test

@testset " forward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(Float64, N, N)
    y = rfft(x)
    y_ref = 0*y
    y_ref[1] = length(x)
    @test y ≈ y_ref
    x = randn(N,N)
    @test rfft(x) ≈ rfft(reshape(x,1,N,N), [2,3])[1,:,:]
    @test rfft(x) ≈ rfft(reshape(x,1,N,N,1), [2,3])[1,:,:,1]
    @test rfft(x) ≈ rfft(reshape(x,1,1,N,N,1), [3,4])[1,1,:,:,1]
    @test size(rfft(x)) == (N÷2+1, N)
end

@testset "allocations" begin
    X = randn(256, 256)
    rfft(X) # compile
    @test (@test_allocations rfft(X)) <= 61
end
