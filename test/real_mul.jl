using FFTA, Test, LinearAlgebra

# `mul!` on real plans (forward: real -> complex, backward: complex -> real),
# including N-d arrays with a `dims` argument, 2D plans, views, and the
# absence of allocations for 1D plans.

@testset "1D real plans, mul!, n=$n" for n in (1, 2, 3, 4, 8, 9, 15, 16, 31, 64, 100, 101, 1000)
    x = randn(n)
    p = plan_rfft(x)
    y = p * x
    y2 = similar(y)
    @test mul!(y2, p, x) === y2
    @test y2 == y
    @test y ≈ naive_1d_fourier_transform(x, FFTA.FFT_FORWARD)[1:(n ÷ 2 + 1)]
    pb = plan_brfft(y, n)
    xb = similar(x)
    @test mul!(xb, pb, y) === xb
    @test xb ≈ n * x
    @test xb == pb * y
    if n < FFTA.DEFAULT_BLUESTEIN_CUTOFF   # Bluestein still allocates scratch per call
        @test (@test_allocations mul!(y2, p, x)) == 0
        @test (@test_allocations mul!(xb, pb, y)) == 0
    end
    # views as input and output
    X = randn(n, 3)
    Y = zeros(ComplexF64, n ÷ 2 + 1, 3)
    mul!(view(Y, :, 2), p, view(X, :, 2))
    @test Y[:, 2] ≈ rfft(X[:, 2])
    Xb = zeros(n, 3)
    mul!(view(Xb, :, 3), pb, view(Y, :, 2))
    @test Xb[:, 3] ≈ n * X[:, 2]
    # wrong sizes
    @test_throws DimensionMismatch mul!(similar(y, n ÷ 2 + 2), p, x)
    @test_throws DimensionMismatch mul!(similar(x, n + 1), pb, y)
    @test_throws ArgumentError mul!(similar(x), p, y)
    @test_throws ArgumentError mul!(similar(y), pb, x)
end

@testset "1D real plans on N-d arrays, mul!, size $sz" for sz in ((7, 4), (8, 5), (6, 7, 8), (5, 8, 9))
    x = randn(sz)
    for d in 1:length(sz)
        n = size(x, d)
        p = plan_rfft(x, d)
        y = p * x
        y2 = similar(y)
        @test mul!(y2, p, x) == y
        @test y ≈ mapslices(rfft, x; dims = d)
        @test (@test_allocations mul!(y2, p, x)) == 0
        pb = plan_brfft(y, n, d)
        xb = similar(x)
        @test mul!(xb, pb, y) ≈ n * x
        @test (@test_allocations mul!(xb, pb, y)) == 0
        @test_throws DimensionMismatch mul!(similar(x, ntuple(i -> i == d ? n + 1 : sz[i], length(sz))), pb, y)
    end
end

@testset "2D real plans, mul!, size $sz" for sz in ((8, 6), (9, 6), (8, 7), (9, 7), (1, 1), (2, 3), (64, 64), (100, 101))
    x = randn(sz)
    p = plan_rfft(x)
    y = p * x
    @test y ≈ naive_2d_fourier_transform(x, FFTA.FFT_FORWARD)[1:(sz[1] ÷ 2 + 1), :]
    y2 = similar(y)
    @test mul!(y2, p, x) == y
    pb = plan_brfft(y, sz[1])
    xb = similar(x)
    @test mul!(xb, pb, y) ≈ prod(sz) * x
    @test xb == pb * y
    @test_throws DimensionMismatch mul!(similar(y, sz[1] ÷ 2 + 1, sz[2] + 1), p, x)
    @test_throws DimensionMismatch mul!(xb, pb, similar(y, sz[1] ÷ 2 + 1, sz[2] + 1))
end

@testset "2D real plans on 3D arrays, mul!, region $r" for r in ((1, 2), (1, 3), (2, 3), [1, 2], 2:3)
    x = randn(6, 7, 8)
    p = plan_rfft(x, r)
    y = p * x
    @test y ≈ mapslices(rfft, x; dims = r)
    y2 = similar(y)
    @test mul!(y2, p, x) == y
    n1 = size(x, first(r))
    pb = plan_brfft(y, n1, r)
    xb = similar(x)
    @test mul!(xb, pb, y) ≈ prod(size(x, i) for i in r) * x
end

@testset "Float32 real plans" begin
    x = randn(Float32, 48)
    y = rfft(x)
    @test eltype(y) == ComplexF32
    @test y ≈ rfft(Float64.(x))
    @test irfft(y, 48) ≈ x
    X = randn(Float32, 12, 10)
    @test irfft(rfft(X), 12) ≈ X
end
