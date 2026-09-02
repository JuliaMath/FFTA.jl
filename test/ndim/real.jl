using FFTA, Test, LinearAlgebra

# Real transforms over three or more dimensions: real along the first region
# dimension, complex along the others.

@testset "rfft/irfft, size $sz, region $r" for sz in ((4, 5, 6), (5, 4, 6), (8, 3, 5, 2), (3, 3, 3, 3)),
                                                r in (1:3, (1, 2, 3), (2, 3), (1, 3), [1, 2, 4], 1:4)
    maximum(r) <= length(sz) || continue
    x = randn(sz)
    y = rfft(x, r)
    ref = fft(complex(x), r)
    d1 = first(r)
    @test size(y) == ntuple(i -> i == d1 ? sz[i] ÷ 2 + 1 : sz[i], length(sz))
    @test y ≈ selectdim(ref, d1, 1:sz[d1] ÷ 2 + 1)
    @test irfft(y, sz[d1], r) ≈ x
    @test brfft(y, sz[d1], r) ≈ prod(sz[i] for i in r) * x
    p = plan_rfft(x, r)
    y2 = similar(y)
    @test mul!(y2, p, x) == y
    @test inv(p) * y ≈ x
    @test p \ y ≈ x
    pb = plan_brfft(y, sz[d1], r)
    xb = similar(x)
    @test mul!(xb, pb, y) ≈ prod(sz[i] for i in r) * x
    x32 = randn(Float32, sz)
    @test irfft(rfft(x32, r), sz[d1], r) ≈ x32
end

@testset "3D real transform, multiple workers" begin
    x = randn(16, 24, 20)
    y1 = plan_rfft(x, 1:3; num_threads = 1) * x
    @test plan_rfft(x, 1:3; num_threads = 3) * x == y1
    @test y1 ≈ selectdim(fft(complex(x)), 1, 1:9)
end
