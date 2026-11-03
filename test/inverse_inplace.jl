using FFTA, Test, LinearAlgebra

# Inverse plans (`inv`, `\`, `ldiv!`, `plan_ifft`, `plan_irfft`) and in-place
# plans (`plan_fft!`, `plan_bfft!`, `plan_ifft!`, `fft!`, `bfft!`, `ifft!`).

_inplace_allocs(w, p) = (mul!(w, p, w); @test_allocations mul!(w, p, w))

@testset "1D, n=$n, $T" for n in (1, 2, 8, 9, 15, 64, 100, 101, 1000), T in (ComplexF64, ComplexF32)
    x = randn(T, n)
    p = plan_fft(x)
    y = p * x
    @test inv(p) * y ≈ x
    @test p \ y ≈ x
    z = similar(x)
    @test ldiv!(z, p, y) === z
    @test z ≈ x
    @test inv(p) === inv(p)            # cached in the plan
    @test inv(inv(p)) * x ≈ y
    @test plan_ifft(x) * y ≈ x
    pb = plan_bfft(x)
    @test inv(pb) * (pb * x) ≈ x

    xr = randn(real(T), n)
    pr = plan_rfft(xr)
    yr = pr * xr
    @test inv(pr) * yr ≈ xr
    @test pr \ yr ≈ xr
    zr = similar(xr)
    @test ldiv!(zr, pr, yr) === zr
    @test zr ≈ xr
    pbr = plan_brfft(yr, n)
    @test inv(pbr) * (pbr * yr) ≈ yr
    @test plan_irfft(yr, n) * yr ≈ xr

    p! = plan_fft!(x)
    w = copy(x)
    @test (p! * w) === w
    @test w ≈ y
    w = copy(x)
    @test mul!(w, p!, w) === w
    @test w ≈ y
    v = similar(x)
    @test mul!(v, p!, x) ≈ y          # out of place use of an in-place plan
    pb! = plan_bfft!(x)
    w = copy(y)
    pb! * w
    @test w ≈ n * x
    w = copy(y)
    @test (inv(p!) * w) === w
    @test w ≈ x
    @test inv(p!) === inv(p!)
    w = copy(x)
    fft!(w)
    @test w ≈ y
    ifft!(w)
    @test w ≈ x
    bfft!(w)
    @test w ≈ bfft(x)
    @test plan_ifft!(x) * copy(y) ≈ x
    if n < FFTA.DEFAULT_BLUESTEIN_CUTOFF
        @test _inplace_allocs(copy(x), p!) == 0
    end
    @test_throws DimensionMismatch mul!(zeros(T, n + 1), p!, zeros(T, n + 1))
end

@testset "N-d complex, region $r" for r in (1, 2, 3, (1, 2), (2, 3), (1, 2, 3), 1:3, [1, 3])
    X = randn(ComplexF64, 6, 7, 8)
    p = plan_fft(X, r)
    Y = p * X
    @test inv(p) * Y ≈ X
    @test p \ Y ≈ X
    Z = similar(X)
    ldiv!(Z, p, Y)
    @test Z ≈ X
    p! = plan_fft!(X, r)
    W = copy(X)
    @test (p! * W) === W
    @test W ≈ Y
    W = copy(X)
    fft!(W, r)
    @test W ≈ Y
    ifft!(W, r)
    @test W ≈ X
    bfft!(W, r)
    @test W ≈ bfft(X, r)
    if r isa Int   # plans over several dimensions still allocate their pencil buffers per call
        @test _inplace_allocs(copy(X), p!) == 0
    end
end

@testset "N-d real, region $r" for r in (1, 2, (1, 2), (2, 3), [1, 2])
    Xr = randn(8, 6, 5)
    p = plan_rfft(Xr, r)
    Y = p * Xr
    @test inv(p) * Y ≈ Xr
    @test p \ Y ≈ Xr
    Z = similar(Xr)
    ldiv!(Z, p, Y)
    @test Z ≈ Xr
    pb = plan_brfft(Y, size(Xr, first(r)), r)
    @test inv(pb) * (pb * Y) ≈ Y
end

@testset "in-place plan, plan and inverse reused (DSP.jl pattern)" begin
    x = randn(ComplexF64, 32)
    p! = plan_fft!(x)
    ip! = inv(p!)
    @test ip!.p isa FFTA.FFTAPlan_inplace
    buf = copy(x)
    p! * buf
    ip! * buf
    @test buf ≈ x
end
