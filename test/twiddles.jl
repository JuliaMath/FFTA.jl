using FFTA, Test, LinearAlgebra

# Twiddle tables are built at plan time (see src/callgraph.jl). These tests
# pin down their layout and accuracy against direct evaluation.

@testset "unit_roots: N=$N, $dir, $T" for N in (1, 2, 3, 5, 8, 12, 16, 24, 64, 73, 1000, 4096),
                                      dir in (FFTA.FFT_FORWARD, FFTA.FFT_BACKWARD),
                                      T in (ComplexF64, ComplexF32, Complex{BigFloat})
    W = FFTA.unit_roots(T, N, dir)
    @test length(W) == N
    @test eltype(W) == T
    ref = [FFTA.twiddle(T, dir, k, N) for k in 0:N-1]
    @test maximum(abs.(W .- ref)) <= 2 * eps(real(T))
    if N > 1
        # w^k · w^(N-k) == 1 for the symmetric construction
        @test maximum(abs.(W[2:end] .* reverse(W[2:end]) .- 1)) <= 4 * eps(real(T))
    end
end

@testset "pow2_twiddles layout, N=$N" for N in (2, 4, 8, 16, 32, 64, 256, 1024)
    for dir in (FFTA.FFT_FORWARD, FFTA.FFT_BACKWARD)
        tw = FFTA.pow2_twiddles(ComplexF64, N, dir)
        off = 0
        M = N
        while M > 4
            m = M ÷ 4
            for k in 0:m-1, j in 1:3
                @test tw[off + 3k + j] ≈ FFTA.twiddle(ComplexF64, dir, j * k, M) atol=1e-15
            end
            off += 3m
            M = m
        end
        @test length(tw) == off
    end
end

@testset "pow3_twiddles layout, N=$N" for N in (3, 9, 27, 81, 729)
    for dir in (FFTA.FFT_FORWARD, FFTA.FFT_BACKWARD)
        tw = FFTA.pow3_twiddles(ComplexF64, N, dir)
        off = 0
        M = N
        while M > 3
            m = M ÷ 3
            for k in 0:m-1, j in 1:2
                @test tw[off + 2k + j] ≈ FFTA.twiddle(ComplexF64, dir, j * k, M) atol=1e-15
            end
            off += 2m
            M = m
        end
        @test length(tw) == off
    end
end

@testset "composite_twiddles layout" begin
    for (N1, N2) in ((4, 5), (5, 7), (8, 125), (3, 3))
        N = N1 * N2
        tw = FFTA.composite_twiddles(ComplexF64, N, N1, N2, FFTA.FFT_FORWARD)
        @test length(tw) == (N1 - 1) * (N2 - 1)
        for j1 in 1:N1-1, k2 in 1:N2-1
            @test tw[(j1 - 1) * (N2 - 1) + k2] ≈ FFTA.twiddle(ComplexF64, FFTA.FFT_FORWARD, j1 * k2, N) atol=1e-15
        end
    end
end

@testset "CallGraph carries tables for its direction" begin
    for dir in (FFTA.FFT_FORWARD, FFTA.FFT_BACKWARD)
        g = FFTA.CallGraph{ComplexF64}(1000, FFTA.DEFAULT_BLUESTEIN_CUTOFF, dir)
        @test g.dir === dir
        @test length(g.twiddles) == length(g.nodes) == length(g.blue_index)
        @test all(iszero, g.blue_index)
        @test isempty(g.bluestein)
        for (i, n) in enumerate(g.nodes)
            n.type === FFTA.COMPOSITE_FFT && @test length(g.twiddles[i]) == (g.nodes[i + n.left].sz - 1) * (g.nodes[i + n.right].sz - 1)
            n.type === FFTA.DFT && @test length(g.twiddles[i]) == n.sz
        end
        gb = FFTA.CallGraph{ComplexF64}(2 * 1009, FFTA.DEFAULT_BLUESTEIN_CUTOFF, dir)
        bi = findfirst(n -> n.type === FFTA.BLUESTEIN, gb.nodes)
        @test bi !== nothing
        @test gb.blue_index[bi] == 1
        s = gb.bluestein[1]
        @test s.N == 1009 && s.pad_len == 2048 && length(s.chirp) == 1009 && length(s.chirp_fft) == 2048
    end
    # the opposite direction is rejected
    g = FFTA.CallGraph{ComplexF64}(8, FFTA.DEFAULT_BLUESTEIN_CUTOFF, FFTA.FFT_FORWARD)
    y = zeros(ComplexF64, 8)
    @test_throws ArgumentError FFTA.fft!(y, ones(ComplexF64, 8), 1, 1, FFTA.FFT_BACKWARD, g[1].type, g, 1)
end

@testset "kernels with tables computed on the fly" begin
    # convenience methods used by the tests and for experimentation
    x = randn(ComplexF64, 64)
    y1 = similar(x); y2 = similar(x)
    FFTA.fft_pow2_radix4!(y1, x, 64, 1, 1, 1, 1, FFTA.FFT_FORWARD)
    @test y1 ≈ fft(x)
    FFTA.fft_dft!(y2, x, 64, 1, 1, 1, 1, FFTA.FFT_FORWARD)
    @test y2 ≈ fft(x)
    x3 = randn(ComplexF64, 27); y3 = similar(x3)
    m120 = cispi(2 / 3)
    FFTA.fft_pow3!(y3, x3, 27, 1, 1, 1, 1, m120, FFTA.FFT_FORWARD)
    @test y3 ≈ fft(x3)
    xb = randn(ComplexF64, 101); yb = similar(xb)
    FFTA.fft_bluestein!(yb, xb, FFTA.FFT_FORWARD, 101, 1, 1, 1, 1)
    @test yb ≈ fft(xb)
end

@testset "planned execution does not allocate, n=$n" for n in (5, 64, 73, 101, 720, 1000, 1009, 4096, 65537)
    x = randn(ComplexF64, n)
    p = plan_fft(x)
    y = p * x
    mul!(y, p, x)
    @test (@test_allocations mul!(y, p, x)) == 0
end
