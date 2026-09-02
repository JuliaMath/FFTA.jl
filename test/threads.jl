using FFTA, Test, LinearAlgebra

# Plans own one worker (call-graph workspace + pencil buffers) per thread they
# may use; pencils of multidimensional and batched transforms are distributed
# over them. Results must not depend on the worker count. The suite normally
# runs with one Julia thread, in which case the parallel tasks simply run one
# after the other; the chunking logic is exercised either way.

@testset "num_threads keyword" begin
    x = randn(ComplexF64, 16)
    @test length(plan_fft(x; num_threads = 1).workers) == 1
    @test length(plan_fft(x; num_threads = 3).workers) == 3
    @test length(plan_fft(x).workers) == Threads.nthreads()
    @test length(plan_rfft(randn(16); num_threads = 2).workers) == 2
    @test_throws ArgumentError plan_fft(x; num_threads = 0)
    p = plan_fft(x; num_threads = 2)
    @test length(inv(p).p.workers) == 2                # inverse keeps the worker count
    @test p.workers[1].callgraph[1] === p.callgraph[1]
    @test p.workers[2].callgraph[1] !== p.callgraph[1]
    @test p.workers[2].callgraph[1].nodes === p.callgraph[1].nodes   # nodes are shared
end

# sizes above and below THREAD_THRESHOLD, contiguous and strided pencils,
# composite lengths (which use the call-graph workspace)
@testset "results independent of worker count, size $sz, region $r" for sz in ((64, 64), (513, 40), (30, 700, 4), (300, 300), (8, 8, 8)),
                                                                          r in (1, 2, (1, 2))
    X = randn(ComplexF64, sz)
    p1 = plan_fft(X, r; num_threads = 1)
    Y1 = p1 * X
    @test Y1 ≈ mapslices(fft, X; dims = r)
    for nt in (2, 3, 5)
        p = plan_fft(X, r; num_threads = nt)
        @test p * X == Y1
        Y = similar(Y1)
        @test mul!(Y, p, X) == Y1
        @test inv(p) * Y1 ≈ X
        pb = plan_bfft(X, r; num_threads = nt)
        @test pb * X == plan_bfft(X, r; num_threads = 1) * X
    end
    Xr = randn(Float64, sz)
    Yr1 = plan_rfft(Xr, r; num_threads = 1) * Xr
    @test Yr1 ≈ mapslices(rfft, Xr; dims = r)
    for nt in (2, 3)
        pr = plan_rfft(Xr, r; num_threads = nt)
        @test pr * Xr == Yr1
        @test inv(pr) * Yr1 ≈ Xr
        pbr = plan_brfft(Yr1, size(Xr, first(r)), r; num_threads = nt)
        @test pbr * Yr1 == plan_brfft(Yr1, size(Xr, first(r)), r; num_threads = 1) * Yr1
    end
end

@testset "single-worker execution does not allocate" begin
    X = randn(ComplexF64, 64, 64)
    for r in (1, 2, (1, 2))
        p = plan_fft(X, r; num_threads = 1)
        Y = p * X
        mul!(Y, p, X)
        @test (@test_allocations mul!(Y, p, X)) == 0
    end
    Xr = randn(64, 64)
    for r in (1, 2)
        p = plan_rfft(Xr, r; num_threads = 1)
        Y = p * Xr
        mul!(Y, p, Xr)
        @test (@test_allocations mul!(Y, p, Xr)) == 0
    end
end
