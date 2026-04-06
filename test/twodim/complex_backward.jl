using FFTA, Test

@testset "backward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(ComplexF64, N, N)
    y = bfft(x)
    y_ref = 0*y
    y_ref[1] = length(x)
    @test y ≈ y_ref
end

@testset "2D plan, 2D array. Size: $n" for n in 1:64
    @testset "size: ($m, $n)" for m in n:(n + 1)
        x = randn(ComplexF64, n)
        # Assuming that fft works since it is tested independently
        y = fft(x)

        @testset "round tripping with ifft" begin
            @test ifft(y) ≈ x
        end

        @testset "allocations" begin
            @test (@test_allocations bfft(y)) <= 116
        end
    end
end

@testset "$(N)D plan, $(N+1)D array" for N in 2:3
    rg = N == 2 ? (1:64) : (1:16)
    dims_lst = [[1,2], [1,3], [2,3]]
    if N == 3
        foreach(v -> push!(v, 4), dims_lst)
    end
    @testset "against $(N)D arrays with mapslices, r=$r" for r in dims_lst
        for n in rg
            x = randn(ComplexF64, ntuple(i -> n + (i - 1), N + 1))

            t = Tuple(r)    # test tuple region argument
            @test bfft(x, t) == bfft(x, r) == mapslices(bfft, x; dims = r)
        end
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch bfft(zeros(ComplexF64, 0, 0))
end
