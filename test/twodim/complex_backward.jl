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

@testset "2D plan, ND array. Size: $n" for n in 1:64
    x = randn(ComplexF64, n, n + 1, n + 2)

    @testset "against 1D array with mapslices, r=$r" for r in [[1,2], [1,3], [2,3]]
        @test bfft(x, r) == mapslices(bfft, x; dims = r)
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch bfft(zeros(ComplexF64, 0, 0))
end
