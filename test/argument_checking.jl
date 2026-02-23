using Test, FFTA
using LinearAlgebra: LinearAlgebra

@testset "Only 1D and 2D FFTs" begin
    xr = zeros(2, 2, 2)
    xc = complex(xr)

    @test_throws ArgumentError("only supports 1D and 2D FFTs") plan_rfft(xr, 1:3)
    @test_throws ArgumentError("only supports 1D and 2D FFTs") plan_brfft(xc, 2, 1:3)
end

@testset "mismatch between plan and array" begin
    @testset "1D plan 1D array" begin
        xr1 = randn(3)
        yr1 = rfft(xr1)
        xc1 = complex.(xr1)

        xr1p = [xr1; 1]
        xc1p = [xc1; 1]
        yr1p = [yr1; 1]

        @test_throws DimensionMismatch plan_fft(xc1) * xc1p
        @test_throws DimensionMismatch plan_bfft(xc1) * xc1p
        @test_throws DimensionMismatch plan_rfft(xr1) * xr1p
        @test_throws DimensionMismatch plan_brfft(yr1, length(xr1)) * yr1p
    end

    @testset "2D array" begin
        xr2 = randn(3, 3)
        xc2 = complex.(xr2)

        xr2p = [[xr2; ones(1, size(xr2, 2))] ones(size(xr2, 1) + 1, 1)]
        xc2p = [[xc2; ones(1, size(xr2, 2))] ones(size(xc2, 1) + 1, 1)]

        @testset "1D plan, region=$(region)" for region in 1:2
            yr2 = rfft(xr2, region)

            yr2p = if region == 1
                [yr2; ones(1, size(yr2, 2))]
            else
                [yr2 ones(size(yr2, 1), 1)]
            end

            @test_throws DimensionMismatch plan_fft(xc2, region) * xc2p
            @test_throws DimensionMismatch plan_bfft(xc2, region) * xc2p
            @test_throws DimensionMismatch plan_rfft(xr2, region) * xr2p
            @test_throws DimensionMismatch plan_brfft(yr2, size(xr2, region), region) * yr2p
        end

        @testset "2D plan" begin
            yr2 = rfft(xr2)

            yr2p = [yr2; ones(1, 3)]

            @test_throws DimensionMismatch plan_fft(xc2) * xc2p
            @test_throws DimensionMismatch plan_bfft(xc2) * xc2p
            @test_throws DimensionMismatch plan_rfft(xr2) * xr2p
            @test_throws DimensionMismatch plan_brfft(yr2, size(xr2, 1)) * yr2p
        end
    end
    @testset "3D array" begin
        xc3 = randn(ComplexF64, 3, 3, 3)
        yc3 = randn(ComplexF64, 5, 5, 5)
        pxc3 = plan_fft(xc3)
        @test_throws DimensionMismatch pxc3 * yc3
        invalid_p = plan_fft(randn(ComplexF64, ntuple(i -> 3, 5)), 3:5)
        xc4 = randn(ComplexF64, (1, ntuple(i -> 5, 3)...))

        ### plan region out of bounds

        # all same dims
        @test_throws DimensionMismatch("Plan region is outside array dimensions.") invalid_p * xc3
        # dim(p) < dim(out) = dim(in)
        @test_throws DimensionMismatch("Plan region is outside array dimensions.") LinearAlgebra.mul!(xc4, invalid_p, xc4)
    end
end

@testset "mismatch between input and output arrays" begin
    @testset "1D plan 1D array" begin
        x1 = randn(ComplexF64, 3)
        y1 = similar(x1, length(x1) + 1)

        @test_throws DimensionMismatch LinearAlgebra.mul!(y1, plan_fft(x1), x1)
    end

    @testset "$(N)D array" for N in 2:3
        xN = randn(ComplexF64, ntuple(i -> 3, N))
        yN = similar(xN, size(xN) .+ 1)

        @testset "1D plan, region=$(region)" for region in 1:N
            @test_throws DimensionMismatch LinearAlgebra.mul!(yN, plan_fft(xN, region), xN)
        end

        @testset "$(N)D plan" begin
            @test_throws DimensionMismatch LinearAlgebra.mul!(yN, plan_fft(xN), xN)
        end
    end
end
