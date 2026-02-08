using Test, FFTA
using LinearAlgebra: LinearAlgebra

@testset "Only 1D and 2D FFTs" begin
    xr = zeros(2, 2)
    xc = complex(xr)
    @test_throws ArgumentError("only supports 1D and 2D FFTs") plan_fft(xc, 1:3)
    @test_throws ArgumentError("only supports 1D and 2D FFTs") plan_bfft(xc, 1:3)
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
end

@testset "mismatch between input and output arrays" begin
    @testset "1D plan 1D array" begin
        x1 = complex(randn(3))
        y1 = similar(x1, length(x1) + 1)

        @test_throws DimensionMismatch LinearAlgebra.mul!(y1, plan_fft(x1), x1)
    end

    @testset "2D array" begin
        x2 = complex.(randn(3, 3), randn(3, 3))
        y2 = similar(x2, size(x2, 1) + 1, size(x2, 2) + 1)

        @testset "1D plan, region=$(region)" for region in [1, 2]
            @test_throws DimensionMismatch LinearAlgebra.mul!(y2, plan_fft(x2, region), x2)
        end

        @testset "2D plan" begin
            @test_throws DimensionMismatch LinearAlgebra.mul!(y2, plan_fft(x2), x2)
        end
    end
end
