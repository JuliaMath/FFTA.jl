using FFTA, Test
test_nums = [8, 11, 15, 16, 27, 100]
@testset verbose = true " forward" begin
    for N in test_nums
        x = ones(Float64, N)
        y = rfft(x)
        y_ref = 0*y
        y_ref[1] = N
        @test y ≈ y_ref atol=1e-12
    end
end

@testset verbose = true "against naive implementation. Size: $n" for n in 1:64
    x = randn(n)
    y = rfft(x)
    @test naive_1d_fourier_transform(x, FFTA.FFT_FORWARD)[1:(n ÷ 2 + 1)] ≈ y

    @testset "temporarily test real dft separately until used by rfft" begin
        y_dft = similar(y)
        FFTA.fft_dft!(y_dft, x, n, 1, 1, 1, 1, cispi(-2/n))
        @test y ≈ y_dft
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch rfft(zeros(0))
end