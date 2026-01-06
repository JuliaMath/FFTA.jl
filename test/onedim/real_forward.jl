using FFTA, Test

@testset verbose = true " forward. N=$N" for N in [8, 11, 15, 16, 27, 100]
    x = ones(Float64, N)
    y = rfft(x)
    y_ref = 0*y
    y_ref[1] = N
    @test y ≈ y_ref atol=1e-12
    @test y == rfft(reshape(x,1,1,N),3)[1,1,:]
    @test y == rfft(reshape(x,N,1),1)[:,1]
end

@testset "More forward tests. Size: $n" for n in 1:64
    x = randn(n)

    @testset "against naive implementation" begin
        y = rfft(x)
        @test naive_1d_fourier_transform(x, FFTA.FFT_FORWARD)[1:(n ÷ 2 + 1)] ≈ y

        @testset "temporarily test real dft separately until used by rfft" begin
            y_dft = similar(y)
            FFTA.fft_dft!(y_dft, x, n, 1, 1, 1, 1, cispi(-2/n))
            @test y ≈ y_dft
        end
    end

    @testset "allocation regression" begin
        @test (@test_allocations rfft(x)) <= 51
    end
end

@testset "error messages" begin
    @test_throws DimensionMismatch rfft(zeros(0))
end