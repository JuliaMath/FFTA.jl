using FFTA, Test, LinearAlgebra
test_nums = [8, 11, 15, 16, 27, 100]
@testset "backward" begin
    for N in test_nums
        x = ones(Float64, N)
        y = brfft(x, 2*(N-1))
        y_ref = 0*y
        y_ref[1] = 2*(N-1)
        if !isapprox(y_ref, y, atol=1e-12)
            println(norm(y_ref - y))
        end
        @test y_ref ≈ y atol=1e-12
    end
end

@testset verbose = true "against naive implementation. Size: $n" for n in 1:64
    x = complex.(randn(n ÷ 2 + 1), randn(n ÷ 2 + 1))
    x[begin] = real(x[begin])
    if iseven(n)
        x[end] = real(x[end])
        xe = [x; conj.(reverse(x[begin + 1:end - 1]))]
    else
        xe = [x; conj.(reverse(x[begin + 1:end]))]
    end
    @test naive_1d_fourier_transform(xe, FFTA.FFT_BACKWARD) ≈ brfft(x, n)
end

@testset "error messages" begin
    @test_throws DimensionMismatch brfft(zeros(ComplexF64, 0), 0)
end
