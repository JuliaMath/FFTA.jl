using FFTA, Test

@testset "Basic ND checks" begin
    for sz in ((3, 5, 7), (4, 14, 9), (103, 5, 13), (26, 33, 35, 4), ntuple(i -> 3, 5))
        x = ones(sz)
        @test fft(x) ≈ setindex!(zeros(sz), prod(sz), 1)
    end

    y = zeros((3, 3, 3))
    y[2, 2, 2] = 1
    w1 = -0.5 - sqrt(3)im / 2
    w2 = conj(w1)
    y_ref = reshape(ComplexF64[
            1 w1 w2;
            w1 w2 1;
            w2 1 w1
            ;;;
            w1 w2 1;
            w2 1 w1;
            1 w1 w2
            ;;;
            w2 1 w1;
            1 w1 w2;
            w1 w2 1
        ], 3, 3, 3)
    @test isapprox(fft(y), y_ref)
end
