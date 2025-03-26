using Test, Random, FFTA

macro test_allocations(args)
    if Base.VERSION >= v"1.9"
        :(@allocations($(esc(args))))
    else
        :(0)
    end
end

function naive_1d_fourier_transform(x::Vector, d::FFTA.Direction)
    n = length(x)
    y = zeros(Complex{Float64}, n)

    for u in 0:(n - 1)
        s = 0.0 + 0.0im
        for v in 0:(n - 1)
            a = FFTA.direction_sign(d) * 2π * u * v / n
            s += x[v + 1] * exp(im * a)
        end
        y[u + 1] = s
    end

    return y
end

function naive_2d_fourier_transform(X::Matrix, d::FFTA.Direction)
    rows, cols = size(X)
    Y = zeros(Complex{Float64}, rows, cols)

    for u in 0:(rows - 1)
        for v in 0:(cols - 1)
            s = 0.0 + 0.0im
            for x in 0:(rows - 1)
                for y in 0:(cols - 1)
                    a = FFTA.direction_sign(d) * 2π * (u * x / rows + v * y / cols)
                    s += X[x + 1, y + 1] * exp(im * a)
                end
            end
            Y[u + 1, v + 1] = s
        end
    end

    return Y
end

Random.seed!(1)
@testset verbose = true "FFTA" begin
    @testset verbose = true "1D" begin
        @testset verbose = false "Complex" begin
            include("onedim/complex_forward.jl")
            include("onedim/complex_backward.jl")
            x = rand(ComplexF64, 100)
            y = fft(x)
            x2 = bfft(y)/length(x)
            @test x ≈ x2 atol=1e-12
        end
        @testset verbose = false "Real" begin
            include("onedim/real_forward.jl")
            include("onedim/real_backward.jl")
            x = rand(Float64, 100)
            y = fft(x)
            x2 = bfft(y)/length(x)
            @test x ≈ x2 atol=1e-12
        end
    end
    @testset verbose = false "2D" begin
        @testset verbose = true "Complex" begin
            include("twodim/complex_forward.jl")
            include("twodim/complex_backward.jl")
            x = rand(ComplexF64, 100, 100)
            y = fft(x)
            x2 = bfft(y)/length(x)
            @test x ≈ x2
        end
        @testset verbose = true "Real" begin
            include("twodim/real_forward.jl")
            include("twodim/real_backward.jl")
            x = rand(Float64, 100, 100)
            y = fft(x)
            x2 = bfft(y)/length(x)
            @test x ≈ x2
        end
    end
end