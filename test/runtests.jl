using Test, Random, FFTA

macro test_allocations(args)
    @static if Base.VERSION >= v"1.9"
        ex = Expr(:macrocall, Symbol("@allocations"), __source__, args)
        return esc(ex)
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
            a = FFTA.direction_sign(d) * (2 * u * v) / n
            s = muladd(x[v + 1], cispi(a), s)
        end
        y[u + 1] = s
    end

    return y
end

function naive_2d_fourier_transform(X::Matrix, d::FFTA.Direction)
    rows, cols = size(X)
    Y = zeros(Complex{Float64}, rows, cols)

    for u in 0:(rows - 1), v in 0:(cols - 1)
        s = 0.0 + 0.0im
        for x in 0:(rows - 1), y in 0:(cols - 1)
            a = FFTA.direction_sign(d) * 2 * ((u * x) / rows + (v * y) / cols)
            s = muladd(X[x + 1, y + 1], cispi(a), s)
        end
        Y[u + 1, v + 1] = s
    end

    return Y
end

Random.seed!(1)
@testset verbose = true "FFTA" begin
    @testset verbose = true "QA" begin
        include("qa/aqua.jl")
        include("qa/explicit_imports.jl")
    end
    @testset verbose = true "Argument checking" begin
        include("argument_checking.jl")
    end
    @testset verbose = true "1D" begin
        @testset verbose = true "Complex" begin
            @testset verbose = false "Forward" begin
                include("onedim/complex_forward.jl")
            end
            @testset verbose = false "Backward" begin
                include("onedim/complex_backward.jl")
            end
        end
        @testset verbose = true "Real" begin
            @testset verbose = false "Forward" begin
                include("onedim/real_forward.jl")
            end
            @testset verbose = false "Backward" begin
                include("onedim/real_backward.jl")
            end
        end
    end
    @testset verbose = true "2D" begin
        @testset verbose = true "Complex" begin
            @testset verbose = false "Forward" begin
                include("twodim/complex_forward.jl")
            end
            @testset verbose = false "Backward" begin
                include("twodim/complex_backward.jl")
            end
        end
        @testset verbose = true "Real" begin
            @testset verbose = false "Forward" begin
                include("twodim/real_forward.jl")
            end
            @testset verbose = false "Backward" begin
                include("twodim/real_backward.jl")
            end
        end
    end
    @testset verbose = true "Custom element types" begin
        include("custom_element_types.jl")
    end
end