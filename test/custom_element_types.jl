using Test, FFTA

@testset "element type: $T" for T in (Float16, BigFloat)
    x = randn(2*3*4*5)
    Tx = T.(x)

    @testset "AbstractFFTs believes that single and double precision is everything." begin
        # Ref https://github.com/JuliaMath/FFTA.jl/issues/77
        @test_broken fft(Tx) isa Vector{Complex{T}}
    end

    # Complex
    cTx = complex(Tx)
    new_cTx = ifft(fft(cTx))
    @test typeof(new_cTx) == typeof(cTx)
    @test cTx ≈ new_cTx

    # Real
    new_Tx = irfft(rfft(Tx), length(Tx))
    @test typeof(new_Tx) == typeof(Tx)
    @test Tx ≈ new_Tx
end
