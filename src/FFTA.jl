module FFTA

using DocStringExtensions: TYPEDEF, TYPEDSIGNATURES
using LinearAlgebra: LinearAlgebra
using MuladdMacro: @muladd
using Primes: Primes

export fft, bfft, ifft, rfft, brfft, irfft
export plan_fft, plan_bfft, plan_rfft, plan_brfft

include("callgraph.jl")
include("singleton_twiddle.jl")
include("algos.jl")
include("plan.jl")
include("main.jl")
end
