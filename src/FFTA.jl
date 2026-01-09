module FFTA

using DocStringExtensions: TYPEDSIGNATURES
using LinearAlgebra: LinearAlgebra
using MuladdMacro: @muladd
using Primes: Primes
using Reexport: @reexport

export fft, bfft, ifft, rfft, brfft, irfft
export plan_fft, plan_bfft, plan_rfft, plan_brfft

include("callgraph.jl")
include("algos.jl")
include("plan.jl")
include("main.jl")

end
