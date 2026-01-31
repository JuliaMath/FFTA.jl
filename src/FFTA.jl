module FFTA

using AbstractFFTs: AbstractFFTs
using DocStringExtensions: TYPEDEF, TYPEDSIGNATURES
using LinearAlgebra: LinearAlgebra
using MuladdMacro: @muladd
using Primes: Primes
using Reexport: @reexport

@reexport using AbstractFFTs

include("callgraph.jl")
include("algos.jl")
include("plan.jl")
end
