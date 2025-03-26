module FFTA

using Primes, DocStringExtensions, Reexport, MuladdMacro, LinearAlgebra
@reexport using AbstractFFTs

import AbstractFFTs: Plan

include("callgraph.jl")
include("algos.jl")
include("plan.jl")


end
