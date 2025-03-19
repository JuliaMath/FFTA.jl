@enum Direction FFT_FORWARD=-1 FFT_BACKWARD=1
@enum Pow24 POW2=2 POW4=1

abstract type AbstractFFTType end

# Represents a Composite Cooley-Tukey FFT
struct CompositeFFT <: AbstractFFTType end

# Represents a Radix-2 Cooley-Tukey FFT
struct Pow2FFT <: AbstractFFTType end

# Represents a Radix-3 Cooley-Tukey FFT
struct Pow3FFT <: AbstractFFTType end

# Represents a Radix-4 Cooley-Tukey FFT
struct Pow4FFT <: AbstractFFTType end

# Represents an O(N²) DFT
struct DFT <: AbstractFFTType end


@enum FFTEnum compositeFFT dft pow2FFT pow3FFT pow4FFT

"""
$(TYPEDSIGNATURES)
Node of a call graph

# Arguments
`left`: Offset to the left child node
`right`: Offset to the right child node
`type`: Object representing the type of FFT
`sz`: Size of this FFT

"""
struct CallGraphNode{T}
    left::Int
    right::Int
    type::FFTEnum
    sz::Int
    s_in::Int
    s_out::Int
    w::T
end

"""
$(TYPEDSIGNATURES)
Object representing a graph of FFT Calls

# Arguments
`nodes`: Nodes keeping track of the graph
`workspace`: Preallocated Workspace

"""
struct CallGraph{T<:Complex}
    nodes::Vector{CallGraphNode{T}}
    workspace::Vector{Vector{T}}
end

# Get the node in the graph at index i
Base.getindex(g::CallGraph{T}, i::Int) where {T} = g.nodes[i]

"""
$(TYPEDSIGNATURES)
Check if `N` is a power of 2 or 4

"""
function _ispow24(N::Int)
    N < 1 && return nothing
    while N & 0b11 == 0
        N >>= 2
    end
    return N < 3 ? Pow24(N) : nothing
end

"""
$(TYPEDSIGNATURES)
Recursively instantiate a set of `CallGraphNode`s

# Arguments
`nodes`: A vector (which gets expanded) of `CallGraphNode`s
`N`: The size of the FFT
`workspace`: A vector (which gets expanded) of preallocated workspaces
`s_in`: The stride of the input
`s_out`: The stride of the output

"""
function CallGraphNode!(nodes::Vector{CallGraphNode{T}}, N::Int, workspace::Vector{Vector{T}}, s_in::Int, s_out::Int)::Int where {T}
    if N == 0
        throw(DimensionMismatch("array has to be non-empty"))
    end
    w = convert(T, cispi(2/N))
    if iseven(N)
        pow = _ispow24(N)
        if !isnothing(pow)
            push!(workspace, T[])
            push!(nodes, CallGraphNode(0, 0, pow == POW2 ? pow2FFT : pow4FFT, N, s_in, s_out, w))
            return 1
        end
    end
    if N % 3 == 0
        if nextpow(3, N) == N
            push!(workspace, T[])
            push!(nodes, CallGraphNode(0, 0, pow3FFT, N, s_in, s_out, w))
            return 1
        end
    end
    if N == 1 || isprime(N)
        push!(workspace, T[])
        push!(nodes, CallGraphNode(0, 0, dft, N, s_in, s_out, w))
        return 1
    end
    Ns = [first(x) for x in collect(factor(N)) for _ in 1:last(x)]
    if Ns[1] == 2
        N1 = prod(Ns[Ns .== 2])
    elseif Ns[1] == 3
        N1 = prod(Ns[Ns .== 3])
    else
        # Greedy search for closest factor of N to sqrt(N)
        Nsqrt = sqrt(N)
        N_cp = cumprod(Ns[end:-1:1])[end:-1:1]
        N_prox = abs.(N_cp .- Nsqrt)
        _,N1_idx = findmin(N_prox)
        N1 = N_cp[N1_idx]
    end
    N2 = N ÷ N1
    push!(nodes, CallGraphNode(0, 0, dft, N, s_in, s_out, w))
    sz = length(nodes)
    push!(workspace, Vector{T}(undef, N))
    left_len = CallGraphNode!(nodes, N1, workspace, N2, N2*s_out)
    right_len = CallGraphNode!(nodes, N2, workspace, N1*s_in, 1)
    nodes[sz] = CallGraphNode(1, 1 + left_len, compositeFFT, N, s_in, s_out, w)
    return 1 + left_len + right_len
end

"""
$(TYPEDSIGNATURES)
Instantiate a CallGraph from a number `N`

"""
function CallGraph{T}(N::Int) where {T}
    nodes = CallGraphNode{T}[]
    workspace = Vector{Vector{T}}()
    CallGraphNode!(nodes, N, workspace, 1, 1)
    CallGraph(nodes, workspace)
end
