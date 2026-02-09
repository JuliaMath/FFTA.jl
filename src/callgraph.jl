@enum Direction FFT_FORWARD=-1 FFT_BACKWARD=1
@enum Pow24 POW2=2 POW4=1
@enum FFTEnum COMPOSITE_FFT DFT POW3_FFT POW2RADIX4_FFT

"""
$(TYPEDEF)
Node of a call graph

# Arguments
- `left`: Offset to the left child node
- `right`: Offset to the right child node
- `type`: Object representing the type of FFT
- `sz`: Size of this FFT

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
$(TYPEDEF)
Object representing a graph of FFT Calls

# Arguments
- `nodes`: Nodes keeping track of the graph
- `workspace`: Preallocated Workspace

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
    if ispow2(N)
        zero_cnt = trailing_zeros(N)
        return iseven(zero_cnt) ? POW4 : POW2
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Recursively instantiate a set of `CallGraphNode`s

# Arguments
- `nodes`: A vector (which gets expanded) of `CallGraphNode`s
- `N`: The size of the FFT
- `workspace`: A vector (which gets expanded) of preallocated workspaces
- `s_in`: The stride of the input
- `s_out`: The stride of the output

"""
function CallGraphNode!(nodes::Vector{CallGraphNode{T}}, N::Int, workspace::Vector{Vector{T}}, s_in::Int, s_out::Int)::Int where {T}
    if N <= 0
        throw(DimensionMismatch("Array length must be strictly positive"))
    end
    w = cispi(T(2) / N)
    if iseven(N)
        pow = _ispow24(N)
        if !isnothing(pow)
            push!(workspace, T[])
            push!(nodes, CallGraphNode(0, 0, POW2RADIX4_FFT, N, s_in, s_out, w))
            return 1
        end
    elseif N % 3 == 0 && nextpow(3, N) == N
        push!(workspace, T[])
        push!(nodes, CallGraphNode(0, 0, POW3_FFT, N, s_in, s_out, w))
        return 1
    elseif N == 1 || Primes.isprime(N)
        push!(workspace, T[])
        push!(nodes, CallGraphNode(0, 0, DFT, N, s_in, s_out, w))
        return 1
    end
    fzn = Primes.factor(N)
    Nf1, Nf1_cnt = first(fzn)
    if Nf1 == 2 || Nf1 == 3
        N1 = Nf1^Nf1_cnt
    else
        Ns = [first(x) for x in fzn for _ in 1:last(x)]
        # Greedy search for closest factor of N to sqrt(N)
        N_fsqrt =  sqrt(N)
        N_isqrt = isqrt(N)
        N_cp = cumprod(Ns)      # reverse(Ns) another choice
        N1_idx = searchsortedlast(N_cp, N_isqrt)
        N1 = N_cp[N1_idx]       # N1 <= N_isqrt <= N_fsqrt
        if N1_idx != lastindex(N_cp) && (abs(N_cp[N1_idx+1] - N_fsqrt) < (N_fsqrt - N1))
            N1 = N_cp[N1_idx+1] # can be >= N_fsqrt
        end
    end
    N2 = N ÷ N1
    push!(nodes, CallGraphNode(0, 0, DFT, N, s_in, s_out, w))
    sz = length(nodes)
    push!(workspace, Vector{T}(undef, N))
    left_len  = CallGraphNode!(nodes, N1, workspace, N2       , N2 * s_out)
    right_len = CallGraphNode!(nodes, N2, workspace, N1 * s_in,          1)
    nodes[sz] = CallGraphNode(1, 1 + left_len, COMPOSITE_FFT, N, s_in, s_out, w)
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
