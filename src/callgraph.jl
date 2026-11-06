@enum Direction FFT_FORWARD=-1 FFT_BACKWARD=1
@enum Pow24 POW2 POW4
@enum FFTEnum COMPOSITE_FFT DFT POW3_FFT POW2RADIX4_FFT BLUESTEIN

@inline function direction_sign(d::Direction)
    Int(d)
end

"""
$(TYPEDEF)
Node of a call graph

# Arguments
- `left`: Offset to the left child node
- `right`: Offset to the right child node
- `type`: Object representing the type of FFT
- `sz`: Size of this FFT
- `s_in`: The stride of the input
- `s_out`: The stride of the output

"""
struct CallGraphNode
    left::Int
    right::Int
    type::FFTEnum
    sz::Int
    s_in::Int
    s_out::Int
end

"""
$(TYPEDEF)
Precomputed data and scratch space for Bluestein's algorithm on a length-`N`
transform in direction `dir` (see `fft_bluestein!`).

# Fields
- `N`: Length of the transform
- `pad_len`: Length ≥ 2N-1 of the padded convolution (see `bluestein_pad_length`)
- `chirp`: `b_n = exp(∓iπ n²/N)` for `n = 0..N-1`
- `chirp_fft`: forward FFT of the zero-padded, periodised chirp, divided by `pad_len`
- `a`, `tmp`: work arrays of length `pad_len`
- `graph`: forward `CallGraph` for the length-`pad_len` transforms
"""
struct BluesteinScratch{T<:Complex,G}
    N::Int
    pad_len::Int
    chirp::Vector{T}
    chirp_fft::Vector{T}
    a::Vector{T}
    tmp::Vector{T}
    graph::G
end

"""
$(TYPEDEF)
Object representing a graph of FFT Calls

# Arguments
- `nodes`: Nodes keeping track of the graph
- `workspace`: Preallocated Workspace
- `twiddles`: Precomputed twiddle factors of each node (see `node_twiddles`)
- `bluestein`: Precomputed data for each `BLUESTEIN` node, indexed through `blue_index`
- `blue_index`: Index into `bluestein` for each node (`0` for other node types)
- `dir`: Direction the twiddle factors were computed for
- `BLUESTEIN_CUTOFF`: Minimum prime that will be FFTed with the
    Bluestein algorithm, below which the O(N^2) DFT is used.

"""
struct CallGraph{T<:Complex}
    nodes::Vector{CallGraphNode}
    workspace::Vector{Vector{T}}
    twiddles::Vector{Vector{T}}
    bluestein::Vector{BluesteinScratch{T,CallGraph{T}}}
    blue_index::Vector{Int}
    dir::Direction
    BLUESTEIN_CUTOFF::Int
end

# Primes below this use the O(N²) DFT with a twiddle table; at and above it
# Bluestein's algorithm is used. Measured crossovers (ComplexF64, planned
# execution): x86-64 AVX2 (Core Ultra 7 165H) n ≈ 23 — Bluestein is 1.7–2.2×
# faster than the DFT leaf for n = 41–47; aarch64 NEON (Neoverse-N1) n ≈ 29
# for the 64-point pad, with the DFT leaf up to 1.45× faster again at n = 37–43
# where the pad grows to 128, and Bluestein ahead from 47. 29 never loses
# much on either.
const DEFAULT_BLUESTEIN_CUTOFF = 29

# Relative per-element cost of a 3-smooth transform length compared to a
# power of two, above which `bluestein_pad_length` keeps the power of two.
# See its docstring for the measurements behind the value.
const BLUESTEIN_SMOOTH_FACTOR = 2.1

# Get the node in the graph at index i
Base.getindex(g::CallGraph{T}, i::Int) where {T} = g.nodes[i]

"""
$(TYPEDSIGNATURES)
Recursively instantiate a set of `CallGraphNode`s

# Arguments
- `nodes`: A vector (which gets expanded) of `CallGraphNode`s
- `N`: The size of the FFT
- `workspace`: A vector (which gets expanded) of preallocated workspaces
- `BLUESTEIN_CUTOFF`: Minimum prime that will be FFTed with the
    Bluestein algorithm, below which the O(N^2) DFT is used.
- `s_in`: The stride of the input
- `s_out`: The stride of the output

"""
function CallGraphNode!(
    nodes::Vector{CallGraphNode},
    N::Int,
    workspace::Vector{Vector{T}},
    BLUESTEIN_CUTOFF::Int,
    s_in::Int, s_out::Int
)::Int where {T}
    if N <= 0
        throw(DimensionMismatch("Array length must be strictly positive"))
    end
    if iseven(N) && ispow2(N)
        push!(workspace, T[])
        push!(nodes, CallGraphNode(0, 0, POW2RADIX4_FFT, N, s_in, s_out))
        return 1
    elseif N % 3 == 0 && nextpow(3, N) == N
        push!(workspace, T[])
        push!(nodes, CallGraphNode(0, 0, POW3_FFT, N, s_in, s_out))
        return 1
    elseif N == 1 || Primes.isprime(N)
        push!(workspace, T[])
        # use Bluestein's algorithm for big primes
        LEAF_ALG = N < BLUESTEIN_CUTOFF ? DFT : BLUESTEIN
        push!(nodes, CallGraphNode(0, 0, LEAF_ALG, N, s_in, s_out))
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
        if N1_idx != lastindex(N_cp) && (N_cp[N1_idx+1] - N_fsqrt < (N_fsqrt - N1))
            N1 = N_cp[N1_idx+1] # can be >= N_fsqrt
        end
    end
    N2 = N ÷ N1
    push!(nodes, CallGraphNode(0, 0, DFT, N, s_in, s_out))
    sz = length(nodes)
    push!(workspace, Vector{T}(undef, N))
    left_len  = CallGraphNode!(nodes, N1, workspace, BLUESTEIN_CUTOFF, N2       , N2 * s_out)
    right_len = CallGraphNode!(nodes, N2, workspace, BLUESTEIN_CUTOFF, N1 * s_in,          1)
    nodes[sz] = CallGraphNode(1, 1 + left_len, COMPOSITE_FFT, N, s_in, s_out)
    return 1 + left_len + right_len
end

# ---------------------------------------------------------------------------
# Twiddle factors
#
# All twiddle factors are computed once at plan time with `cispi`, which is
# accurate to a rounding error, and stored per node in a layout chosen so that
# the kernels read them sequentially. `k` is reduced modulo `N` before the
# conversion to floating point so that large exponents lose no accuracy.
# ---------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)
The twiddle factor `exp(dir · 2πi · k / N)` as an element of type `T`.
"""
@inline function twiddle(::Type{T}, dir::Direction, k::Integer, N::Integer) where {T<:Complex}
    R = real(T)
    T(cispi(2 * R(direction_sign(dir) * mod(k, N)) / R(N)))
end

"""
$(TYPEDSIGNATURES)
The `N` unit roots `W[k + 1] = exp(dir · 2πi k / N)`, `k = 0..N-1`. When `8 | N`
only the first octant is evaluated with `sincospi` and the rest follows from
symmetry, which makes planning large power-of-two transforms cheap.
"""
function unit_roots(::Type{T}, N::Int, dir::Direction) where {T<:Complex}
    W = Vector{T}(undef, N)
    if N % 8 != 0
        for k in 0:N-1
            W[k + 1] = twiddle(T, dir, k, N)
        end
        return W
    end
    R = real(T)
    sgn = R(direction_sign(dir))
    q = N ÷ 4
    for k in 0:N÷8
        s, c = sincospi(2 * R(k) / R(N))
        s *= sgn
        # with s = sgn·sin θ, c = cos θ and the imaginary part carrying `sgn`:
        W[k + 1]          = T(c, s)                 # θ
        W[q - k + 1]      = T(sgn * s, sgn * c)     # π/2 - θ
        W[q + k + 1]      = T(-sgn * s, sgn * c)    # π/2 + θ
        W[2q - k + 1]     = T(-c, s)                # π - θ
        W[2q + k + 1]     = T(-c, -s)               # π + θ
        W[3q - k + 1]     = T(-sgn * s, -sgn * c)   # 3π/2 - θ
        W[3q + k + 1]     = T(sgn * s, -sgn * c)    # 3π/2 + θ
        k > 0 && (W[N - k + 1] = T(c, -s))          # 2π - θ
    end
    return W
end

"""
$(TYPEDSIGNATURES)
Twiddle table for the O(N²) DFT: `W[k + 1] = exp(dir · 2πi k / N)`, `k = 0..N-1`.
"""
dft_twiddles(::Type{T}, N::Int, dir::Direction) where {T} = unit_roots(T, N, dir)

"""
$(TYPEDSIGNATURES)
Twiddle table for the Cooley-Tukey step of a composite `N = N1 · N2` transform:
the factor for output row `j1` and inner index `k2` is stored at
`(j1 - 1) * (N2 - 1) + k2` for `j1 = 1..N1-1`, `k2 = 1..N2-1`.
"""
function composite_twiddles(::Type{T}, N::Int, N1::Int, N2::Int, dir::Direction) where {T}
    W = unit_roots(T, N, dir)
    tw = Vector{T}(undef, (N1 - 1) * (N2 - 1))
    i = 1
    for j1 in 1:N1-1
        idx = 0
        for k2 in 1:N2-1
            idx += j1
            idx >= N && (idx -= N)
            tw[i] = W[idx + 1]
            i += 1
        end
    end
    return tw
end

"""
$(TYPEDSIGNATURES)
Twiddle table for the radix-4 power-of-two kernel. For every recursion level of
size `M = N, N/4, N/16, …` (down to, but excluding, the 4- and 2-point base
cases) the table holds the triplets `(w^k, w^2k, w^3k)`, `w = exp(dir · 2πi/M)`,
for `k = 0..M/4-1`, one level after the other. A level of size `M` therefore
occupies `3M/4` entries and the next level starts `3M/4` entries later.
"""
function pow2_twiddles(::Type{T}, N::Int, dir::Direction) where {T}
    N > 4 || return T[]
    W = unit_roots(T, N, dir)
    tw = Vector{T}(undef, (N - 2) )   # 3N/4 + 3N/16 + ... < N
    i = 1
    M = N
    while M > 4
        m = M ÷ 4
        s = N ÷ M         # w_M^k = w_N^(s k)
        for k in 0:m-1
            tw[i]     = W[s * k + 1]
            tw[i + 1] = W[2 * s * k + 1]
            tw[i + 2] = W[3 * s * k + 1]
            i += 3
        end
        M = m
    end
    resize!(tw, i - 1)
    return tw
end

"""
$(TYPEDSIGNATURES)
Twiddle table for the radix-3 kernel, laid out like `pow2_twiddles`: for every
level of size `M = N, N/3, …` (excluding the 3-point base case) the pairs
`(w^k, w^2k)`, `w = exp(dir · 2πi/M)`, for `k = 0..M/3-1`.
"""
function pow3_twiddles(::Type{T}, N::Int, dir::Direction) where {T}
    N > 3 || return T[]
    W = unit_roots(T, N, dir)
    tw = Vector{T}(undef, N)   # 2N/3 + 2N/9 + ... < N
    i = 1
    M = N
    while M > 3
        m = M ÷ 3
        s = N ÷ M
        for k in 0:m-1
            tw[i]     = W[s * k + 1]
            tw[i + 1] = W[2 * s * k + 1]
            i += 2
        end
        M = m
    end
    resize!(tw, i - 1)
    return tw
end

"""
$(TYPEDSIGNATURES)
Twiddle table of the node at index `idx` of `nodes`, see `dft_twiddles`,
`composite_twiddles`, `pow2_twiddles` and `pow3_twiddles`. `BLUESTEIN` nodes
keep their data in a `BluesteinScratch` instead and get an empty table.
"""
function node_twiddles(::Type{T}, nodes::Vector{CallGraphNode}, idx::Int, dir::Direction) where {T}
    node = nodes[idx]
    N = node.sz
    if node.type === COMPOSITE_FFT
        N1 = nodes[idx + node.left].sz
        N2 = nodes[idx + node.right].sz
        return composite_twiddles(T, N, N1, N2, dir)
    elseif node.type === DFT
        return dft_twiddles(T, N, dir)
    elseif node.type === POW2RADIX4_FFT
        return pow2_twiddles(T, N, dir)
    elseif node.type === POW3_FFT
        return pow3_twiddles(T, N, dir)
    else
        return T[]
    end
end

"""
$(TYPEDSIGNATURES)
Length of the padded convolution in Bluestein's algorithm for a length-`N`
transform: the smallest power of two ≥ 2N-1, unless a 3-smooth length
`2^a 3^b ≥ 2N-1` is enough smaller to be cheaper.

How much cheaper it has to be is machine dependent. Measured per `n log n`
against a power of two (ComplexF64, planned execution), a length with
factors of 3 costs 1.3–2.3× on aarch64 NEON (Neoverse-N1) but 2.0–3.1× on
x86-64 AVX2 (Core Ultra 7 165H); with a threshold of 1.9 the 3-smooth pad
was a 1.58× win at N = 8443 on aarch64 and a 1.19× loss on x86-64. The
threshold (`factor`) is therefore set to 2.1, at which the chooser never
loses on either machine — which, since the smallest admissible length is
always more than half the power of two, means it currently never picks a
3-smooth length: the mechanism is kept (and the constant should be
lowered) for when the composite/radix-3 path becomes competitive with the
radix-4 kernel. Below 2048 the fixed overhead of the composite step
dominates and the power of two is always used.
"""
function bluestein_pad_length(N::Int; factor::Real = BLUESTEIN_SMOOTH_FACTOR)
    m = 2N - 1
    p = nextpow(2, m)
    best = p
    best_cost = p * log2(p)
    pow3 = 3
    while pow3 < p
        c = pow3 * nextpow(2, cld(m, pow3))   # smallest 2^a·3^b with this power of 3
        cost = factor * c * log2(c)
        if c >= max(m, 2048) && cost < best_cost
            best, best_cost = c, cost
        end
        pow3 *= 3
    end
    return best
end

"""
$(TYPEDSIGNATURES)
Precompute the chirp, its transform, the work arrays and the call graph of the
padded transform used by `fft_bluestein!` for a length-`N` transform in
direction `d`.
"""
function BluesteinScratch{T}(N::Int, d::Direction) where {T<:Complex}
    pad_len = bluestein_pad_length(N)
    R = real(T)

    # chirp b_n = exp(sgn · iπ n²/N); n² is tracked modulo 2N so that the
    # argument of `cispi` stays small.
    chirp = Vector{T}(undef, N)
    sgn = -direction_sign(d)
    p = 0   # n^2 mod 2N, kept in (-N, N]
    for i in 1:N
        chirp[i] = T(cispi(R(sgn * p) / R(N)))
        p += (2i - 1)   # prevents overflow unless N is absolutely massive
        p > N && (p -= 2N)
    end

    # the padded length is 3-smooth, so its graph has no Bluestein node
    graph = CallGraph{T}(pad_len, 2, FFT_FORWARD)

    # periodised, zero-padded chirp and its forward transform, scaled by
    # 1/pad_len so that the inverse transform in `fft_bluestein!` needs no
    # further normalisation
    a = zeros(T, pad_len)
    copyto!(a, 1, chirp, 1, N)
    for j in 0:N-2
        a[pad_len - j] = chirp[2 + j]
    end
    chirp_fft = Vector{T}(undef, pad_len)
    fft!(chirp_fft, a, 1, 1, FFT_FORWARD, graph[1].type, graph, 1)
    chirp_fft ./= R(pad_len)

    return BluesteinScratch{T,CallGraph{T}}(N, pad_len, chirp, chirp_fft, a, Vector{T}(undef, pad_len), graph)
end

"""
$(TYPEDSIGNATURES)
Instantiate a CallGraph from a number `N`, with twiddle factors for direction `dir`

"""
function CallGraph{T}(N::Int, BLUESTEIN_CUTOFF::Int, dir::Direction=FFT_FORWARD) where {T}
    nodes = CallGraphNode[]
    workspace = Vector{Vector{T}}()
    CallGraphNode!(nodes, N, workspace, BLUESTEIN_CUTOFF, 1, 1)
    twiddles = [node_twiddles(T, nodes, idx, dir) for idx in eachindex(nodes)]
    bluestein = BluesteinScratch{T,CallGraph{T}}[]
    blue_index = zeros(Int, length(nodes))
    for (idx, node) in enumerate(nodes)
        if node.type === BLUESTEIN
            push!(bluestein, BluesteinScratch{T}(node.sz, dir))
            blue_index[idx] = length(bluestein)
        end
    end
    CallGraph(nodes, workspace, twiddles, bluestein, blue_index, dir, BLUESTEIN_CUTOFF)
end
