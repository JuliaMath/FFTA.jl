# Plans

abstract type FFTAPlan{T,N} <: AbstractFFTs.Plan{T} end

const RegionTypes{N} = Union{Int,AbstractVector{Int},NTuple{N,Int}}

# The plan types are mutable so that `AbstractFFTs.inv` can cache the inverse
# plan in the initially undefined `pinv` field (see `plan_inv`), as FFTW.jl does.

"""
$(TYPEDEF)
Per-thread state of a plan: call graphs (sharing the nodes of the plan's call
graphs but with their own workspace) and pencil buffers. A plan owns one
`Worker` per thread it may use (see the `num_threads` keyword of the `plan_*`
functions); pencils of a multidimensional or batched transform are
distributed over them.
"""
struct Worker{T,N,S<:Real}
    callgraph::NTuple{N,CallGraph{T}}
    obuf::Vector{T}   # transform output of one pencil before it is copied back
    buf::Vector{T}    # real<->complex packing scratch, see `_re_buflen` (real plans only)
    rbuf::Vector{S}   # contiguous copy of a strided real pencil, see `_re_pencil_loop!` (real plans only)
    cbuf::Vector{T}   # contiguous copy of a strided complex pencil (real plans only)
end

# A call graph that shares the (immutable) nodes but has its own workspace.
_clone_workspace(g::CallGraph{T}) where {T} =
    CallGraph{T}(g.nodes, [Vector{T}(undef, length(w)) for w in g.workspace], g.BLUESTEIN_CUTOFF)

# `buflen` is the packing scratch length and `rlen` the real length of a real
# plan (0 for complex plans, which need neither).
function _workers(cg::Tuple{CallGraph{T},Vararg{CallGraph{T}}}, ::Val{N}, num_threads::Int, buflen::Int, rlen::Int = 0) where {T,N}
    num_threads >= 1 || throw(ArgumentError("num_threads must be at least 1"))
    S = real(T)
    obuflen = maximum(g -> first(g.nodes).sz, cg)
    clen = rlen == 0 ? 0 : rlen ÷ 2 + 1
    mk(graphs) = Worker{T,N,S}(graphs, Vector{T}(undef, obuflen), Vector{T}(undef, buflen), Vector{S}(undef, rlen), Vector{T}(undef, clen))
    workers = [mk(cg)]
    for _ in 2:num_threads
        push!(workers, mk(map(_clone_workspace, cg)))
    end
    return workers
end

# Transforms with fewer elements than this are not split over threads.
const THREAD_THRESHOLD = 1 << 15

mutable struct FFTAPlan_cx{T,N,R<:RegionTypes{N},S<:Real} <: FFTAPlan{T,N}
    const callgraph::NTuple{N,CallGraph{T}}
    const region::R
    const dir::Direction
    const workers::Vector{Worker{T,N,S}}
    pinv::AbstractFFTs.ScaledPlan
    FFTAPlan_cx{T,N,R}(cg::NTuple{N,CallGraph{T}}, r::R, dir::Direction, workers::Vector{Worker{T,N,S}}) where {T,N,R<:RegionTypes{N},S} =
        new{T,N,R,S}(cg, r, dir, workers)
end
function FFTAPlan_cx{T,N}(
    cg::NTuple{N,CallGraph{T}}, r::R,
    dir::Direction; num_threads::Int=Threads.nthreads()
) where {T,N,R<:RegionTypes{N}}
    FFTAPlan_cx{T,N,R}(cg, r, dir, _workers(cg, Val(N), num_threads, 0))
end

mutable struct FFTAPlan_re{T,N,R<:RegionTypes{N},S<:Real} <: FFTAPlan{T,N}
    const callgraph::NTuple{N,CallGraph{T}}
    const region::R
    const dir::Direction
    const flen::Int
    const workers::Vector{Worker{T,N,S}}
    pinv::AbstractFFTs.ScaledPlan
    FFTAPlan_re{T,N,R}(cg::NTuple{N,CallGraph{T}}, r::R, dir::Direction, flen::Int, workers::Vector{Worker{T,N,S}}) where {T,N,R<:RegionTypes{N},S} =
        new{T,N,R,S}(cg, r, dir, flen, workers)
end
function FFTAPlan_re{T,N}(
    cg::NTuple{N,CallGraph{T}}, r::R,
    dir::Direction, flen::Int; num_threads::Int=Threads.nthreads()
) where {T,N,R<:RegionTypes{N}}
    FFTAPlan_re{T,N,R}(cg, r, dir, flen, _workers(cg, Val(N), num_threads, _re_buflen(flen, dir), flen))
end

_num_threads(p::Union{FFTAPlan_cx,FFTAPlan_re}) = length(p.workers)

# Scratch length needed by the real-transform pencil kernels for a plan of
# real length `n` (see `_rfft_pencil!` / `_brfft_pencil!`).
function _re_buflen(n::Int, dir::Direction)
    if iseven(n)
        dir === FFT_FORWARD ? n >> 1 : n
    else
        dir === FFT_FORWARD ? n : 2n
    end
end

function Base.size(p::FFTAPlan{<:Any,N}, i::Int) where N
    if i < 1
        throw(DomainError(i, "No non-positive dimensions"))
    elseif i > N
        1
    elseif p isa FFTAPlan_re && i == 1
        p.flen
    else
        first(p.callgraph[i].nodes).sz
    end
end
Base.size(p::FFTAPlan{<:Any,N}) where N = ntuple(Base.Fix1(size, p), Val{N}())

function _sort(region::T)::T where {N,T<:NTuple{N,Int}}
    @static if VERSION >= v"1.12"
        sort(region)
    else
        if N == 2
            minmax(region[1], region[2])
        elseif N == 3
            t1, t2, t3 = region
            t1, t2 = minmax(t1, t2)
            t2, t3 = minmax(t2, t3)
            t1, t2 = minmax(t1, t2)
            (t1, t2, t3)
        else
            NTuple{N}(sort!(collect(region)))
        end
    end
end

_sort(region::T) where T<:RegionTypes = issorted(region) ? copy(region) : sort(region)

AbstractFFTs.plan_fft(x::AbstractArray{T,N}, region; kwargs...) where {T<:Complex,N} =
    _plan_fft(x, region, FFT_FORWARD; kwargs...)

AbstractFFTs.plan_bfft(x::AbstractArray{T,N}, region; kwargs...) where {T<:Complex,N} =
    _plan_fft(x, region, FFT_BACKWARD; kwargs...)

function _plan_fft(
    x::AbstractArray{T,N},
    region::RegionTypes,
    dir::Direction;
    BLUESTEIN_CUTOFF=DEFAULT_BLUESTEIN_CUTOFF, num_threads::Int=Threads.nthreads(), _kwargs...
) where {T<:Complex,N}
    M = length(region)
    if M == 1
        R1 = Int(region[1])
        g = CallGraph{T}(size(x, R1), BLUESTEIN_CUTOFF)
        return FFTAPlan_cx{T,1}((g,), R1, dir; num_threads)
    elseif M == 2
        R2 = _sort(region)
        g1 = CallGraph{T}(size(x, R2[1]), BLUESTEIN_CUTOFF)
        g2 = CallGraph{T}(size(x, R2[2]), BLUESTEIN_CUTOFF)
        return FFTAPlan_cx{T,2}((g1, g2), R2, dir; num_threads)
    else
        RM = _sort(region)
        return FFTAPlan_cx{T,M}(
            ntuple(i -> CallGraph{T}(size(x, RM[i]), BLUESTEIN_CUTOFF), Val(M)),
            RM, dir; num_threads
        )
    end
end

# The `AbstractFFTs` entry points deliberately leave `region` unannotated: an
# annotation such as `region::RegionTypes` makes these methods ambiguous with
# FFTW.jl's `plan_rfft(::StridedArray{Float64}, region)` when both packages
# are loaded (neither is more specific in every argument), which turns
# `rfft(x)` into a `MethodError`. With the annotation on the array only,
# FFTW's methods are strictly more specific and win, as AbstractFFTs intends.
AbstractFFTs.plan_rfft(x::AbstractArray{T,N}, region; kwargs...) where {T<:Real,N} =
    _plan_rfft(x, _region(region); kwargs...)

AbstractFFTs.plan_brfft(x::AbstractArray{T,N}, len::Integer, region; kwargs...) where {T,N} =
    _plan_brfft(x, Int(len), _region(region); kwargs...)

_region(r::RegionTypes) = r
_region(r) = collect(Int, r)

function _plan_rfft(
    x::AbstractArray{T,N},
    region::RegionTypes;
    BLUESTEIN_CUTOFF=DEFAULT_BLUESTEIN_CUTOFF, num_threads::Int=Threads.nthreads(), _kwargs...
) where {T<:Real,N}
    M = length(region)
    if M == 1
        R1 = Int(region[1])
        n = size(x, R1)
        # For even length problems, we solve the real problem with
        # two n/2 complex FFTs followed by a butterfly. For odd size
        # problems, we just solve the problem as a single complex
        nn = iseven(n) ? n >> 1 : n
        g = CallGraph{Complex{T}}(nn, BLUESTEIN_CUTOFF)
        return FFTAPlan_re{Complex{T},1}((g,), R1, FFT_FORWARD, n; num_threads)
    else
        # real transform along the first region dimension (half length when
        # even), complex transforms along the others
        RM = _sort(region)
        n = size(x, RM[1])
        nn = iseven(n) ? n >> 1 : n
        cg = (CallGraph{Complex{T}}(nn, BLUESTEIN_CUTOFF),
              ntuple(i -> CallGraph{Complex{T}}(size(x, RM[i + 1]), BLUESTEIN_CUTOFF), Val(M - 1))...)
        return FFTAPlan_re{Complex{T},M}(cg, RM, FFT_FORWARD, n; num_threads)
    end
end

function _plan_brfft(
    x::AbstractArray{T,N},
    len::Int,
    region::RegionTypes;
    BLUESTEIN_CUTOFF=DEFAULT_BLUESTEIN_CUTOFF, num_threads::Int=Threads.nthreads(), _kwargs...
) where {T,N}
    M = length(region)
    if M == 1
        # For even length problems, we solve the real problem with
        # two n/2 complex FFTs followed by a butterfly. For odd size
        # problems, we just solve the problem as a single complex
        R1 = Int(region[1])
        nn = iseven(len) ? len >> 1 : len
        g = CallGraph{T}(nn, BLUESTEIN_CUTOFF)
        return FFTAPlan_re{T,1}((g,), R1, FFT_BACKWARD, len; num_threads)
    else
        RM = _sort(region)
        nn = iseven(len) ? len >> 1 : len
        cg = (CallGraph{T}(nn, BLUESTEIN_CUTOFF),
              ntuple(i -> CallGraph{T}(size(x, RM[i + 1]), BLUESTEIN_CUTOFF), Val(M - 1))...)
        return FFTAPlan_re{T,M}(cg, RM, FFT_BACKWARD, len; num_threads)
    end
end


# Multiplication
## mul!
### Complex
#### 1D plan 1D array
function LinearAlgebra.mul!(y::AbstractVector{U}, p::FFTAPlan_cx{T,1}, x::AbstractVector{T}) where {T,U}
    Base.require_one_based_indexing(x, y)
    if axes(x) != axes(y)
        throw(DimensionMismatch("input array has axes $(axes(x)), but output array has axes $(axes(y))"))
    end
    if size(p) != size(x)
        throw(DimensionMismatch("plan has axes $(size(p)), but input array has axes $(size(x))"))
    end
    fft_kernel!(y, x, 1, 1, p.dir, p.callgraph[1][1].type, p.callgraph[1], 1)
    return y
end

#### 1D plan ND array
function LinearAlgebra.mul!(y::AbstractArray{U,N}, p::FFTAPlan_cx{T,1}, x::AbstractArray{T,N}) where {T,U,N}
    Base.require_one_based_indexing(x, y)

    ax_x, ax_y = axes(x), axes(y)
    if ax_x != ax_y
        throw(DimensionMismatch("input array has axes $ax_x, but output array has axes $ax_y"))
    end

    R1 = only(p.region)
    plen, xlen = size(p, 1), size(x, R1)
    if plen != xlen
        throw(DimensionMismatch("plan has size $plen, but input array has size $xlen along region $R1"))
    end

    if @generated
        quote
            Base.Cartesian.@nif $N d -> (d == R1) dim -> (_mul_loop!(y, x, p, Val(dim)))
        end
    else
        _mul_loop!(y, x, p, Val(R1))
    end
    return y
end

"""
$(TYPEDSIGNATURES)
Call `f(worker, Ipre, Ipost)` for every pencil `A[Ipre, :, Ipost]` along
dimension `dim` of `A`. When the plan has several workers and
the transform is large enough (`THREAD_THRESHOLD` elements), the pencils are
split into one contiguous chunk per worker and the chunks run as parallel
tasks, each using its own worker (so its own call-graph workspace and
buffers). Results do not depend on the number of threads.
"""
function _foreach_pencil(f::F, A::AbstractArray{<:Any,N}, ::Val{dim}, workers::Vector{W}) where {F,N,dim,W}
    Rpre  = CartesianIndices(ntuple(Base.Fix1(size, A), Val(dim - 1)))
    Rpost = CartesianIndices(ntuple(i -> size(A, dim + i), Val(N - dim)))
    npre = length(Rpre)
    total = npre * length(Rpost)
    nt = min(length(workers), total ÷ 2)
    if nt > 1 && total * size(A, dim) >= THREAD_THRESHOLD
        Base.@sync for c in 1:nt
            lo = (c - 1) * total ÷ nt + 1
            hi = c * total ÷ nt
            # (the task's worker gets its own name: a variable shared with the
            # serial branch below would be captured by every task)
            wc = workers[c]
            Threads.@spawn for k in lo:hi
                ipost, ipre = divrem(k - 1, npre)
                f(wc, Rpre[ipre + 1], Rpost[ipost + 1])
            end
        end
    else
        w1 = workers[1]
        for Ipost in Rpost, Ipre in Rpre
            f(w1, Ipre, Ipost)
        end
    end
    return nothing
end

function _mul_loop!(
    y::AbstractArray{U,N},
    x::AbstractArray{T,N},
    p::FFTAPlan_cx{T,1},
    ::Val{R}
) where {T,U,N,R}
    dir = p.dir
    _foreach_pencil(x, Val(R), p.workers) do w, Ipre, Ipost
        cg = w.callgraph[1]
        @views fft_kernel!(y[Ipre,:,Ipost], x[Ipre,:,Ipost], 1, 1, dir, cg[1].type, cg, 1)
    end
end

#### ND plan ND array
function LinearAlgebra.mul!(
    out::AbstractArray{U,N},
    p::FFTAPlan_cx{T,N},
    X::AbstractArray{T,N}
) where {T,U,N}
    Base.require_one_based_indexing(out, X)
    if size(out) != size(X)
        throw(DimensionMismatch("input array has axes $(axes(X)), but output array has axes $(axes(out))"))
    elseif size(p) != size(X)
        throw(DimensionMismatch("plan has size $(size(p)), but input array has size $(size(X))"))
    elseif !region_isvalid(p.region, N)
        throw(DimensionMismatch("Plan region is outside array dimensions."))
    end

    dir = p.dir
    workers = p.workers

    copyto!(out, X) # operate in-place on output array

    if @generated
        quote
            Base.Cartesian.@nexprs $N dim -> fft_along_dim!(out, workers, dim, dir, Val(dim))
        end
    else
        for dim in 1:N
            fft_along_dim!(out, workers, dim, dir, Val(dim))
        end
    end

    return out
end

#### MD plan ND array (M<N)
function LinearAlgebra.mul!(
    out::AbstractArray{U,N},
    p::FFTAPlan_cx{T,M},
    X::AbstractArray{T,N}
) where {T,U,N,M}
    Base.require_one_based_indexing(out, X)
    if size(out) != size(X)
        throw(DimensionMismatch("input array has axes $(axes(X)), but output array has axes $(axes(out))"))
    elseif !region_isvalid(p.region, M, N)
        throw(DimensionMismatch("Region is invalid."))
    elseif M > N || first(p.region) < 1 || last(p.region) > N
        throw(DimensionMismatch("Plan region is outside array dimensions."))
    end

    copyto!(out, X) # operate in-place on output array

    _execute_mdfft!(out, p.workers, p.dir, p.region)

    return out
end

@noinline function _execute_mdfft!(
    out::AbstractArray{U,N},
    workers::Vector{<:Worker{T,M}},
    dir::Direction,
    @nospecialize(region::RegionTypes),
) where {T,U,N,M}
    if @generated
        quote
            k = 1
            # region is assumed to be pre-sorted during planning
            Base.Cartesian.@nexprs $N dim -> begin
                if region[k] == dim
                    fft_along_dim!(out, workers, k, dir, Val(dim))
                    k = min(k + 1, M)
                end
            end
            return nothing
        end
    else
        for k in 1:M
            fft_along_dim!(out, workers, k, dir, Val(Int(region[k])))
        end
    end
end

"""
$(TYPEDSIGNATURES)
Transform every pencil of `A` along dimension `dim` in place, using the `k`-th
call graph of the workers. The kernels read the (strided) pencil directly and
write the worker's output buffer, which is then copied back.
"""
function fft_along_dim!(
    A::AbstractArray{U,N},
    workers::Vector{<:Worker{T,M}},
    k::Int, d::Direction,
    ::Val{dim}
) where {T <: Complex{<:AbstractFloat}, U, N, M, dim}
    n = size(A, dim)
    _foreach_pencil(A, Val(dim), workers) do w, Ipre, Ipost
        cg = w.callgraph[k]
        obuf = w.obuf
        pencil = @view A[Ipre, :, Ipost]
        fft_kernel!(obuf, pencil, 1, 1, d, cg[1].type, cg, 1)
        @inbounds for j in 1:n
            pencil[j] = obuf[j]
        end
    end
end

region_isvalid(r::Int, N::Int, _::Int=0) = r == N == 1
region_isvalid(r::AbstractVector{Int}, N::Int) = r == 1:N
region_isvalid(r::AbstractRange{Int}, M::Int, _::Int) = issorted(r) && length(r) == M
function region_isvalid(r::NTuple{M,Int}, N::Int) where M
    isvalid = M == N
    for i in 1:M
        isvalid &= (r[i] == i)
    end
    isvalid
end
function region_isvalid(r::Union{AbstractVector{Int},NTuple{<:Any,Int}}, M::Int, _::Int)
    isvalid = length(r) == M
    maybe_p = Iterators.peel(r)
    isnothing(maybe_p) && return isvalid
    p, rest = maybe_p
    for n in rest
        isvalid = isvalid && (p < n)
        p = n
    end
    isvalid
end

## *
### Complex
function Base.:*(p::FFTAPlan_cx{T,1}, x::AbstractVector{T}) where {T<:Complex}
    y = similar(x)
    LinearAlgebra.mul!(y, p, x)
    y
end

function Base.:*(p::FFTAPlan_cx{T,N1}, x::AbstractArray{T,N2}) where {T<:Complex,N1,N2}
    y = similar(x)
    LinearAlgebra.mul!(y, p, x)
    y
end

### Real
# Real transforms are computed from complex transforms of half (even `n`) or
# full (odd `n`) length; see `_rfft_pencil!`/`_brfft_pencil!`. All entry points
# funnel into `mul!`, which is allocation-free for 1D plans (the scratch space
# lives in the plan) and applies the same 1D kernel to every pencil along the
# transform dimension of an N-d array.

function _check_re_dims(y, p::FFTAPlan_re, x, d1::Int, fwd::Bool)
    # `x` is the input, `y` the output; the real-length dimension is `d1`.
    rlen, clen = p.flen, p.flen ÷ 2 + 1
    xr, yr = fwd ? (rlen, clen) : (clen, rlen)
    if ndims(x) != ndims(y)
        throw(DimensionMismatch("input has $(ndims(x)) dimensions, output has $(ndims(y))"))
    elseif size(x, d1) != xr
        throw(DimensionMismatch("real 1D plan has size $rlen. Dimension of input array along region $d1 should have size $xr, but has size $(size(x, d1))"))
    elseif size(y, d1) != yr
        throw(DimensionMismatch("output array should have size $yr along region $d1, but has size $(size(y, d1))"))
    end
    for i in 1:ndims(x)
        if i != d1 && size(x, i) != size(y, i)
            throw(DimensionMismatch("input array has size $(size(x)), but output array has size $(size(y))"))
        end
    end
    return nothing
end

# Forward real-to-complex transform of one pencil: `x` real of length `n`,
# `y` complex of length `n ÷ 2 + 1`.
function _rfft_pencil!(y::AbstractVector{T}, x::AbstractVector{<:Real}, w::Worker{T}, n::Int) where {T<:Complex}
    R = real(T)
    cg = w.callgraph[1]
    buf = w.buf
    if iseven(n)
        # Solve the rfft problem by splitting the input into even and odd parts
        # and solving them simultaneously as a single (complex) fft of half
        # the size, see equations (6)-(8) of Sorensen, H. V., D. Jones, Michael
        # Heideman, and C. Burrus. "Real-valued fast Fourier transform
        # algorithms." IEEE Transactions on acoustics, speech, and signal
        # processing 35, no. 6 (2003): 849-863.
        m = n >> 1
        @inbounds for j in 1:m
            buf[j] = T(x[2j - 1], x[2j])
        end
        fft_kernel!(view(y, 1:m), buf, 1, 1, FFT_FORWARD, cg[1].type, cg, 1)

        # Construct the result by first constructing the elements of the
        # real and imaginary part, followed by the usual radix-2 assembly,
        # see eq (9). The twiddle is for `n`, not `m`, so it is recomputed.
        z1 = singleton_params(-one(R) / n)
        wj = cispi(-R(2) / n)
        @inbounds begin
            y1 = y[1]
            y[1]     = real(y1) + imag(y1)
            y[m + 1] = real(y1) - imag(y1)
            for j in 2:((m >> 1) + 1)
                yj  = y[j]
                ymj = y[m - j + 2]
                XX = R(0.5) * ( yj + conj(ymj))
                XY = R(0.5) * (-yj + conj(ymj)) * im
                y[j]         =      XX + wj * XY
                y[m - j + 2] = conj(XX - wj * XY)
                wj = singleton_step(wj, z1)
            end
        end
    else
        # Odd length: run the full transform on the real input (the kernels
        # accept real input; the DFT leaf exploits its symmetry) and keep the
        # first half.
        fft_kernel!(buf, x, 1, 1, FFT_FORWARD, cg[1].type, cg, 1)
        @inbounds for j in 1:(n ÷ 2 + 1)
            y[j] = buf[j]
        end
    end
    return y
end

# Backward complex-to-real transform of one pencil: `y` complex of length
# `n ÷ 2 + 1`, `x` real of length `n`.
function _brfft_pencil!(x::AbstractVector{<:Real}, y::AbstractVector{T}, w::Worker{T}, n::Int) where {T<:Complex}
    R = real(T)
    cg = w.callgraph[1]
    buf = w.buf
    if iseven(n)
        # Inverse of the even-length trick in `_rfft_pencil!`.
        m = n >> 1
        tmp = view(buf, 1:m)
        out = view(buf, m + 1:2m)
        z1 = singleton_params(one(R) / n)
        wj = cispi(R(2) / n)
        @inbounds begin
            tmp[1] = T(real(y[1]) + real(y[m + 1]), real(y[1]) - real(y[m + 1]))
            for j in 2:((m >> 1) + 1)
                XX =       y[j] + conj(y[m - j + 2])
                XY = wj * (y[j] - conj(y[m - j + 2]))
                tmp[j]         =      XX + im * XY
                tmp[m - j + 2] = conj(XX - im * XY)
                wj = singleton_step(wj, z1)
            end
        end
        fft_kernel!(out, tmp, 1, 1, FFT_BACKWARD, cg[1].type, cg, 1)
        @inbounds for j in 1:m
            x[2j - 1] = real(out[j])
            x[2j]     = imag(out[j])
        end
    else
        # Odd length: rebuild the conjugate-symmetric spectrum and transform.
        h = n ÷ 2 + 1
        tmp = view(buf, 1:n)
        out = view(buf, n + 1:2n)
        @inbounds for j in 1:h
            tmp[j] = y[j]
        end
        @inbounds for j in h + 1:n
            tmp[j] = conj(y[n - j + 2])
        end
        fft_kernel!(out, tmp, 1, 1, FFT_BACKWARD, cg[1].type, cg, 1)
        @inbounds for j in 1:n
            x[j] = real(out[j])
        end
    end
    return x
end

# Unit-stride pencils are handed to the kernels directly; strided ones are
# copied to the worker's contiguous buffers first. The copy costs one pass over
# the pencil but keeps the kernel's scattered reads within cache lines — for a
# stride of a cache line or more the kernel is otherwise 1.2–1.3× slower. The
# test is sufficient rather than exact (a 1×N array's dim-2 pencils are
# contiguous but still take the copy): a wasted copy, never a wrong answer.
_unit_stride(v::StridedArray) = stride(v, 1) == 1
_unit_stride(v) = false

# Apply `kernel!(y_pencil, x_pencil, worker, flen)` along dimension `R` of `x` and `y`.
function _re_pencil_loop!(kernel!::F, y::AbstractArray{<:Any,N}, x::AbstractArray{<:Any,N}, p::FFTAPlan_re, ::Val{R}) where {F,N,R}
    n = p.flen
    if R == 1 && _unit_stride(x) && _unit_stride(y)
        _foreach_pencil(x, Val(R), p.workers) do w, Ipre, Ipost
            @views kernel!(y[Ipre, :, Ipost], x[Ipre, :, Ipost], w, n)
        end
    else
        _foreach_pencil(x, Val(R), p.workers) do w, Ipre, Ipost
            xin, yout = _pencil_buffers(kernel!, w)   # per worker: no sharing between tasks
            copyto!(xin, @view x[Ipre, :, Ipost])
            kernel!(yout, xin, w, n)
            copyto!(@view(y[Ipre, :, Ipost]), yout)
        end
    end
    return y
end
# forward: real in, complex out; backward: complex in, real out
_pencil_buffers(::typeof(_rfft_pencil!), w::Worker) = (w.rbuf, w.cbuf)
_pencil_buffers(::typeof(_brfft_pencil!), w::Worker) = (w.cbuf, w.rbuf)
# forward: real in, complex out; backward: complex in, real out
_pencil_buffers(::typeof(_rfft_pencil!), p::FFTAPlan_re) = (p.rbuf, p.cbuf)
_pencil_buffers(::typeof(_brfft_pencil!), p::FFTAPlan_re) = (p.cbuf, p.rbuf)

# Dispatch on the transform dimension so that the loop above is type-stable.
function _re_along_dim!(kernel!::F, y::AbstractArray{<:Any,N}, x::AbstractArray{<:Any,N}, p::FFTAPlan_re, d::Int) where {F,N}
    if @generated
        quote
            Base.Cartesian.@nif $N dim -> (d == dim) dim -> (_re_pencil_loop!(kernel!, y, x, p, Val(dim)))
        end
    else
        _re_pencil_loop!(kernel!, y, x, p, Val(d))
    end
end

function _cx_along_dim!(A::AbstractArray{<:Any,N}, workers::Vector{<:Worker{T,M}}, k::Int, dir::Direction, d::Int) where {N,T,M}
    if @generated
        quote
            Base.Cartesian.@nif $N dim -> (d == dim) dim -> (fft_along_dim!(A, workers, k, dir, Val(dim)))
        end
    else
        fft_along_dim!(A, workers, k, dir, Val(d))
    end
end

## mul!
#### 1D plan
##### Forward
function LinearAlgebra.mul!(y::AbstractArray{T,N}, p::FFTAPlan_re{T,1}, x::AbstractArray{<:Real,N}) where {T<:Complex,N}
    if p.dir !== FFT_FORWARD
        throw(ArgumentError("only FFT_FORWARD supported for real $(N == 1 ? "vectors" : "arrays")"))
    end
    Base.require_one_based_indexing(x, y)
    d1 = only(p.region)
    _check_re_dims(y, p, x, d1, true)
    if N == 1
        _rfft_pencil!(y, x, p.workers[1], p.flen)
    else
        _re_along_dim!(_rfft_pencil!, y, x, p, d1)
    end
    return y
end

##### Backward
function LinearAlgebra.mul!(y::AbstractArray{<:Real,N}, p::FFTAPlan_re{T,1}, x::AbstractArray{T,N}) where {T<:Complex,N}
    if p.dir !== FFT_BACKWARD
        throw(ArgumentError("only FFT_BACKWARD supported for complex $(N == 1 ? "vectors" : "arrays")"))
    end
    Base.require_one_based_indexing(x, y)
    d1 = only(p.region)
    _check_re_dims(y, p, x, d1, false)
    if N == 1
        _brfft_pencil!(y, x, p.workers[1], p.flen)
    else
        _re_along_dim!(_brfft_pencil!, y, x, p, d1)
    end
    return y
end

#### M-d plan (M >= 2)
# The real transform is taken along the first region dimension, then complex
# transforms along the others (forward), or the reverse (backward).
function _check_re_other_dims(p::FFTAPlan_re{T,M}, x, fwd::Bool) where {T,M}
    for k in 2:M
        d = p.region[k]
        if size(x, d) != size(p, k)
            expected = ntuple(i -> i == 1 && !fwd ? size(p, 1) ÷ 2 + 1 : size(p, i), Val(M))
            actual = ntuple(i -> size(x, p.region[i]), Val(M))
            throw(DimensionMismatch("real $(M)D plan has size $(size(p)). Transform dimensions of input array are $actual but should be $expected"))
        end
    end
end

##### Forward
function LinearAlgebra.mul!(y::AbstractArray{T,N}, p::FFTAPlan_re{T,M}, x::AbstractArray{<:Real,N}) where {T<:Complex,N,M}
    if p.dir !== FFT_FORWARD
        throw(ArgumentError("only FFT_FORWARD supported for real arrays"))
    end
    Base.require_one_based_indexing(x, y)
    d1 = p.region[1]
    _check_re_dims(y, p, x, d1, true)
    _check_re_other_dims(p, x, true)
    _re_along_dim!(_rfft_pencil!, y, x, p, d1)
    for k in 2:M
        _cx_along_dim!(y, p.workers, k, FFT_FORWARD, Int(p.region[k]))
    end
    return y
end

##### Backward
function LinearAlgebra.mul!(y::AbstractArray{<:Real,N}, p::FFTAPlan_re{T,M}, x::AbstractArray{T,N}) where {T<:Complex,N,M}
    if p.dir !== FFT_BACKWARD
        throw(ArgumentError("only FFT_BACKWARD supported for complex arrays"))
    end
    Base.require_one_based_indexing(x, y)
    d1 = p.region[1]
    _check_re_dims(y, p, x, d1, false)
    _check_re_other_dims(p, x, false)
    tmp = copy(x)   # the complex passes must not modify the input
    for k in M:-1:2
        _cx_along_dim!(tmp, p.workers, k, FFT_BACKWARD, Int(p.region[k]))
    end
    _re_along_dim!(_brfft_pencil!, y, tmp, p, d1)
    return y
end

## *
function Base.:*(p::FFTAPlan_re{T,M}, x::AbstractArray{<:Real,N}) where {T<:Complex,M,N}
    if p.dir !== FFT_FORWARD
        throw(ArgumentError("only FFT_FORWARD supported for real $(N == 1 ? "vectors" : "arrays")"))
    end
    d1 = first(p.region)
    y = similar(x, T, ntuple(i -> i == d1 ? p.flen ÷ 2 + 1 : size(x, i), Val(N)))
    return LinearAlgebra.mul!(y, p, x)
end

function Base.:*(p::FFTAPlan_re{T,M}, x::AbstractArray{T,N}) where {T<:Complex,M,N}
    if p.dir !== FFT_BACKWARD
        throw(ArgumentError("only FFT_BACKWARD supported for complex $(N == 1 ? "vectors" : "arrays")"))
    end
    d1 = first(p.region)
    y = similar(x, real(T), ntuple(i -> i == d1 ? p.flen : size(x, i), Val(N)))
    return LinearAlgebra.mul!(y, p, x)
end


# Inverse plans
# `AbstractFFTs.inv(p)` calls `plan_inv(p)` once and caches it in `p.pinv`;
# `p \\ x` and `ldiv!(y, p, x)` go through `inv(p)`.

function AbstractFFTs.plan_inv(p::FFTAPlan_cx{T,N,R}) where {T,N,R}
    dir = p.dir === FFT_FORWARD ? FFT_BACKWARD : FFT_FORWARD
    cutoff = p.callgraph[1].BLUESTEIN_CUTOFF
    cg = ntuple(i -> CallGraph{T}(size(p, i), cutoff), Val(N))
    q = FFTAPlan_cx{T,N,R}(cg, p.region, dir, _workers(cg, Val(N), _num_threads(p), 0))
    AbstractFFTs.ScaledPlan(q, _normalization(p))
end

function AbstractFFTs.plan_inv(p::FFTAPlan_re{T,N,R,S}) where {T,N,R,S}
    dir = p.dir === FFT_FORWARD ? FFT_BACKWARD : FFT_FORWARD
    cutoff = p.callgraph[1].BLUESTEIN_CUTOFF
    n = p.flen
    nn = iseven(n) ? n >> 1 : n
    cg = (CallGraph{T}(nn, cutoff), ntuple(i -> CallGraph{T}(size(p, i + 1), cutoff), Val(N - 1))...)
    q = FFTAPlan_re{T,N,R}(cg, p.region, dir, n, _workers(cg, Val(N), _num_threads(p), _re_buflen(n, dir), n))
    AbstractFFTs.ScaledPlan(q, _normalization(p))
end

# 1 / (product of the transform lengths); `size(p)` lists exactly those.
_normalization(p::FFTAPlan{T,N}) where {T,N} =
    AbstractFFTs.normalization(real(T), size(p), ntuple(identity, Val(N)))

# In-place plans
# (`AbstractFFTs.fft!`/`bfft!`/`ifft!` build these; the internal kernel is
# `fft_kernel!` so that it does not shadow `AbstractFFTs.fft!`.)
# FFTA's kernels are out of place, so an in-place plan wraps an ordinary plan
# and, when the input and output alias, transforms a copy of the input held
# in the plan's buffer.

mutable struct FFTAPlan_inplace{T,N,M,P<:FFTAPlan_cx{T,M}} <: FFTAPlan{T,M}
    const p::P
    const buf::Array{T,N}   # copy of the input when it aliases the output
    pinv::AbstractFFTs.ScaledPlan
    FFTAPlan_inplace(p::P, buf::Array{T,N}) where {T,N,M,P<:FFTAPlan_cx{T,M}} = new{T,N,M,P}(p, buf)
end
Base.size(p::FFTAPlan_inplace) = size(p.p)
Base.size(p::FFTAPlan_inplace, i::Int) = size(p.p, i)

AbstractFFTs.plan_fft!(x::AbstractArray{T,N}, region; kwargs...) where {T<:Complex,N} =
    FFTAPlan_inplace(_plan_fft(x, region, FFT_FORWARD; kwargs...), Array{T,N}(undef, size(x)))
AbstractFFTs.plan_bfft!(x::AbstractArray{T,N}, region; kwargs...) where {T<:Complex,N} =
    FFTAPlan_inplace(_plan_fft(x, region, FFT_BACKWARD; kwargs...), Array{T,N}(undef, size(x)))

function LinearAlgebra.mul!(y::AbstractArray{T,N}, ip::FFTAPlan_inplace{T,N}, x::AbstractArray{T,N}) where {T,N}
    if y === x
        if size(ip.buf) != size(x)
            throw(DimensionMismatch("in-place plan was created for size $(size(ip.buf)), input has size $(size(x))"))
        end
        copyto!(ip.buf, x)
        LinearAlgebra.mul!(y, ip.p, ip.buf)
    else
        LinearAlgebra.mul!(y, ip.p, x)
    end
    return y
end
Base.:*(ip::FFTAPlan_inplace{T}, x::AbstractArray{T}) where {T} = LinearAlgebra.mul!(x, ip, x)

function AbstractFFTs.plan_inv(ip::FFTAPlan_inplace)
    s = AbstractFFTs.plan_inv(ip.p)
    AbstractFFTs.ScaledPlan(FFTAPlan_inplace(s.p, ip.buf), s.scale)
end
