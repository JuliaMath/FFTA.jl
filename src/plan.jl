# Plans

abstract type FFTAPlan{T,N} <: AbstractFFTs.Plan{T} end

const RegionTypes{N} = Union{Int,AbstractVector{Int},NTuple{N,Int}}

# The plan types are mutable so that `AbstractFFTs.inv` can cache the inverse
# plan in the initially undefined `pinv` field (see `plan_inv`), as FFTW.jl does.

mutable struct FFTAPlan_cx{T,N,R<:RegionTypes{N}} <: FFTAPlan{T,N}
    const callgraph::NTuple{N,CallGraph{T}}
    const region::R
    const dir::Direction
    pinv::AbstractFFTs.ScaledPlan
    FFTAPlan_cx{T,N,R}(cg::NTuple{N,CallGraph{T}}, r::R, dir::Direction) where {T,N,R<:RegionTypes{N}} =
        new{T,N,R}(cg, r, dir)
end
function FFTAPlan_cx{T,N}(
    cg::NTuple{N,CallGraph{T}}, r::R,
    dir::Direction
) where {T,N,R<:RegionTypes{N}}
    FFTAPlan_cx{T,N,R}(cg, r, dir)
end

mutable struct FFTAPlan_re{T,N,R<:RegionTypes{N},S<:Real} <: FFTAPlan{T,N}
    const callgraph::NTuple{N,CallGraph{T}}
    const region::R
    const dir::Direction
    const flen::Int
    const buf::Vector{T}    # scratch for the real<->complex packing, see `_re_buflen`
    const rbuf::Vector{S}   # contiguous copy of a strided real pencil, see `_re_pencil_loop!`
    const cbuf::Vector{T}   # contiguous copy of a strided complex pencil
    pinv::AbstractFFTs.ScaledPlan
    FFTAPlan_re{T,N,R,S}(cg::NTuple{N,CallGraph{T}}, r::R, dir::Direction, flen::Int, buf::Vector{T}, rbuf::Vector{S}, cbuf::Vector{T}) where {T,N,R<:RegionTypes{N},S<:Real} =
        new{T,N,R,S}(cg, r, dir, flen, buf, rbuf, cbuf)
end
function FFTAPlan_re{T,N}(
    cg::NTuple{N,CallGraph{T}}, r::R,
    dir::Direction, flen::Int
) where {T,N,R<:RegionTypes{N}}
    S = real(T)
    buf = Vector{T}(undef, _re_buflen(flen, dir))
    FFTAPlan_re{T,N,R,S}(cg, r, dir, flen, buf, Vector{S}(undef, flen), Vector{T}(undef, flen ÷ 2 + 1))
end

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
    BLUESTEIN_CUTOFF=DEFAULT_BLUESTEIN_CUTOFF, _kwargs...
) where {T<:Complex,N}
    M = length(region)
    if M == 1
        R1 = Int(region[1])
        g = CallGraph{T}(size(x, R1), BLUESTEIN_CUTOFF)
        return FFTAPlan_cx{T,1}((g,), R1, dir)
    elseif M == 2
        R2 = _sort(region)
        g1 = CallGraph{T}(size(x, R2[1]), BLUESTEIN_CUTOFF)
        g2 = CallGraph{T}(size(x, R2[2]), BLUESTEIN_CUTOFF)
        return FFTAPlan_cx{T,2}((g1, g2), R2, dir)
    else
        RM = _sort(region)
        return FFTAPlan_cx{T,M}(
            ntuple(i -> CallGraph{T}(size(x, RM[i]), BLUESTEIN_CUTOFF), Val(M)),
            RM, dir
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
    BLUESTEIN_CUTOFF=DEFAULT_BLUESTEIN_CUTOFF, _kwargs...
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
        return FFTAPlan_re{Complex{T},1}((g,), R1, FFT_FORWARD, n)
    elseif M == 2
        R2 = _sort(region)
        n = size(x, R2[1])
        nn = iseven(n) ? n >> 1 : n
        g1 = CallGraph{Complex{T}}(nn, BLUESTEIN_CUTOFF)
        g2 = CallGraph{Complex{T}}(size(x, R2[2]), BLUESTEIN_CUTOFF)
        return FFTAPlan_re{Complex{T},2}((g1, g2), R2, FFT_FORWARD, n)
    else
        throw(ArgumentError("only supports 1D and 2D FFTs"))
    end
end

function _plan_brfft(
    x::AbstractArray{T,N},
    len::Int,
    region::RegionTypes;
    BLUESTEIN_CUTOFF=DEFAULT_BLUESTEIN_CUTOFF, _kwargs...
) where {T,N}
    M = length(region)
    if M == 1
        # For even length problems, we solve the real problem with
        # two n/2 complex FFTs followed by a butterfly. For odd size
        # problems, we just solve the problem as a single complex
        R1 = Int(region[1])
        nn = iseven(len) ? len >> 1 : len
        g = CallGraph{T}(nn, BLUESTEIN_CUTOFF)
        return FFTAPlan_re{T,1}((g,), R1, FFT_BACKWARD, len)
    elseif M == 2
        R2 = _sort(region)
        nn = iseven(len) ? len >> 1 : len
        g1 = CallGraph{T}(nn, BLUESTEIN_CUTOFF)
        g2 = CallGraph{T}(size(x, R2[2]), BLUESTEIN_CUTOFF)
        return FFTAPlan_re{T,2}((g1, g2), R2, FFT_BACKWARD, len)
    else
        throw(ArgumentError("only supports 1D and 2D FFTs"))
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

function _mul_loop!(
    y::AbstractArray{U,N},
    x::AbstractArray{T,N},
    p::FFTAPlan_cx{T,1},
    ::Val{R}
) where {T,U,N,R}
    Rpre  = CartesianIndices(ntuple(Base.Fix1(size, x),  Val(R - 1)))
    Rpost = CartesianIndices(ntuple(i -> size(x, R + i), Val(N - R)))
    cg = p.callgraph[1]
    t = cg[1].type
    for Ipost in Rpost, Ipre in Rpre
        @views fft_kernel!(y[Ipre,:,Ipost], x[Ipre,:,Ipost], 1, 1, p.dir, t, cg, 1)
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

    sz = size(X)
    max_sz = maximum(sz)
    obuf = Vector{T}(undef, max_sz)
    ibuf = Vector{T}(undef, max_sz)
    sizehint!(obuf, max_sz) # not guaranteed but hopefully prevents allocations
    sizehint!(ibuf, max_sz)
    dir = p.dir

    copyto!(out, X) # operate in-place on output array

    if @generated
        quote
            Base.Cartesian.@nexprs $N dim -> begin
                n = size(out, dim)
                resize!(obuf, n)
                resize!(ibuf, n)
                cg = p.callgraph[dim]

                fft_along_dim!(out, ibuf, obuf, cg, dir, Val(dim))
            end
        end
    else
        for dim in 1:N
            n = size(out, dim)
            resize!(obuf, n)
            resize!(ibuf, n)
            cg = p.callgraph[dim]

            fft_along_dim!(out, ibuf, obuf, cg, dir, Val(dim))
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

    max_sz = maximum(Base.Fix1(size, out), p.region)
    obuf = Vector{T}(undef, max_sz)
    ibuf = Vector{T}(undef, max_sz)
    sizehint!(obuf, max_sz) # not guaranteed but hopefully prevents allocations
    sizehint!(ibuf, max_sz)

    copyto!(out, X) # operate in-place on output array

    _execute_mdfft!(out, ibuf, obuf, p.dir, p.region, p.callgraph)

    return out
end

@noinline function _execute_mdfft!(
    out::AbstractArray{U,N},
    ibuf::Vector{T}, obuf::Vector{T},
    dir::Direction,
    @nospecialize(region::RegionTypes),
    @nospecialize(callgraphs::NTuple)
) where {T,U,N}

    M = length(region)
    if @generated
        quote
            k = 1
            # region is assumed to be pre-sorted during planning
            Base.Cartesian.@nexprs $N dim -> begin
                if region[k] == dim
                    n = size(out, dim)
                    resize!(obuf, n)
                    resize!(ibuf, n)
                    cg = callgraphs[k]

                    fft_along_dim!(out, ibuf, obuf, cg, dir, Val(dim))

                    k = min(k + 1, M)
                end
            end
            return nothing
        end
    else
        for dim in 1:M
            pdim = region[dim]
            n = size(out, pdim)
            resize!(obuf, n)
            resize!(ibuf, n)
            cg = callgraphs[dim]

            fft_along_dim!(out, ibuf, obuf, cg, dir, Val(pdim))
        end
    end
end

function fft_along_dim!(
    A::AbstractArray{U,N},
    ibuf::Vector{T}, obuf::Vector{T},
    cg::CallGraph{T}, d::Direction,
    ::Val{dim}
) where {T <: Complex{<:AbstractFloat}, U, N, dim}

    Rpre  = CartesianIndices(ntuple(Base.Fix1(size, A),    Val(dim - 1)))
    Rpost = CartesianIndices(ntuple(i -> size(A, dim + i), Val(N - dim)))
    t = cg[1].type
    cols = eachindex(axes(A, dim), ibuf, obuf)

    for Ipost in Rpost, Ipre in Rpre
        for j in cols
            ibuf[j] = A[Ipre, j, Ipost]
        end
        fft_kernel!(obuf, ibuf, 1, 1, d, t, cg, 1)
        for j in cols
            A[Ipre, j, Ipost] = obuf[j]
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
function _rfft_pencil!(y::AbstractVector{T}, x::AbstractVector{<:Real}, p::FFTAPlan_re{T}) where {T<:Complex}
    n = p.flen
    R = real(T)
    cg = p.callgraph[1]
    buf = p.buf
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
function _brfft_pencil!(x::AbstractVector{<:Real}, y::AbstractVector{T}, p::FFTAPlan_re{T}) where {T<:Complex}
    n = p.flen
    R = real(T)
    cg = p.callgraph[1]
    buf = p.buf
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
# copied to the plan's contiguous buffers first. The copy costs one pass over
# the pencil but keeps the kernel's scattered reads within cache lines — for a
# stride of a cache line or more the kernel is otherwise 1.2–1.3× slower.
_unit_stride(v::StridedArray) = stride(v, 1) == 1
_unit_stride(v) = false

# Apply `kernel!(y_pencil, x_pencil, p)` along dimension `R` of `x` and `y`.
function _re_pencil_loop!(kernel!::F, y::AbstractArray{<:Any,N}, x::AbstractArray{<:Any,N}, p::FFTAPlan_re, ::Val{R}) where {F,N,R}
    Rpre  = CartesianIndices(ntuple(Base.Fix1(size, x),  Val(R - 1)))
    Rpost = CartesianIndices(ntuple(i -> size(x, R + i), Val(N - R)))
    contiguous = R == 1 && _unit_stride(x) && _unit_stride(y)
    if contiguous
        for Ipost in Rpost, Ipre in Rpre
            @views kernel!(y[Ipre, :, Ipost], x[Ipre, :, Ipost], p)
        end
    else
        xin, yout = _pencil_buffers(kernel!, p)
        for Ipost in Rpost, Ipre in Rpre
            xv = @view x[Ipre, :, Ipost]
            yv = @view y[Ipre, :, Ipost]
            copyto!(xin, xv)
            kernel!(yout, xin, p)
            copyto!(yv, yout)
        end
    end
    return y
end
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

function _cx_along_dim!(A::AbstractArray{<:Any,N}, cg::CallGraph{T}, dir::Direction, d::Int) where {N,T}
    n = size(A, d)
    ibuf = Vector{T}(undef, n)
    obuf = Vector{T}(undef, n)
    if @generated
        quote
            Base.Cartesian.@nif $N dim -> (d == dim) dim -> (fft_along_dim!(A, ibuf, obuf, cg, dir, Val(dim)))
        end
    else
        fft_along_dim!(A, ibuf, obuf, cg, dir, Val(d))
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
        _rfft_pencil!(y, x, p)
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
        _brfft_pencil!(y, x, p)
    else
        _re_along_dim!(_brfft_pencil!, y, x, p, d1)
    end
    return y
end

#### 2D plan
# The real transform is taken along the first region dimension, then a complex
# transform along the second (forward), or the reverse (backward).
##### Forward
function LinearAlgebra.mul!(y::AbstractArray{T,N}, p::FFTAPlan_re{T,2}, x::AbstractArray{<:Real,N}) where {T<:Complex,N}
    if p.dir !== FFT_FORWARD
        throw(ArgumentError("only FFT_FORWARD supported for real arrays"))
    end
    Base.require_one_based_indexing(x, y)
    d1, d2 = p.region
    _check_re_dims(y, p, x, d1, true)
    if size(x, d2) != size(p, 2)
        throw(DimensionMismatch("real 2D plan has size $(size(p)). Transform dimensions of input array are $((size(x, d1), size(x, d2))) but should be $(size(p))"))
    end
    _re_along_dim!(_rfft_pencil!, y, x, p, d1)
    _cx_along_dim!(y, p.callgraph[2], FFT_FORWARD, d2)
    return y
end

##### Backward
function LinearAlgebra.mul!(y::AbstractArray{<:Real,N}, p::FFTAPlan_re{T,2}, x::AbstractArray{T,N}) where {T<:Complex,N}
    if p.dir !== FFT_BACKWARD
        throw(ArgumentError("only FFT_BACKWARD supported for complex arrays"))
    end
    Base.require_one_based_indexing(x, y)
    d1, d2 = p.region
    _check_re_dims(y, p, x, d1, false)
    if size(x, d2) != size(p, 2)
        throw(DimensionMismatch("real 2D plan has size $(size(p)). Transform dimensions of input array are $((size(x, d1), size(x, d2))) but should be $((size(p, 1) ÷ 2 + 1, size(p, 2)))"))
    end
    tmp = copy(x)   # the complex pass must not modify the input
    _cx_along_dim!(tmp, p.callgraph[2], FFT_BACKWARD, d2)
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
    q = FFTAPlan_cx{T,N,R}(cg, p.region, dir)
    AbstractFFTs.ScaledPlan(q, _normalization(p))
end

function AbstractFFTs.plan_inv(p::FFTAPlan_re{T,N,R,S}) where {T,N,R,S}
    dir = p.dir === FFT_FORWARD ? FFT_BACKWARD : FFT_FORWARD
    cutoff = p.callgraph[1].BLUESTEIN_CUTOFF
    n = p.flen
    nn = iseven(n) ? n >> 1 : n
    cg = (CallGraph{T}(nn, cutoff), ntuple(i -> CallGraph{T}(size(p, i + 1), cutoff), Val(N - 1))...)
    q = FFTAPlan_re{T,N,R,S}(cg, p.region, dir, n, Vector{T}(undef, _re_buflen(n, dir)), Vector{S}(undef, n), Vector{T}(undef, n ÷ 2 + 1))
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
