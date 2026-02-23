# Plans

abstract type FFTAPlan{T,N} <: AbstractFFTs.Plan{T} end

struct FFTAInvPlan{T,N} <: FFTAPlan{T,N} end

struct FFTAPlan_cx{T,N,R<:Union{Int,AbstractVector{Int}}} <: FFTAPlan{T,N}
    callgraph::NTuple{N,CallGraph{T}}
    region::R
    dir::Direction
    pinv::FFTAInvPlan{T,N}
end
function FFTAPlan_cx{T,N}(
    cg::NTuple{N,CallGraph{T}}, r::R,
    dir::Direction, pinv::FFTAInvPlan{T,N}
) where {T,N,R<:Union{Int,AbstractVector{Int}}}
    FFTAPlan_cx{T,N,R}(cg, r, dir, pinv)
end

struct FFTAPlan_re{T,N,R<:Union{Int,AbstractVector{Int}}} <: FFTAPlan{T,N}
    callgraph::NTuple{N,CallGraph{T}}
    region::R
    dir::Direction
    pinv::FFTAInvPlan{T,N}
    flen::Int
end
function FFTAPlan_re{T,N}(
    cg::NTuple{N,CallGraph{T}}, r::R,
    dir::Direction, pinv::FFTAInvPlan{T,N}, flen::Int
) where {T,N,R<:Union{Int,AbstractVector{Int}}}
    FFTAPlan_re{T,N,R}(cg, r, dir, pinv, flen)
end

Base.size(p::FFTAPlan_cx, i::Int) = i <= length(p.callgraph) ? first(p.callgraph[i].nodes).sz : 1
function Base.size(p::FFTAPlan_re{<:Any,1}, i::Int)
    if i == p.region[]
        p.flen
    elseif i <= length(p.callgraph)
        first(p.callgraph[i].nodes).sz
    else
        1
    end
end
function Base.size(p::FFTAPlan_re{<:Any,2}, i::Int)
    if i == 1
        return p.flen
    elseif i == 2
        first(p.callgraph[2].nodes).sz
    else
        1
    end
end
Base.size(p::FFTAPlan{<:Any,N}) where N = ntuple(Base.Fix1(size, p), Val{N}())

Base.complex(p::FFTAPlan_re{T,N,R}) where {T,N,R} = FFTAPlan_cx{T,N,R}(p.callgraph, p.region, p.dir, p.pinv)

AbstractFFTs.plan_fft(x::AbstractArray{T,N}, region::R; kwargs...) where {T<:Complex,N,R} =
    _plan_fft(x, region, FFT_FORWARD; kwargs...)

AbstractFFTs.plan_bfft(x::AbstractArray{T,N}, region::R; kwargs...) where {T<:Complex,N,R} =
    _plan_fft(x, region, FFT_BACKWARD; kwargs...)

function _plan_fft(x::AbstractArray{T,N}, region::R, dir::Direction; BLUESTEIN_CUTOFF=DEFAULT_BLUESTEIN_CUTOFF, kwargs...) where {T<:Complex,N,R}
    FFTN = length(region)
    if FFTN == 1
        R1 = Int(region[])
        g = CallGraph{T}(size(x, R1), BLUESTEIN_CUTOFF)
        pinv = FFTAInvPlan{T,1}()
        return FFTAPlan_cx{T,1,Int}((g,), R1, dir, pinv)
    elseif FFTN == 2
        sort!(region)
        g1 = CallGraph{T}(size(x, region[1]), BLUESTEIN_CUTOFF)
        g2 = CallGraph{T}(size(x, region[2]), BLUESTEIN_CUTOFF)
        pinv = FFTAInvPlan{T,2}()
        return FFTAPlan_cx{T,2,R}((g1, g2), region, dir, pinv)
    else
        sort!(region)
        return FFTAPlan_cx{T,FFTN,R}(
            ntuple(i -> CallGraph{T}(size(x, region[i]), BLUESTEIN_CUTOFF), Val(FFTN)),
            region, dir, FFTAInvPlan{T,FFTN}()
        )
    end
end

function AbstractFFTs.plan_rfft(x::AbstractArray{T,N}, region::R; BLUESTEIN_CUTOFF=DEFAULT_BLUESTEIN_CUTOFF, kwargs...) where {T<:Real,N,R}
    FFTN = length(region)
    if FFTN == 1
        R1 = Int(region[])
        n = size(x, R1)
        # For even length problems, we solve the real problem with
        # two n/2 complex FFTs followed by a butterfly. For odd size
        # problems, we just solve the problem as a single complex
        nn = iseven(n) ? n >> 1 : n
        g = CallGraph{Complex{T}}(nn, BLUESTEIN_CUTOFF)
        pinv = FFTAInvPlan{Complex{T},1}()
        return FFTAPlan_re{Complex{T},1,Int}((g,), R1, FFT_FORWARD, pinv, n)
    elseif FFTN == 2
        sort!(region)
        g1 = CallGraph{Complex{T}}(size(x, region[1]), BLUESTEIN_CUTOFF)
        g2 = CallGraph{Complex{T}}(size(x, region[2]), BLUESTEIN_CUTOFF)
        pinv = FFTAInvPlan{Complex{T},2}()
        return FFTAPlan_re{Complex{T},2,R}((g1, g2), region, FFT_FORWARD, pinv, size(x, region[1]))
    else
        throw(ArgumentError("only supports 1D and 2D FFTs"))
    end
end

function AbstractFFTs.plan_brfft(x::AbstractArray{T,N}, len, region::R; BLUESTEIN_CUTOFF=DEFAULT_BLUESTEIN_CUTOFF, kwargs...) where {T,N,R}
    FFTN = length(region)
    if FFTN == 1
        # For even length problems, we solve the real problem with
        # two n/2 complex FFTs followed by a butterfly. For odd size
        # problems, we just solve the problem as a single complex
        R1 = Int(region[])
        nn = iseven(len) ? len >> 1 : len
        g = CallGraph{T}(nn, BLUESTEIN_CUTOFF)
        pinv = FFTAInvPlan{T,1}()
        return FFTAPlan_re{T,1,Int}((g,), R1, FFT_BACKWARD, pinv, len)
    elseif FFTN == 2
        sort!(region)
        g1 = CallGraph{T}(len, BLUESTEIN_CUTOFF)
        g2 = CallGraph{T}(size(x, region[2]), BLUESTEIN_CUTOFF)
        pinv = FFTAInvPlan{T,2}()
        return FFTAPlan_re{T,2,R}((g1, g2), region, FFT_BACKWARD, pinv, len)
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
    fft!(y, x, 1, 1, p.dir, p.callgraph[1][1].type, p.callgraph[1], 1)
    return y
end

#### 1D plan ND array
function LinearAlgebra.mul!(y::AbstractArray{U,N}, p::FFTAPlan_cx{T,1}, x::AbstractArray{T,N}) where {T,U,N}
    Base.require_one_based_indexing(x, y)
    if axes(x) != axes(y)
        throw(DimensionMismatch("input array has axes $(axes(x)), but output array has axes $(axes(y))"))
    end
    if size(p, 1) != size(x, p.region[])
        throw(DimensionMismatch("plan has size $(size(p, 1)), but input array has size $(size(x, p.region[])) along region $(p.region[])"))
    end
    Rpre  = CartesianIndices(size(x)[1:p.region[]-1])
    Rpost = CartesianIndices(size(x)[p.region[]+1:end])
    _mul_loop!(y, x, Rpre, Rpost, p)
    return y
end

function _mul_loop!(
    y::AbstractArray{U,N},
    x::AbstractArray{T,N},
    Rpre::CartesianIndices,
    Rpost::CartesianIndices,
    p::FFTAPlan_cx{T,1}) where {T,U,N}
    for Ipost in Rpost, Ipre in Rpre
        @views fft!(y[Ipre,:,Ipost], x[Ipre,:,Ipost], 1, 1, p.dir, p.callgraph[1][1].type, p.callgraph[1], 1)
    end
end

#### ND plan ND array
@generated function LinearAlgebra.mul!(
    out::AbstractArray{U,N},
    p::FFTAPlan_cx{T,N},
    X::AbstractArray{T,N}
) where {T,U,N}

    quote
        Base.require_one_based_indexing(out, X)
        if size(out) != size(X)
            throw(DimensionMismatch("input array has axes $(axes(X)), but output array has axes $(axes(out))"))
        elseif size(p) != size(X)
            throw(DimensionMismatch("plan has size $(size(p)), but input array has size $(size(X))"))
        elseif !(p.region == N || p.region == 1:N)
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

        Base.Cartesian.@nexprs $N dim -> begin
            n = size(out, dim)
            resize!(obuf, n)
            resize!(ibuf, n)
            cg = p.callgraph[dim]

            Rpre_{dim}  = CartesianIndices(sz[1:dim-1])
            Rpost_{dim} = CartesianIndices(sz[dim+1:N])

            fft_along_dim!(out, ibuf, obuf, cg, dir, Rpre_{dim}, Rpost_{dim})
        end

        return out
    end
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
    elseif M > N || p.region[1] < 1 || p.region[end] > N
        throw(DimensionMismatch("Plan region is outside array dimensions."))
    end

    sz = size(X)
    max_sz = maximum(Base.Fix1(size, out), p.region)
    obuf = Vector{T}(undef, max_sz)
    ibuf = Vector{T}(undef, max_sz)
    sizehint!(obuf, max_sz) # not guaranteed but hopefully prevents allocations
    sizehint!(ibuf, max_sz)
    dir = p.dir

    copyto!(out, X) # operate in-place on output array

    # don't use generated functions because this cannot be type-stable anyway
    for dim in 1:M
        pdim = p.region[dim]
        n = size(out, pdim)
        resize!(obuf, n)
        resize!(ibuf, n)
        cg = p.callgraph[dim]

        Rpre  = CartesianIndices(sz[1:pdim-1])
        Rpost = CartesianIndices(sz[pdim+1:N])

        fft_along_dim!(out, ibuf, obuf, cg, dir, Rpre, Rpost)
    end

    return out
end

function fft_along_dim!(
    A::AbstractArray,
    ibuf::Vector{T}, obuf::Vector{T},
    cg::CallGraph{T}, d::Direction,
    Rpre::CartesianIndices{M}, Rpost::CartesianIndices
) where {T <: Complex{<:AbstractFloat}, M}

    t = cg[1].type
    dim = M + 1
    cols = eachindex(axes(A, dim), ibuf, obuf)

    for Ipost in Rpost, Ipre in Rpre
        for j in cols
            ibuf[j] = A[Ipre, j, Ipost]
        end
        fft!(obuf, ibuf, 1, 1, d, t, cg, 1)
        for j in cols
            A[Ipre, j, Ipost] = obuf[j]
        end
    end
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
# By converting the problem to complex and back to real
#### 1D plan 1D array
##### Forward
function Base.:*(p::FFTAPlan_re{Complex{T},1}, x::AbstractVector{T}) where {T<:Real}
    Base.require_one_based_indexing(x)
    if p.dir === FFT_FORWARD
        n = p.flen
        if iseven(n)
            # For problems of even size, we solve the rfft problem by splitting the
            # problem into the even and odd part and solving the simultanously as
            # a single (complex) fft of half the size, see equations (6)-(8) of
            # Sorensen, H. V., D. Jones, Michael Heideman, and C. Burrus.
            # "Real-valued fast Fourier transform algorithms."
            # IEEE Transactions on acoustics, speech, and signal processing 35, no. 6 (2003): 849-863.
            if x isa Vector && isbitstype(T)
                # For a vector of bits, we can just reintepret the bits to get the
                # approciate representation of even (zero based) elements as the real
                # part and the odd as the complex part
                x_c = reinterpret(Complex{T}, x)
            else
                # for non-bits, we'd have to copy to a new array
                x_c = complex.(view(x, 1:2:n), view(x, 2:2:n))
            end

            m = n >> 1
            # Allocate complex result vector of half the input size plus one
            y = similar(x_c, m + 1)
            # Solve the complex fft of half the size
            LinearAlgebra.mul!(view(y, 1:m), complex(p), x_c)

            # The w stored in the plan is for m, not n, so probably cheapest to
            # just recompute it instead of taking a square root
            wj = w = cispi(-T(2) / n)

            # Construct the result by first constructing the elements of the
            # real and imaginary part, followed by the usual radix-2 assembly,
            # see eq (9)
            y1     = y[1]
            y[1]   = real(y1) + imag(y1)
            y[end] = real(y1) - imag(y1)

            @inbounds for j in 2:((m >> 1) + 1)
                yj  = y[j]
                ymj = y[m-j+2]
                XX = T(0.5) * ( yj + conj(ymj))
                XY = T(0.5) * (-yj + conj(ymj)) * im
                y[j]     =      XX + wj * XY
                y[m-j+2] = conj(XX - wj * XY)
                wj *= w
            end
            return y
        else
            # when the problem cannot be split in two equal size chunks we
            # convert the problem to a complex fft and truncate the redundant
            # part of the result vector
            x_c = similar(x, Complex{T})
            y = similar(x_c)
            copyto!(x_c, x)
            LinearAlgebra.mul!(y, complex(p), x_c)
            return y[1:end÷2+1]
        end
    end
    throw(ArgumentError("only FFT_FORWARD supported for real vectors"))
end

##### Backward
function Base.:*(p::FFTAPlan_re{T,1}, x::AbstractVector{T}) where {T<:Complex}
    Base.require_one_based_indexing(x)
    if p.dir === FFT_BACKWARD
        n = p.flen
        # See explantion of this approach in the method for the FORWARD transform
        if iseven(n)
            m = n >> 1
            wj = w = cispi(T(2) / n)
            x_tmp = similar(x, length(x) - 1)
            x_tmp[1] = complex(
                (real(x[1]) + real(x[end])),
                (real(x[1]) - real(x[end]))
            )
            for j in 2:((m >> 1) + 1)
                XX =       x[j] + conj(x[m-j+2])
                XY = wj * (x[j] - conj(x[m-j+2]))
                x_tmp[j]     =      XX + im * XY
                x_tmp[m-j+2] = conj(XX - im * XY)
                wj *= w
            end
            y_c = complex(p) * x_tmp
            if isbitstype(T)
                return copy(reinterpret(real(T), y_c))
            else
                return mapreduce(t -> [real(t); imag(t)], vcat, y_c)
            end
        else
            x_tmp = similar(x, n)
            x_tmp[1:end÷2+1] .= x
            x_tmp[end÷2+2:end] .= iseven(n) ? conj.(x[end-1:-1:2]) : conj.(x[end:-1:2])
            y = similar(x_tmp)
            LinearAlgebra.mul!(y, complex(p), x_tmp)
            return real(y)
        end
    end
    throw(ArgumentError("only FFT_BACKWARD supported for complex vectors"))
end

#### 1D plan ND array
##### Forward
function Base.:*(p::FFTAPlan_re{Complex{T},1}, x::AbstractArray{T,N}) where {T<:Real,N}
    Base.require_one_based_indexing(x)
    if p.dir === FFT_FORWARD
        return mapslices(Base.Fix1(*, p), x; dims=p.region[1])
    end
    throw(ArgumentError("only FFT_FORWARD supported for real arrays"))
end

##### Backward
function Base.:*(p::FFTAPlan_re{T,1}, x::AbstractArray{T,N}) where {T<:Complex,N}
    Base.require_one_based_indexing(x)
    if p.flen ÷ 2 + 1 != size(x, p.region[])
        throw(DimensionMismatch("real 1D plan has size $(p.flen). Dimension of input array along region $(p.region[]) should have size $(size(p, p.region[]) ÷ 2 + 1), but has size $(size(x, p.region[]))"))
    end
    if p.dir === FFT_BACKWARD
        return mapslices(Base.Fix1(*, p), x; dims=p.region[1])
    end
    throw(ArgumentError("only FFT_BACKWARD supported for complex arrays"))
end

#### 2D plan ND array
##### Forward
function Base.:*(p::FFTAPlan_re{Complex{T},2}, x::AbstractArray{T,N}) where {T<:Real,N}
    Base.require_one_based_indexing(x)
    if p.dir === FFT_FORWARD
        half_1 = 1:(p.flen÷2+1)
        x_c = similar(x, Complex{T})
        copy!(x_c, x)
        y = similar(x_c)
        LinearAlgebra.mul!(y, complex(p), x_c)
        return copy(selectdim(y, p.region[1], half_1))
    end
    throw(ArgumentError("only FFT_FORWARD supported for real arrays"))
end

##### Backward
function Base.:*(p::FFTAPlan_re{T,2}, x::AbstractArray{T,N}) where {T<:Complex,N}
    Base.require_one_based_indexing(x)
    if size(p, 1) ÷ 2 + 1 != size(x, p.region[1])
        throw(DimensionMismatch("real 2D plan has size $(size(p)). First transform dimension of input array should have size ($(size(p, 1) ÷ 2 + 1)), but has size $(size(x, p.region[1]))"))
    end
    if p.dir === FFT_BACKWARD
        res_size = ntuple(i -> ifelse(i == p.region[1], p.flen, size(x, i)), Val(N))
        # for the inverse transformation we have to reconstruct the full array
        half_1 = 1:(p.flen ÷ 2 + 1)
        half_2 = half_1[end]+1:p.flen
        x_full = similar(x, res_size)
        # use first half as is
        copy!(selectdim(x_full, p.region[1], half_1), x)

        # the second half in the first transform dimension is reversed and conjugated
        x_half_2 = selectdim(x_full, p.region[1], half_2) # view to the second half of x
        start_reverse = size(x, p.region[1]) - iseven(p.flen)

        map!(conj, x_half_2, selectdim(x, p.region[1], start_reverse:-1:2))
        # for the 2D transform we have to reverse index 2:end of the same block in the second transform dimension as well
        reverse!(selectdim(x_half_2, p.region[2], 2:size(x, p.region[2])), dims=p.region[2])

        y = similar(x_full)
        LinearAlgebra.mul!(y, complex(p), x_full)
        return real(y)
    end
    throw(ArgumentError("only FFT_BACKWARD supported for complex arrays"))
end
