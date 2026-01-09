# Plans

abstract type FFTAPlan{T,N} <: AbstractFFTs.Plan{T} end

struct FFTAInvPlan{T,N} <: FFTAPlan{T,N} end

struct FFTAPlan_cx{T,N} <: FFTAPlan{T,N}
    callgraph::NTuple{N, CallGraph{T}}
    region::Union{Int,AbstractVector{<:Int}}
    dir::Direction
    pinv::FFTAInvPlan{T}
end

struct FFTAPlan_re{T,N} <: FFTAPlan{T,N}
    callgraph::NTuple{N, CallGraph{T}}
    region::Union{Int,AbstractVector{<:Int}}
    dir::Direction
    pinv::FFTAInvPlan{T}
    flen::Int
end

Base.size(p::FFTAPlan, i::Int) = i <= length(p.callgraph) ? first(p.callgraph[i].nodes).sz : 1
Base.size(p::FFTAPlan{<:Any,N}) where N = ntuple(Base.Fix1(size, p), Val{N}())

Base.complex(p::FFTAPlan_re{T,N}) where {T,N} = FFTAPlan_cx{T,N}(p.callgraph, p.region, p.dir, p.pinv)

function AbstractFFTs.plan_fft(x::AbstractArray{T,N}, region; kwargs...)::FFTAPlan_cx{T} where {T <: Complex, N}
    FFTN = length(region)
    if FFTN == 1
        g = CallGraph{T}(size(x,region[]))
        pinv = FFTAInvPlan{T,FFTN}()
        return FFTAPlan_cx{T,FFTN}((g,), region, FFT_FORWARD, pinv)
    elseif FFTN == 2
        sort!(region)
        g1 = CallGraph{T}(size(x,region[1]))
        g2 = CallGraph{T}(size(x,region[2]))
        pinv = FFTAInvPlan{T,FFTN}()
        return FFTAPlan_cx{T,FFTN}((g1,g2), region, FFT_FORWARD, pinv)
    else
        throw(ArgumentError("only supports 1D and 2D FFTs"))
    end
end

function AbstractFFTs.plan_bfft(x::AbstractArray{T,N}, region; kwargs...)::FFTAPlan_cx{T} where {T <: Complex,N}
    FFTN = length(region)
    if FFTN == 1
        g = CallGraph{T}(size(x,region[]))
        pinv = FFTAInvPlan{T,FFTN}()
        return FFTAPlan_cx{T,FFTN}((g,), region, FFT_BACKWARD, pinv)
    elseif FFTN == 2
        sort!(region)
        g1 = CallGraph{T}(size(x,region[1]))
        g2 = CallGraph{T}(size(x,region[2]))
        pinv = FFTAInvPlan{T,FFTN}()
        return FFTAPlan_cx{T,FFTN}((g1,g2), region, FFT_BACKWARD, pinv)
    else
        throw(ArgumentError("only supports 1D and 2D FFTs"))
    end
end

function AbstractFFTs.plan_rfft(x::AbstractArray{T,N}, region; kwargs...)::FFTAPlan_re{Complex{T}} where {T <: Real,N}
    FFTN = length(region)
    if FFTN == 1
        g = CallGraph{Complex{T}}(size(x,region[]))
        pinv = FFTAInvPlan{Complex{T},FFTN}()
        return FFTAPlan_re{Complex{T},FFTN}(tuple(g), region, FFT_FORWARD, pinv, size(x,region[]))
    elseif FFTN == 2
        if N !== 2
            throw(ArgumentError("2D real FFT only supported for 2D arrays"))
        end
        sort!(region)
        g1 = CallGraph{Complex{T}}(size(x,region[1]))
        g2 = CallGraph{Complex{T}}(size(x,region[2]))
        pinv = FFTAInvPlan{Complex{T},FFTN}()
        return FFTAPlan_re{Complex{T},FFTN}(tuple(g1,g2), region, FFT_FORWARD, pinv, size(x,region[1]))
    else
        throw(ArgumentError("only supports 1D and 2D FFTs"))
    end
end

function AbstractFFTs.plan_brfft(x::AbstractArray{T,N}, len, region; kwargs...)::FFTAPlan_re{T} where {T,N}
    FFTN = length(region)
    if FFTN == 1
        g = CallGraph{T}(len)
        pinv = FFTAInvPlan{T,FFTN}()
        return FFTAPlan_re{T,FFTN}((g,), region, FFT_BACKWARD, pinv, len)
    elseif FFTN == 2
        if N !== 2
            throw(ArgumentError("2D real FFT only supported for 2D arrays"))
        end
        sort!(region)
        g1 = CallGraph{T}(len)
        g2 = CallGraph{T}(size(x,region[2]))
        pinv = FFTAInvPlan{T,FFTN}()
        return FFTAPlan_re{T,FFTN}((g1,g2), region, FFT_BACKWARD, pinv, len)
    else
        throw(ArgumentError("only supports 1D and 2D FFTs"))
    end
end


# Multiplication
## mul!
### Complex
function LinearAlgebra.mul!(y::AbstractVector{U}, p::FFTAPlan_cx{T,1}, x::AbstractVector{T}) where {T,U}
    if axes(x) != axes(y)
        throw(DimensionMismatch("input array has axes $(axes(x)), but output array has axes $(axes(y))"))
    end
    if size(p) != size(x)
        throw(DimensionMismatch("plan has axes $(size(p)), but input array has axes $(size(x))"))
    end
    fft!(y, x, 1, 1, p.dir, p.callgraph[1][1].type, p.callgraph[1], 1)
end

function LinearAlgebra.mul!(y::AbstractArray{U,N}, p::FFTAPlan_cx{T,1}, x::AbstractArray{T,N}) where {T,U,N}
    Base.require_one_based_indexing(x)
    if axes(x) != axes(y)
        throw(DimensionMismatch("input array has axes $(axes(x)), but output array has axes $(axes(y))"))
    end
    if size(p, 1) != size(x, p.region[])
        throw(DimensionMismatch("plan has size $(size(p, 1)), but input array has size $(size(x, p.region[])) along region $(p.region[])"))
    end
    Rpre = CartesianIndices(size(x)[1:p.region[]-1])
    Rpost = CartesianIndices(size(x)[p.region[]+1:end])
    for Ipre in Rpre
        for Ipost in Rpost
            @views fft!(y[Ipre,:,Ipost], x[Ipre,:,Ipost], 1, 1, p.dir, p.callgraph[1][1].type, p.callgraph[1], 1)
        end
    end
end

function LinearAlgebra.mul!(y::AbstractArray{U,N}, p::FFTAPlan_cx{T,2}, x::AbstractArray{T,N}) where {T,U,N}
    Base.require_one_based_indexing(x)
    if axes(x) != axes(y)
        throw(DimensionMismatch("input array has axes $(axes(x)), but output array has axes $(axes(y))"))
    end
    if N < 2
        throw(DimensionMismatch("array dimension $N cannot be smaller than the plan size 2"))
    end
    if size(p) != (size(x, p.region[1]), size(x, p.region[2]))
        throw(DimensionMismatch("plan has size $(size(p)), but input array has size $((size(x, p.region[1]), size(x, p.region[2]))) along regions $(p.region)"))
    end
    R1 = CartesianIndices(size(x)[1:p.region[1]-1])
    R2 = CartesianIndices(size(x)[p.region[1]+1:p.region[2]-1])
    R3 = CartesianIndices(size(x)[p.region[2]+1:end])
    y_tmp = similar(y, axes(y)[p.region])
    rows,cols = size(x)[p.region]
    # Introduce function barrier here since the variables used in the loop ranges aren't inferred. This
    # is partly because the region field of the plan is abstractly typed but even if that wasn't the case,
    # it might be a bit tricky to construct the Rxs in an inferred way.
    _mul_loop(y_tmp, y, x, p, R1, R2, R3, rows, cols)
end

function _mul_loop(
    y_tmp::AbstractArray,
    y::AbstractArray,
    x::AbstractArray,
    p::FFTAPlan,
    R1::CartesianIndices,
    R2::CartesianIndices,
    R3::CartesianIndices,
    rows::Int,
    cols::Int
)
    for I1 in R1
        for I2 in R2
            for I3 in R3
                for k in 1:cols
                    @views fft!(y_tmp[:,k],  x[I1,:,I2,k,I3], 1, 1, p.dir, p.callgraph[1][1].type, p.callgraph[1], 1)
                end

                for k in 1:rows
                    @views fft!(y[I1,k,I2,:,I3], y_tmp[k,:], 1, 1, p.dir, p.callgraph[2][1].type, p.callgraph[2], 1)
                end
            end
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

function Base.:*(p::FFTAPlan_cx{T,N1}, x::AbstractArray{T,N2}) where {T<:Complex, N1, N2}
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
        x_c = similar(x, Complex{T})
        copy!(x_c, x)
        y = similar(x_c)
        LinearAlgebra.mul!(y, complex(p), x_c)
        return y[1:end÷2 + 1]
    end
    throw(ArgumentError("only FFT_FORWARD supported for real vectors"))
end

##### Backward
function Base.:*(p::FFTAPlan_re{T,1}, x::AbstractVector{T}) where {T<:Complex}
    Base.require_one_based_indexing(x)
    if p.dir === FFT_BACKWARD
        x_tmp = similar(x, p.flen)
        x_tmp[1:end÷2 + 1] .= x
        x_tmp[end÷2 + 2:end] .= iseven(p.flen) ? conj.(x[end-1:-1:2]) : conj.(x[end:-1:2])
        y = similar(x_tmp)
        LinearAlgebra.mul!(y, complex(p), x_tmp)
        return real(y)
    end
    throw(ArgumentError("only FFT_BACKWARD supported for complex vectors"))
end

#### 1D plan ND array
##### Forward
function Base.:*(p::FFTAPlan_re{Complex{T},1}, x::AbstractArray{T,N}) where {T<:Real, N}
    Base.require_one_based_indexing(x)
    if p.dir === FFT_FORWARD
        half_1 = 1:(p.flen ÷ 2 + 1)
        x_c = similar(x, Complex{T})
        copy!(x_c, x)
        y = similar(x_c)
        LinearAlgebra.mul!(y, complex(p), x_c)
        return copy(selectdim(y, p.region[1], half_1))
    end
    throw(ArgumentError("only FFT_FORWARD supported for real arrays"))
end

##### Backward
function Base.:*(p::FFTAPlan_re{T,1}, x::AbstractArray{T,N}) where {T<:Complex, N}
    Base.require_one_based_indexing(x)
    if p.flen ÷ 2 + 1 != size(x, p.region[])
        throw(DimensionMismatch("real 1D plan has size $(p.flen). Dimension of input array along region $(p.region[]) should have size $(size(p, p.region[]) ÷ 2 + 1), but has size $(size(x, p.region[]))"))
    end
    if p.dir === FFT_BACKWARD
        # # for the inverse transformation we have to reconstruct the full array
        half_1 = 1:(p.flen ÷ 2 + 1)
        half_2 = half_1[end]+1:p.flen
        res_size = ntuple(i->ifelse(i==p.region[1], p.flen, size(x,i)), ndims(x))
        # for the inverse transformation we have to reconstruct the full array
        x_full = similar(x, res_size)
        # use first half as is
        copy!(selectdim(x_full, p.region[1], half_1), x)
        start_reverse = size(x, p.region[1]) - iseven(p.flen)
        half_reverse = (start_reverse:-1:2)
        # the second half is reversed and conjugated
        map!(conj, selectdim(x_full, p.region[1], half_2), selectdim(x, p.region[1], half_reverse))
        y = similar(x_full)
        LinearAlgebra.mul!(y, complex(p), x_full)
        return real(y)
    end
    throw(ArgumentError("only FFT_BACKWARD supported for complex arrays"))
end

#### 2D plan 2D array
##### Forward
function Base.:*(p::FFTAPlan_re{Complex{T},2}, x::AbstractArray{T,2}) where {T<:Real}
    Base.require_one_based_indexing(x)
    if p.dir === FFT_FORWARD
        half_1 = 1:(p.flen ÷ 2 + 1)
        x_c = similar(x, Complex{T})
        copy!(x_c, x)
        y = similar(x_c)
        LinearAlgebra.mul!(y, complex(p), x_c)
        return y[half_1, :]
    end
    throw(ArgumentError("only FFT_FORWARD supported for real arrays"))
end

##### Backward
function Base.:*(p::FFTAPlan_re{T,2}, x::AbstractArray{T,2}) where {T<:Complex}
    Base.require_one_based_indexing(x)
    if size(p, 1) ÷ 2 + 1 != size(x, 1)
        throw(DimensionMismatch("real 2D plan has size $(size(p)). First dimension of input array should have size ($(size(p, 1) ÷ 2 + 1)), but has size $(size(x, 1))"))
    end
    if p.dir === FFT_BACKWARD
        # for the inverse transformation we have to reconstruct the full array
        m, n = size(x)
        half_1 = 1:(p.flen ÷ 2 + 1)
        half_2 = half_1[end]+1:p.flen
        x_full = similar(x, p.flen, n)
        x_full[1:m, :] = x
        start_reverse = m - iseven(p.flen)
        map!(conj, view(x_full, (m + 1):p.flen, 1), view(x, start_reverse:-1:2, 1))
        map!(conj, view(x_full, half_2, 2:n), view(x, start_reverse:-1:2, n:-1:2))

        y = similar(x_full)
        LinearAlgebra.mul!(y, complex(p), x_full)
        return real(y)
    end
    throw(ArgumentError("only FFT_BACKWARD supported for complex arrays"))
end
