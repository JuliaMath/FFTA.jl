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

function AbstractFFTs.plan_fft(x::AbstractArray{T}, region; kwargs...)::FFTAPlan_cx{T} where {T <: Complex}
    N = length(region)
    @assert N <= 2 "Only supports vectors and matrices"
    if N == 1
        g = CallGraph{T}(size(x,region[]))
        pinv = FFTAInvPlan{T,N}()
        return FFTAPlan_cx{T,N}((g,), region, FFT_FORWARD, pinv)
    else
        sort!(region)
        g1 = CallGraph{T}(size(x,region[1]))
        g2 = CallGraph{T}(size(x,region[2]))
        pinv = FFTAInvPlan{T,N}()
        return FFTAPlan_cx{T,N}((g1,g2), region, FFT_FORWARD, pinv)
    end
end

function AbstractFFTs.plan_bfft(x::AbstractArray{T}, region; kwargs...)::FFTAPlan_cx{T} where {T <: Complex}
    N = length(region)
    @assert N <= 2 "Only supports vectors and matrices"
    if N == 1
        g = CallGraph{T}(size(x,region[]))
        pinv = FFTAInvPlan{T,N}()
        return FFTAPlan_cx{T,N}((g,), region, FFT_BACKWARD, pinv)
    else
        sort!(region)
        g1 = CallGraph{T}(size(x,region[1]))
        g2 = CallGraph{T}(size(x,region[2]))
        pinv = FFTAInvPlan{T,N}()
        return FFTAPlan_cx{T,N}((g1,g2), region, FFT_BACKWARD, pinv)
    end
end

function AbstractFFTs.plan_rfft(x::AbstractArray{T}, region; kwargs...)::FFTAPlan_re{Complex{T}} where {T <: Real}
    N = length(region)
    @assert N <= 2 "Only supports vectors and matrices"
    if N == 1
        g = CallGraph{Complex{T}}(size(x,region[]))
        pinv = FFTAInvPlan{Complex{T},N}()
        return FFTAPlan_re{Complex{T},N}(tuple(g), region, FFT_FORWARD, pinv, size(x,region[]))
    else
        sort!(region)
        g1 = CallGraph{Complex{T}}(size(x,region[1]))
        g2 = CallGraph{Complex{T}}(size(x,region[2]))
        pinv = FFTAInvPlan{Complex{T},N}()
        return FFTAPlan_re{Complex{T},N}(tuple(g1,g2), region, FFT_FORWARD, pinv, size(x,region[1]))
    end
end

function AbstractFFTs.plan_brfft(x::AbstractArray{T}, len, region; kwargs...)::FFTAPlan_re{T} where {T}
    N = length(region)
    @assert N <= 2 "Only supports vectors and matrices"
    if N == 1
        g = CallGraph{T}(len)
        pinv = FFTAInvPlan{T,N}()
        return FFTAPlan_re{T,N}((g,), region, FFT_BACKWARD, pinv, len)
    else
        sort!(region)
        g1 = CallGraph{T}(len)
        g2 = CallGraph{T}(size(x,region[2]))
        pinv = FFTAInvPlan{T,N}()
        return FFTAPlan_re{T,N}((g1,g2), region, FFT_BACKWARD, pinv, len)
    end
end

function AbstractFFTs.plan_bfft(p::FFTAPlan_cx{T,N}) where {T,N}
    return FFTAPlan_cx{T,N}(p.callgraph, p.region, -p.dir, p.pinv)
end

function AbstractFFTs.plan_brfft(p::FFTAPlan_re{T,N}) where {T,N}
    return FFTAPlan_re{T,N}(p.callgraph, p.region, -p.dir, p.pinv, p.flen)
end

function LinearAlgebra.mul!(y::AbstractVector{U}, p::FFTAPlan{T,1}, x::AbstractVector{T}) where {T,U}
    fft!(y, x, 1, 1, p.dir, p.callgraph[1][1].type, p.callgraph[1], 1)
end

function LinearAlgebra.mul!(y::AbstractArray{U,N}, p::FFTAPlan{T,1}, x::AbstractArray{T,N}) where {T,U,N}
    Rpre = CartesianIndices(size(x)[1:p.region-1])
    Rpost = CartesianIndices(size(x)[p.region+1:end])
    for Ipre in Rpre
        for Ipost in Rpost
            @views fft!(y[Ipre,:,Ipost], x[Ipre,:,Ipost], 1, 1, p.dir, p.callgraph[1][1].type, p.callgraph[1], 1)
        end
    end
end

function LinearAlgebra.mul!(y::AbstractArray{U,N}, p::FFTAPlan{T,2}, x::AbstractArray{T,N}) where {T,U,N}
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

function Base.:*(p::FFTAPlan{T,1}, x::AbstractVector{T}) where {T<:Union{Real,Complex}}
    y = similar(x, T <: Real ? Complex{T} : T)
    LinearAlgebra.mul!(y, p, x)
    y
end

function Base.:*(p::FFTAPlan{T,N1}, x::AbstractArray{T,N2}) where {T<:Union{Real, Complex}, N1, N2}
    y = similar(x, T <: Real ? Complex{T} : T)
    LinearAlgebra.mul!(y, p, x)
    y
end

function Base.:*(p::FFTAPlan_re{T,1}, x::AbstractVector{T}) where {T<:Union{Real, Complex}}
    if p.dir == FFT_FORWARD
        y = similar(x, T <: Real ? Complex{T} : T)
        LinearAlgebra.mul!(y, p, x)
        return y[1:end÷2 + 1]
    else
        x_tmp = similar(x, p.flen)
        x_tmp[1:end÷2 + 1] .= x
        x_tmp[end÷2 + 2:end] .= iseven(p.flen) ? conj.(x[end-1:-1:2]) : conj.(x[end:-1:2])
        y = similar(x_tmp)
        LinearAlgebra.mul!(y, p, x_tmp)
        return y
    end
end

function Base.:*(p::FFTAPlan_re{T,N}, x::AbstractArray{T,2}) where {T<:Union{Real, Complex}, N}
    if p.dir == FFT_FORWARD
        y = similar(x, T <: Real ? Complex{T} : T)
        LinearAlgebra.mul!(y, p, x)
        return y[1:end÷2 + 1,:]
    else
        x_tmp = similar(x, p.flen, size(x)[2])
        x_tmp[1:end÷2 + 1,:] .= x
        x_tmp[end÷2 + 2:end,:] .= iseven(p.flen) ? conj.(x[end-1:-1:2,:]) : conj.(x[end:-1:2,:])
        y = similar(x_tmp)
        LinearAlgebra.mul!(y, p, x_tmp)
        return y
    end
end