fft!(::AbstractVector{T}, ::AbstractVector{T}, ::Int, ::Int, ::Direction, ::AbstractFFTType, ::CallGraph{T}, ::Int) where {T} = nothing

@inline function direction_sign(d::Direction)
    Int(d)
end

function (g::CallGraph{T})(out::AbstractVector{T}, in::AbstractVector{U}, start_out::Int, start_in::Int, v::Direction, t::AbstractFFTType, idx::Int) where {T,U}
    fft!(out, in, start_out, start_in, v, t, g, idx)
end

"""
$(TYPEDSIGNATURES)
Cooley-Tukey composite FFT, with a pre-computed call graph

# Arguments
`out`: Output vector
`in`: Input vector
`start_out`: Index of the first element of the output vector
`start_in`: Index of the first element of the input vector
`d`: Direction of the transform
`g`: Call graph for this transform
`idx`: Index of the current transform in the call graph

"""
function fft!(out::AbstractVector{T}, in::AbstractVector{U}, start_out::Int, start_in::Int, d::Direction, ::CompositeFFT, g::CallGraph{T}, idx::Int) where {T,U}
    root = g[idx]
    left_idx = idx + root.left
    right_idx = idx + root.right
    left = g[left_idx]
    right = g[right_idx]
    N  = root.sz
    N1 = left.sz
    N2 = right.sz
    s_in = root.s_in
    s_out = root.s_out

    w1 = convert(T, cispi(direction_sign(d)*2/N))
    wj1 = one(T)
    tmp = g.workspace[idx]
    @inbounds for j1 in 0:N1-1
        wk2 = wj1
        g(tmp, in, N2*j1+1, start_in + j1*s_in, d, right.type, right_idx)
        j1 > 0 && @inbounds for k2 in 1:N2-1
            tmp[N2*j1 + k2 + 1] *= wk2
            wk2 *= wj1
        end
        wj1 *= w1
    end

    @inbounds for k2 in 0:N2-1
        g(out, tmp, start_out + k2*s_out, k2+1, d, left.type, left_idx)
    end
end

"""
$(TYPEDSIGNATURES)
Discrete Fourier Transform, O(N^2) algorithm, in place.

# Arguments
`out`: Output vector
`in`: Input vector
`N`: Size of the transform
`start_out`: Index of the first element of the output vector
`stride_out`: Stride of the output vector
`start_in`: Index of the first element of the input vector
`stride_in`: Stride of the input vector
`d`: Direction of the transform

"""
function fft_dft!(out::AbstractVector{T}, in::AbstractVector{T}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, d::Direction) where {T}
    tmp = in[start_in]
    @inbounds for j in 1:N-1
        tmp += in[start_in + j*stride_in]
    end
    out[start_out] = tmp

    wk = wkn = w = convert(T, cispi(direction_sign(d)*2/N))
    @inbounds for d in 1:N-1
        tmp = in[start_in]
        @inbounds for k in 1:N-1
            tmp += wkn*in[start_in + k*stride_in]
            wkn *= wk
        end
        out[start_out + d*stride_out] = tmp
        wk *= w
        wkn = wk
    end
end

function fft_dft!(out::AbstractVector{Complex{T}}, in::AbstractVector{T}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, d::Direction) where {T<:Real}
    halfN = N÷2

    tmp = Complex{T}(in[start_in])
    @inbounds for j in 1:N-1
        tmp += in[start_in + j*stride_in]
    end
    out[start_out] = tmp

    wk = wkn = w = convert(Complex{T}, cispi(direction_sign(d)*2/N))
    @inbounds for d in 1:halfN
        tmp = Complex{T}(in[start_in])
        @inbounds for k in 1:N-1
            tmp += wkn*in[start_in + k*stride_in]
            wkn *= wk
        end
        out[start_out + d*stride_out] = tmp
        wk *= w
        wkn = wk
    end
end

function fft!(out::AbstractVector{T}, in::AbstractVector{U}, start_out::Int, start_in::Int, d::Direction, ::DFT, g::CallGraph{T}, idx::Int) where {T,U}
    root = g[idx]
    fft_dft!(out, in, root.sz, start_out, root.s_out, start_in, root.s_in, d)
end

"""
$(TYPEDSIGNATURES)
Power of 2 FFT, in place

# Arguments
`out`: Output vector
`in`: Input vector
`N`: Size of the transform
`start_out`: Index of the first element of the output vector
`stride_out`: Stride of the output vector
`start_in`: Index of the first element of the input vector
`stride_in`: Stride of the input vector
`d`: Direction of the transform

"""
function fft_pow2!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, d::Direction) where {T, U}
    if N == 2
        out[start_out]              = in[start_in] + in[start_in + stride_in]
        out[start_out + stride_out] = in[start_in] - in[start_in + stride_in]
        return
    end
    m = N ÷ 2

    fft_pow2!(out, in, m, start_out               , stride_out, start_in            , stride_in*2, d)
    fft_pow2!(out, in, m, start_out + m*stride_out, stride_out, start_in + stride_in, stride_in*2, d)

    w1 = convert(T, cispi(direction_sign(d)*2/N))
    wj = one(T)
    @inbounds for j in 0:m-1
        j1_out = start_out + j*stride_out
        j2_out = start_out + (j+m)*stride_out
        out_j    = out[j1_out]
        out[j1_out] = out_j + wj*out[j2_out]
        out[j2_out] = out_j - wj*out[j2_out]
        wj *= w1
    end
end

function fft!(out::AbstractVector{T}, in::AbstractVector{U}, start_out::Int, start_in::Int, d::Direction, ::Pow2FFT, g::CallGraph{T}, idx::Int) where {T,U}
    root = g[idx]
    N = root.sz
    s_in = root.s_in
    s_out = root.s_out
    fft_pow2!(out, in, N, start_out, s_out, start_in, s_in, d)
end

"""
$(TYPEDSIGNATURES)
Power of 4 FFT, in place

# Arguments
`out`: Output vector
`in`: Input vector
`N`: Size of the transform
`start_out`: Index of the first element of the output vector
`stride_out`: Stride of the output vector
`start_in`: Index of the first element of the input vector
`stride_in`: Stride of the input vector
`d`: Direction of the transform

"""
function fft_pow4!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, d::Direction) where {T, U}
    ds = direction_sign(d)
    plusi = ds*1im
    minusi = ds*-1im
    if N == 4
        out[start_out + 0]            = in[start_in] + in[start_in + stride_in]        + in[start_in + 2*stride_in] + in[start_in + 3*stride_in]
        out[start_out +   stride_out] = in[start_in] + in[start_in + stride_in]*plusi  - in[start_in + 2*stride_in] + in[start_in + 3*stride_in]*minusi
        out[start_out + 2*stride_out] = in[start_in] - in[start_in + stride_in]        + in[start_in + 2*stride_in] - in[start_in + 3*stride_in]
        out[start_out + 3*stride_out] = in[start_in] + in[start_in + stride_in]*minusi - in[start_in + 2*stride_in] + in[start_in + 3*stride_in]*plusi
        return
    end
    m = N ÷ 4

    @muladd fft_pow4!(out, in, m, start_out                 , stride_out, start_in              , stride_in*4, d)
    @muladd fft_pow4!(out, in, m, start_out +   m*stride_out, stride_out, start_in +   stride_in, stride_in*4, d)
    @muladd fft_pow4!(out, in, m, start_out + 2*m*stride_out, stride_out, start_in + 2*stride_in, stride_in*4, d)
    @muladd fft_pow4!(out, in, m, start_out + 3*m*stride_out, stride_out, start_in + 3*stride_in, stride_in*4, d)

    w1 = convert(T, cispi(direction_sign(d)*2/N))
    wj = one(T)

    w1 = convert(T, cispi(ds*2/N))
    w2 = convert(T, cispi(ds*4/N))
    w3 = convert(T, cispi(ds*6/N))
    wk1 = wk2 = wk3 = one(T)

    @inbounds for k in 0:m-1
        @muladd k0 = start_out + k*stride_out
        @muladd k1 = start_out + (k+m)*stride_out
        @muladd k2 = start_out + (k+2*m)*stride_out
        @muladd k3 = start_out + (k+3*m)*stride_out
        y_k0, y_k1, y_k2, y_k3 = out[k0], out[k1], out[k2], out[k3]
        @muladd out[k0] = (y_k0 + y_k2*wk2) + (y_k1*wk1 + y_k3*wk3)
        @muladd out[k1] = (y_k0 - y_k2*wk2) + (y_k1*wk1 - y_k3*wk3) * plusi
        @muladd out[k2] = (y_k0 + y_k2*wk2) - (y_k1*wk1 + y_k3*wk3)
        @muladd out[k3] = (y_k0 - y_k2*wk2) + (y_k1*wk1 - y_k3*wk3) * minusi
        wk1 *= w1
        wk2 *= w2
        wk3 *= w3
    end
end

function fft!(out::AbstractVector{T}, in::AbstractVector{U}, start_out::Int, start_in::Int, d::Direction, ::Pow4FFT, g::CallGraph{T}, idx::Int) where {T,U}
    root = g[idx]
    N = root.sz
    s_in = root.s_in
    s_out = root.s_out
    fft_pow4!(out, in, N, start_out, s_out, start_in, s_in, d)
end

"""
$(TYPEDSIGNATURES)
Power of 3 FFT, in place

# Arguments
out: Output vector
in: Input vector
N: Size of the transform
start_out: Index of the first element of the output vector
stride_out: Stride of the output vector
start_in: Index of the first element of the input vector
stride_in: Stride of the input vector
d: Direction of the transform
plus120: Depending on direction, perform either ±120° rotation
minus120: Depending on direction, perform either ∓120° rotation

"""
function fft_pow3!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, d::Direction, plus120::T, minus120::T) where {T, U}
    if N == 3
        @muladd out[start_out + 0]            = in[start_in] + in[start_in + stride_in]          + in[start_in + 2*stride_in]
        @muladd out[start_out +   stride_out] = in[start_in] + in[start_in + stride_in]*plus120  + in[start_in + 2*stride_in]*minus120
        @muladd out[start_out + 2*stride_out] = in[start_in] + in[start_in + stride_in]*minus120 + in[start_in + 2*stride_in]*plus120
        return
    end

    # Size of subproblem
    Nprime = N ÷ 3

    ds = direction_sign(d)

    # Dividing into subproblems
    fft_pow3!(out, in, Nprime, start_out, stride_out, start_in, stride_in*3, d, plus120, minus120)
    fft_pow3!(out, in, Nprime, start_out + Nprime*stride_out, stride_out, start_in + stride_in, stride_in*3, d, plus120, minus120)
    fft_pow3!(out, in, Nprime, start_out + 2*Nprime*stride_out, stride_out, start_in + 2*stride_in, stride_in*3, d, plus120, minus120)

    w1 = convert(T, cispi(ds*2/N))
    w2 = convert(T, cispi(ds*4/N))
    wk1 = wk2 = one(T)
    for k in 0:Nprime-1
        @muladd k0 = start_out + k*stride_out
        @muladd k1 = start_out + (k+Nprime)*stride_out
        @muladd k2 = start_out + (k+2*Nprime)*stride_out
        y_k0, y_k1, y_k2 = out[k0], out[k1], out[k2]
        @muladd out[k0] = y_k0 + y_k1*wk1 + y_k2*wk2
        @muladd out[k1] = y_k0 + y_k1*wk1*plus120 + y_k2*wk2*minus120
        @muladd out[k2] = y_k0 + y_k1*wk1*minus120 + y_k2*wk2*plus120
        wk1 *= w1
        wk2 *= w2
    end
end

function fft!(out::AbstractVector{T}, in::AbstractVector{U}, start_out::Int, start_in::Int, d::Direction, ::Pow3FFT, g::CallGraph{T}, idx::Int) where {T,U}
    root = g[idx]
    N = root.sz
    s_in = root.s_in
    s_out = root.s_out
    p_120 = convert(T, cispi(2/3))
    m_120 = convert(T, cispi(4/3))
    if d == FFT_FORWARD
        fft_pow3!(out, in, N, start_out, s_out, start_in, s_in, d, m_120, p_120)
    else
        fft_pow3!(out, in, N, start_out, s_out, start_in, s_in, d, p_120, m_120)
    end
end