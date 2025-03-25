@inline function direction_sign(d::Direction)
    Int(d)
end

@inline _conj(w::Complex, d::Direction) = ifelse(direction_sign(d) === 1, w, conj(w))

function fft!(out::AbstractVector{T}, in::AbstractVector{T}, start_out::Int, start_in::Int, d::Direction, t::FFTEnum, g::CallGraph{T}, idx::Int) where T
    if t === compositeFFT
        fft_composite!(out, in, start_out, start_in, d, g, idx)
    else
        root = g[idx]
        if t == dft
            fft_dft!(out, in, root.sz, start_out, root.s_out, start_in, root.s_in, _conj(root.w, d))
        else
            N = root.sz
            s_in = root.s_in
            s_out = root.s_out
            if t === pow2FFT
                fft_pow2!(out, in, N, start_out, s_out, start_in, s_in, _conj(root.w, d))
            elseif t === pow3FFT
                p_120 = convert(T, cispi(2/3))
                m_120 = convert(T, cispi(4/3))
                _p_120, _m_120 = d == FFT_FORWARD ? (p_120, m_120) : (m_120, p_120)
                fft_pow3!(out, in, N, start_out, s_out, start_in, s_in, _conj(root.w, d), _m_120, _p_120)
            elseif t === pow4FFT
                fft_pow4!(out, in, N, start_out, s_out, start_in, s_in, _conj(root.w, d))
            else
                throw(ArgumentError("kernel not implemented"))
            end
        end
    end
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
function fft_composite!(out::AbstractVector{T}, in::AbstractVector{U}, start_out::Int, start_in::Int, d::Direction, g::CallGraph{T}, idx::Int) where {T,U}
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

    w1 = _conj(root.w, d)
    wj1 = one(T)
    tmp = g.workspace[idx]
    @inbounds for j1 in 0:N1-1
        wk2 = wj1
        fft!(tmp, in, N2*j1+1, start_in + j1*s_in, d, right.type, g, right_idx)
        j1 > 0 && @inbounds for k2 in 1:N2-1
            tmp[N2*j1 + k2 + 1] *= wk2
            wk2 *= wj1
        end
        wj1 *= w1
    end

    @inbounds for k2 in 0:N2-1
        fft!(out, tmp, start_out + k2*s_out, k2+1, d, left.type, g, left_idx)
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
`w`: The value `cispi(direction_sign(d) * 2 / N)`

"""
function fft_dft!(out::AbstractVector{T}, in::AbstractVector{T}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T) where {T}
    tmp = in[start_in]
    @inbounds for j in 1:N-1
        tmp += in[start_in + j*stride_in]
    end
    out[start_out] = tmp

    wk = wkn = w
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

function fft_dft!(out::AbstractVector{Complex{T}}, in::AbstractVector{T}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::Complex{T}) where {T<:Real}
    halfN = N÷2

    tmp = Complex{T}(in[start_in])
    @inbounds for j in 1:N-1
        tmp += in[start_in + j*stride_in]
    end
    out[start_out] = tmp

    wk = wkn = w
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
`w`: The value `cispi(direction_sign(d) * 2 / N)`

"""
function fft_pow2!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T) where {T, U}
    if N == 2
        out[start_out]              = in[start_in] + in[start_in + stride_in]
        out[start_out + stride_out] = in[start_in] - in[start_in + stride_in]
        return
    end
    m = N ÷ 2

    fft_pow2!(out, in, m, start_out               , stride_out, start_in            , stride_in*2, w*w)
    fft_pow2!(out, in, m, start_out + m*stride_out, stride_out, start_in + stride_in, stride_in*2, w*w)

    wj = one(T)
    @inbounds for j in 0:m-1
        j1_out = start_out + j*stride_out
        j2_out = start_out + (j+m)*stride_out
        out_j    = out[j1_out]
        out[j1_out] = out_j + wj*out[j2_out]
        out[j2_out] = out_j - wj*out[j2_out]
        wj *= w
    end
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
`w`: The value `cispi(direction_sign(d) * 2 / N)`

"""
function fft_pow4!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T) where {T, U}
    plusi = sign(imag(w))*im
    minusi = -sign(imag(w))*im
    if N == 4
        out[start_out + 0]            = in[start_in] + in[start_in + stride_in]        + in[start_in + 2*stride_in] + in[start_in + 3*stride_in]
        out[start_out +   stride_out] = in[start_in] + in[start_in + stride_in]*plusi  - in[start_in + 2*stride_in] + in[start_in + 3*stride_in]*minusi
        out[start_out + 2*stride_out] = in[start_in] - in[start_in + stride_in]        + in[start_in + 2*stride_in] - in[start_in + 3*stride_in]
        out[start_out + 3*stride_out] = in[start_in] + in[start_in + stride_in]*minusi - in[start_in + 2*stride_in] + in[start_in + 3*stride_in]*plusi
        return
    end
    m = N ÷ 4

    @muladd fft_pow4!(out, in, m, start_out                 , stride_out, start_in              , stride_in*4, w^4)
    @muladd fft_pow4!(out, in, m, start_out +   m*stride_out, stride_out, start_in +   stride_in, stride_in*4, w^4)
    @muladd fft_pow4!(out, in, m, start_out + 2*m*stride_out, stride_out, start_in + 2*stride_in, stride_in*4, w^4)
    @muladd fft_pow4!(out, in, m, start_out + 3*m*stride_out, stride_out, start_in + 3*stride_in, stride_in*4, w^4)

    w1 = w
    w2 = w*w1
    w3 = w*w2
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
w: The value `cispi(direction_sign(d) * 2 / N)`
plus120: Depending on direction, perform either ±120° rotation
minus120: Depending on direction, perform either ∓120° rotation

"""
function fft_pow3!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T, plus120::T, minus120::T) where {T, U}
    if N == 3
        @muladd out[start_out + 0]            = in[start_in] + in[start_in + stride_in]          + in[start_in + 2*stride_in]
        @muladd out[start_out +   stride_out] = in[start_in] + in[start_in + stride_in]*plus120  + in[start_in + 2*stride_in]*minus120
        @muladd out[start_out + 2*stride_out] = in[start_in] + in[start_in + stride_in]*minus120 + in[start_in + 2*stride_in]*plus120
        return
    end

    # Size of subproblem
    Nprime = N ÷ 3

    # Dividing into subproblems
    fft_pow3!(out, in, Nprime, start_out, stride_out, start_in, stride_in*3, w^3, plus120, minus120)
    fft_pow3!(out, in, Nprime, start_out + Nprime*stride_out, stride_out, start_in + stride_in, stride_in*3, w^3, plus120, minus120)
    fft_pow3!(out, in, Nprime, start_out + 2*Nprime*stride_out, stride_out, start_in + 2*stride_in, stride_in*3, w^3, plus120, minus120)

    w1 = w
    w2 = w*w1
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
