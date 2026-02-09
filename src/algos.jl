@inline function direction_sign(d::Direction)
    Int(d)
end

@inline _conj(w::Complex, d::Direction) = ifelse(direction_sign(d) === 1, w, conj(w))

function fft!(out::AbstractVector{T}, in::AbstractVector{T}, start_out::Int, start_in::Int, d::Direction, t::FFTEnum, g::CallGraph{T}, idx::Int) where T
    if t === COMPOSITE_FFT
        fft_composite!(out, in, start_out, start_in, d, g, idx)
    else
        root = g[idx]
        if t == DFT
            fft_dft!(out, in, root.sz, start_out, root.s_out, start_in, root.s_in, _conj(root.w, d))
        else
            s_in = root.s_in
            s_out = root.s_out
            if t === POW2RADIX4_FFT
                fft_pow2_radix4!(out, in, root.sz, start_out, s_out, start_in, s_in, _conj(root.w, d))
            elseif t === POW3_FFT
                p_120 = cispi(T(2)/3)
                m_120 = cispi(T(4)/3)
                _p_120, _m_120 = d == FFT_FORWARD ? (p_120, m_120) : (m_120, p_120)
                fft_pow3!(out, in, root.sz, start_out, s_out, start_in, s_in, _conj(root.w, d), _m_120, _p_120)
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
- `out`: Output vector
- `in`: Input vector
- `start_out`: Index of the first element of the output vector
- `start_in`: Index of the first element of the input vector
- `d`: Direction of the transform
- `g`: Call graph for this transform
- `idx`: Index of the current transform in the call graph

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
- `out`: Output vector
- `in`: Input vector
- `N`: Size of the transform
- `start_out`: Index of the first element of the output vector
- `stride_out`: Stride of the output vector
- `start_in`: Index of the first element of the input vector
- `stride_in`: Stride of the input vector
- `w`: The value `cispi(direction_sign(d) * 2 / N)`

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
Radix-4 FFT for powers of 2, in place

# Arguments
- `out`: Output vector
- `in`: Input vector
- `N`: Size of the transform
- `start_out`: Index of the first element of the output vector
- `stride_out`: Stride of the output vector
- `start_in`: Index of the first element of the input vector
- `stride_in`: Stride of the input vector
- `w`: The value `cispi(direction_sign(d) * 2 / N)`

"""
function fft_pow2_radix4!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T) where {T, U}
    # If N is 2, compute the size two DFT
    @inbounds if N == 2
        out[start_out]              = in[start_in] + in[start_in + stride_in]
        out[start_out + stride_out] = in[start_in] - in[start_in + stride_in]
        return
    end

    # If N is 4, compute an unrolled radix-2 FFT and return
    minusi = -sign(imag(w)) * im
    @inbounds if N == 4
        xee = in[start_in]
        xoe = in[start_in +   stride_in]
        xeo = in[start_in + 2*stride_in]
        xoo = in[start_in + 3*stride_in]
        xee_p_xeo = xee + xeo
        xee_m_xeo = xee - xeo
        xoe_p_xoo = xoe + xoo
        xoe_m_xoo = -(xoe - xoo) * minusi
        out[start_out]                = xee_p_xeo + xoe_p_xoo
        out[start_out +   stride_out] = xee_m_xeo + xoe_m_xoo
        out[start_out + 2*stride_out] = xee_p_xeo - xoe_p_xoo
        out[start_out + 3*stride_out] = xee_m_xeo - xoe_m_xoo
        return
    end

    # ...othersize split the problem in four and recur
    m = N ÷ 4

    w1 = w
    w2 = w * w1
    w3 = w * w2
    w4 = w * w3

    fft_pow2_radix4!(out, in, m, start_out                 , stride_out, start_in              , stride_in*4, w4)
    fft_pow2_radix4!(out, in, m, start_out +   m*stride_out, stride_out, start_in +   stride_in, stride_in*4, w4)
    fft_pow2_radix4!(out, in, m, start_out + 2*m*stride_out, stride_out, start_in + 2*stride_in, stride_in*4, w4)
    fft_pow2_radix4!(out, in, m, start_out + 3*m*stride_out, stride_out, start_in + 3*stride_in, stride_in*4, w4)

    wkoe = wkeo = wkoo = one(T)

    @inbounds for k in 0:m-1
        kee = start_out +  k          * stride_out
        koe = start_out + (k +     m) * stride_out
        keo = start_out + (k + 2 * m) * stride_out
        koo = start_out + (k + 3 * m) * stride_out
        y_kee, y_koe, y_keo, y_koo = out[kee], out[koe], out[keo], out[koo]
        ỹ_keo = y_keo * wkeo
        ỹ_koe = y_koe * wkoe
        ỹ_koo = y_koo * wkoo
        y_kee_p_y_keo = y_kee + ỹ_keo
        y_kee_m_y_keo = y_kee - ỹ_keo
        ỹ_koe_p_ỹ_koo = ỹ_koe + ỹ_koo
        ỹ_koe_m_ỹ_koo = -(ỹ_koe - ỹ_koo) * minusi
        out[kee] = y_kee_p_y_keo + ỹ_koe_p_ỹ_koo
        out[koe] = y_kee_m_y_keo + ỹ_koe_m_ỹ_koo
        out[keo] = y_kee_p_y_keo - ỹ_koe_p_ỹ_koo
        out[koo] = y_kee_m_y_keo - ỹ_koe_m_ỹ_koo
        wkoe *= w1
        wkeo *= w2
        wkoo *= w3
    end
end


"""
$(TYPEDSIGNATURES)
Power of 3 FFT, in place

# Arguments
- `out`: Output vector
- `in`: Input vector
- `N`: Size of the transform
- `start_out`: Index of the first element of the output vector
- `stride_out`: Stride of the output vector
- `start_in`: Index of the first element of the input vector
- `stride_in`: Stride of the input vector
- `w`: The value `cispi(direction_sign(d) * 2 / N)`
- `plus120`: Depending on direction, perform either ±120° rotation
- `minus120`: Depending on direction, perform either ∓120° rotation

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
    w2 = w * w1
    wk1 = wk2 = one(T)
    for k in 0:Nprime-1
        @muladd k0 = start_out + stride_out * k
        @muladd k1 = start_out + stride_out * (k + Nprime)
        @muladd k2 = start_out + stride_out * (k + 2 * Nprime)
        y_k0, y_k1, y_k2 = out[k0], out[k1], out[k2]
        @muladd out[k0] = y_k0 + y_k1 * wk1            + y_k2 * wk2
        @muladd out[k1] = y_k0 + y_k1 * wk1 * plus120  + y_k2 * wk2 * minus120
        @muladd out[k2] = y_k0 + y_k1 * wk1 * minus120 + y_k2 * wk2 * plus120
        wk1 *= w1
        wk2 *= w2
    end
end
