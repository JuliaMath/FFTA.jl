@inline function direction_sign(d::Direction)
    Int(d)
end

@inline _conj(w::Complex, d::Direction) = ifelse(direction_sign(d) === 1, w, conj(w))

function fft!(
    out::AbstractVector{T}, in::AbstractVector{T},
    start_out::Int, start_in::Int,
    d::Direction,
    t::FFTEnum,
    g::CallGraph{T},
    idx::Int
    ) where T
    if t === COMPOSITE_FFT
        fft_composite!(out, in, start_out, start_in, d, g, idx)
    else
        root = g[idx]
        s_in = root.s_in
        s_out = root.s_out
        N = root.sz
        w = _conj(root.w, d)
        if t === DFT
            fft_dft!(out, in, N, start_out, s_out, start_in, s_in, w)
        elseif t === POW2RADIX4_FFT
            fft_pow2_radix4!(out, in, N, start_out, s_out, start_in, s_in, w)
        elseif t === POW3_FFT
            _m_120 = cispi(T(2) / 3)
            m_120 = d === FFT_FORWARD ? _m_120 : conj(_m_120)
            fft_pow3!(out, in, N, start_out, s_out, start_in, s_in, w, m_120)
        elseif t === BLUESTEIN
            fft_bluestein!(out, in, d, N, start_out, s_out, start_in, s_in)
        else
            throw(ArgumentError("kernel not implemented"))
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

    Rt = right.type
    Lt = left.type

    w1 = _conj(root.w, d)
    Rtype = real(T)
    # The composite twiddle at position (j1, k2) is `cispi(dir · 2 j1 k2 / N)`.
    # Singleton's recurrence advances `wk2 = cispi(dir · 2 j1 k2 / N)` in k2
    # for fixed j1; (α, β) depend on j1 so we reset them at each outer step.
    dir = twiddle_direction(w1)
    tmp = g.workspace[idx]

    if Rt === BLUESTEIN
        R_bluestein_scratchspace = prealloc_blue(N2, d, T)
    end
    for j1 in 0:N1-1
        R_start_in  = start_in + j1 * s_in
        R_start_out = 1 + N2 * j1

        if @isdefined R_bluestein_scratchspace
            R_s_in  = right.s_in
            R_s_out = right.s_out
            fft_bluestein!(tmp, in, d, N2, R_start_out, R_s_out, R_start_in, R_s_in, R_bluestein_scratchspace)
        else
            fft!(tmp, in, R_start_out, R_start_in, d, Rt, g, right_idx)
        end

        if j1 > 0
            αi, βi = singleton_params(dir * Rtype(2 * j1) / Rtype(N))
            ci, si = one(Rtype), zero(Rtype)
            @inbounds for k2 in 1:N2-1
                ci, si = singleton_step(ci, si, αi, βi)
                tmp[R_start_out + k2] *= Complex(ci, si)
            end
        end
    end

    if Lt === BLUESTEIN
        L_bluestein_scratchspace = prealloc_blue(N1, d, T)
    end
    for k2 in 0:N2-1
        L_start_out = start_out + k2 * s_out
        L_start_in  = 1 + k2
        if @isdefined L_bluestein_scratchspace
            L_s_in  = left.s_in
            L_s_out = left.s_out
            fft_bluestein!(out, tmp, d, N1, L_start_out, L_s_out, L_start_in, L_s_in, L_bluestein_scratchspace)
        else
            fft!(out, tmp, L_start_out, L_start_in, d, Lt, g, left_idx)
        end
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

    Rtype = real(T)
    dir = twiddle_direction(w)
    @inbounds for d in 1:N-1
        t = in[start_in]
        αk, βk = singleton_params(dir * Rtype(2 * d) / Rtype(N))
        ck, sk = one(Rtype), zero(Rtype)
        @inbounds for k in 1:N-1
            ck, sk = singleton_step(ck, sk, αk, βk)
            t += Complex(ck, sk) * in[start_in + k*stride_in]
        end
        out[start_out + d*stride_out] = t
    end
end

function fft_dft!(out::AbstractVector{Complex{T}}, in::AbstractVector{T}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::Complex{T}) where {T<:Real}
    halfN = N÷2

    tmp = Complex{T}(in[start_in])
    @inbounds for j in 1:N-1
        tmp += in[start_in + j*stride_in]
    end
    out[start_out] = tmp

    dir = twiddle_direction(w)
    @inbounds for d in 1:halfN
        t = Complex{T}(in[start_in])
        αk, βk = singleton_params(dir * T(2 * d) / T(N))
        ck, sk = one(T), zero(T)
        @inbounds for k in 1:N-1
            ck, sk = singleton_step(ck, sk, αk, βk)
            t += Complex{T}(ck, sk) * in[start_in + k*stride_in]
        end
        out[start_out + d*stride_out] = t
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

    Rtype = real(T)
    dir = twiddle_direction(w)
    # Recursive sub-problem step `cispi(dir · 2 / m) = w^4`; use `cispi`
    # directly so the sub-tree gets a < 1 ULP starting phase.
    w_sub = cispi(dir * Rtype(2) / Rtype(m))

    fft_pow2_radix4!(out, in, m, start_out                 , stride_out, start_in              , stride_in*4, w_sub)
    fft_pow2_radix4!(out, in, m, start_out +   m*stride_out, stride_out, start_in +   stride_in, stride_in*4, w_sub)
    fft_pow2_radix4!(out, in, m, start_out + 2*m*stride_out, stride_out, start_in + 2*stride_in, stride_in*4, w_sub)
    fft_pow2_radix4!(out, in, m, start_out + 3*m*stride_out, stride_out, start_in + 3*stride_in, stride_in*4, w_sub)

    # Singleton recurrence for the three running twiddles `w^k`, `w^2k`, `w^3k`.
    α1, β1 = singleton_params(dir * Rtype(2) / Rtype(N))
    α2, β2 = singleton_params(dir * Rtype(4) / Rtype(N))
    α3, β3 = singleton_params(dir * Rtype(6) / Rtype(N))
    c1, s1 = one(Rtype), zero(Rtype)
    c2, s2 = one(Rtype), zero(Rtype)
    c3, s3 = one(Rtype), zero(Rtype)

    @inbounds for k in 0:m-1
        kee = start_out +  k          * stride_out
        koe = start_out + (k +     m) * stride_out
        keo = start_out + (k + 2 * m) * stride_out
        koo = start_out + (k + 3 * m) * stride_out
        y_kee, y_koe, y_keo, y_koo = out[kee], out[koe], out[keo], out[koo]
        t_keo = y_keo * Complex(c2, s2)
        t_koe = y_koe * Complex(c1, s1)
        t_koo = y_koo * Complex(c3, s3)
        y_kee_p_y_keo = y_kee + t_keo
        y_kee_m_y_keo = y_kee - t_keo
        t_koe_p_t_koo = t_koe + t_koo
        t_koe_m_t_koo = -(t_koe - t_koo) * minusi
        out[kee] = y_kee_p_y_keo + t_koe_p_t_koo
        out[koe] = y_kee_m_y_keo + t_koe_m_t_koo
        out[keo] = y_kee_p_y_keo - t_koe_p_t_koo
        out[koo] = y_kee_m_y_keo - t_koe_m_t_koo
        c1, s1 = singleton_step(c1, s1, α1, β1)
        c2, s2 = singleton_step(c2, s2, α2, β2)
        c3, s3 = singleton_step(c3, s3, α3, β3)
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
function fft_pow3!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T, minus120::T) where {T, U}
    plus120 = conj(minus120)
    if N == 3
        @muladd out[start_out + 0]            = in[start_in] + in[start_in + stride_in]          + in[start_in + 2*stride_in]
        @muladd out[start_out +   stride_out] = in[start_in] + in[start_in + stride_in]*plus120  + in[start_in + 2*stride_in]*minus120
        @muladd out[start_out + 2*stride_out] = in[start_in] + in[start_in + stride_in]*minus120 + in[start_in + 2*stride_in]*plus120
        return
    end

    # Size of subproblem
    Nprime = N ÷ 3

    Rtype = real(T)
    dir = twiddle_direction(w)
    # Recursive sub-problem step cispi(dir · 2 / Nprime) = w^3.
    w_sub = cispi(dir * Rtype(2) / Rtype(Nprime))

    # Dividing into subproblems
    fft_pow3!(out, in, Nprime, start_out,                       stride_out, start_in,               stride_in*3, w_sub, minus120)
    fft_pow3!(out, in, Nprime, start_out +   Nprime*stride_out, stride_out, start_in +   stride_in, stride_in*3, w_sub, minus120)
    fft_pow3!(out, in, Nprime, start_out + 2*Nprime*stride_out, stride_out, start_in + 2*stride_in, stride_in*3, w_sub, minus120)

    α1, β1 = singleton_params(dir * Rtype(2) / Rtype(N))
    α2, β2 = singleton_params(dir * Rtype(4) / Rtype(N))
    c1, s1 = one(Rtype), zero(Rtype)
    c2, s2 = one(Rtype), zero(Rtype)
    for k in 0:Nprime-1
        k0 = start_out + stride_out * k
        k1 = start_out + stride_out * (k + Nprime)
        k2 = start_out + stride_out * (k + 2 * Nprime)
        y_k0, y_k1, y_k2 = out[k0], out[k1], out[k2]
        wk1 = Complex(c1, s1)
        wk2 = Complex(c2, s2)
        @muladd out[k0] = y_k0 + y_k1*wk1 + y_k2*wk2
        @muladd out[k1] = y_k0 + y_k1*wk1*plus120 + y_k2*wk2*minus120
        @muladd out[k2] = y_k0 + y_k1*wk1*minus120 + y_k2*wk2*plus120
        c1, s1 = singleton_step(c1, s1, α1, β1)
        c2, s2 = singleton_step(c2, s2, α2, β2)
    end
end


function prealloc_blue(N::Int, d::Direction, ::Type{T}) where T<:Number
    pad_len = nextpow(2, 2N - 1)

    b_series = Vector{T}(undef, pad_len)
    a_series = Vector{T}(undef, pad_len)
    tmp      = Vector{T}(undef, pad_len)

    b_series[N+1:end] .= zero(T)

    sgn = -direction_sign(d)
    p = 0   # n^2
    for i in 1:N
        b_series[i] = cispi(sgn * p / N)
        p += (2i - 1)   # prevents overflow unless N is absolutely massive
        p > N && (p -= 2N)
    end

    # enforce periodic boundaries for b_n
    for j in 0:N-1
        b_series[pad_len-j] = b_series[2+j]
    end

    return (tmp, a_series, b_series, pad_len)
end

"""
$(TYPEDSIGNATURES)
Bluestein's algorithm, still O(N * log(N)) for large primes,
but with a big constant factor.
Zero-pads two sequences derived from the DFT formula to a
power of 2 length greater than `2N-1` and computes their convolution
with a power 2 FFT.

# Arguments
- `out`: Output vector
- `in`: Input vector
- `d`: Direction of the transform
- `N`: Size of the transform
- `start_out`: Index of the first element of the output vector
- `stride_out`: Stride of the output vector
- `start_in`: Index of the first element of the input vector
- `stride_in`: Stride of the input vector
- `w`: The value `cispi(direction_sign(d) * 2 / N)`

"""
function fft_bluestein!(
    out::AbstractVector{T}, in::AbstractVector{T},
    d::Direction,
    N::Int,
    start_out::Int, stride_out::Int,
    start_in::Int,  stride_in::Int,
    scratch::Tuple{Vector{T},Vector{T},Vector{T},Int}=prealloc_blue(N, d, T)
) where T<:Number

    (tmp, a_series, b_series, pad_len) = scratch

    a_series[N+1:end] .= zero(T)
    tmp[N+1:end]      .= zero(T)

    for i in 1:N
        a_series[i] = in[start_in+(i-1)*stride_in] * conj(b_series[i])
    end

    w_pad = cispi(T(2) / pad_len)
    # leave b_n vector alone for last step
    fft_pow2_radix4!(tmp,      a_series, pad_len, 1, 1, 1, 1, w_pad)    # Fa
    fft_pow2_radix4!(a_series, b_series, pad_len, 1, 1, 1, 1, w_pad)    # Fb

    tmp .*= a_series
    # convolution theorem ifft
    fft_pow2_radix4!(a_series, tmp, pad_len, 1, 1, 1, 1, conj(w_pad))
    conv_a_b = a_series

    Xk = tmp
    for i in 1:N
        Xk[i] = conj(b_series[i]) * conv_a_b[i] / pad_len
    end

    out_inds = range(start_out; step=stride_out, length=N)
    copyto!(out, CartesianIndices((out_inds,)), Xk, CartesianIndices((N,)))
    return nothing
end
