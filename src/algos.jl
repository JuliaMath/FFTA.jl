function fft!(
    out::AbstractVector{T}, in::AbstractVector{<:Number},
    start_out::Int, start_in::Int,
    d::Direction,
    t::FFTEnum,
    g::CallGraph{T},
    idx::Int
) where T
    if d !== g.dir
        throw(ArgumentError("call graph was planned for direction $(g.dir), not $d"))
    end
    if t === COMPOSITE_FFT
        fft_composite!(out, in, start_out, start_in, d, g, idx)
    else
        root = g[idx]
        s_in = root.s_in
        s_out = root.s_out
        N = root.sz
        tw = g.twiddles[idx]
        if t === DFT
            fft_dft!(out, in, N, start_out, s_out, start_in, s_in, tw)
        elseif t === POW2RADIX4_FFT
            fft_pow2_radix4!(out, in, N, start_out, s_out, start_in, s_in, d, tw, 0)
        elseif t === POW3_FFT
            _m_120 = cispi(T(2) / 3)
            m_120 = d === FFT_FORWARD ? _m_120 : conj(_m_120)
            fft_pow3!(out, in, N, start_out, s_out, start_in, s_in, m_120, d, tw, 0)
        elseif t === BLUESTEIN
            fft_bluestein!(out, in, d, N, start_out, s_out, start_in, s_in, g.bluestein[g.blue_index[idx]])
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
function fft_composite!(
    out::AbstractVector{T}, in::AbstractVector{U},
    start_out::Int, start_in::Int,
    d::Direction,
    g::CallGraph{T},
    idx::Int
) where {T,U}
    root = g[idx]
    left_idx = idx + root.left
    right_idx = idx + root.right
    left = g[left_idx]
    right = g[right_idx]
    N1 = left.sz
    N2 = right.sz
    s_in = root.s_in
    s_out = root.s_out

    Rt = right.type
    Lt = left.type

    tmp = g.workspace[idx]
    tw = g.twiddles[idx]   # see `composite_twiddles`

    for j1 in 0:N1-1
        R_start_in  = start_in + j1 * s_in
        R_start_out = 1 + N2 * j1

        fft!(tmp, in, R_start_out, R_start_in, d, Rt, g, right_idx)

        if j1 > 0
            base = (j1 - 1) * (N2 - 1)
            @inbounds for k2 in 1:N2-1
                tmp[R_start_out + k2] *= tw[base + k2]
            end
        end
    end

    for k2 in 0:N2-1
        L_start_out = start_out + k2 * s_out
        L_start_in  = 1 + k2
        fft!(out, tmp, L_start_out, L_start_in, d, Lt, g, left_idx)
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
- `W`: Twiddle table, see `dft_twiddles` (or `d`, the direction, in which
  case the table is computed on the fly)

"""
function fft_dft!(
    out::AbstractVector{T}, in::AbstractVector{T},
    N::Int,
    start_out::Int, stride_out::Int,
    start_in::Int, stride_in::Int,
    W::AbstractVector{T}
) where {T<:Complex}
    @inbounds begin
        tmp = in[start_in]
        for j in 1:N-1
            tmp += in[start_in + j*stride_in]
        end
        out[start_out] = tmp

        for j in 1:N-1
            tmp = in[start_in]
            idx = 0   # j * k mod N
            for k in 1:N-1
                idx += j
                idx >= N && (idx -= N)
                tmp += W[idx + 1] * in[start_in + k*stride_in]
            end
            out[start_out + j*stride_out] = tmp
        end
    end
end

function fft_dft!(
    out::AbstractVector{Complex{T}}, in::AbstractVector{T},
    N::Int,
    start_out::Int, stride_out::Int,
    start_in::Int, stride_in::Int,
    W::AbstractVector{Complex{T}}
) where {T<:Real}
    halfN = N÷2

    @inbounds begin
        tmp = Complex{T}(in[start_in])
        for j in 1:N-1
            tmp += in[start_in + j*stride_in]
        end
        out[start_out] = tmp

        for j in 1:halfN
            tmp = Complex{T}(in[start_in])
            idx = 0
            for k in 1:N-1
                idx += j
                idx >= N && (idx -= N)
                tmp += W[idx + 1] * in[start_in + k*stride_in]
            end
            out[start_out + j*stride_out] = tmp
            out[start_out + (N-j)*stride_out] = conj(tmp)
        end
    end
end

fft_dft!(out::AbstractVector{T}, in::AbstractVector, N::Int, start_out::Int, stride_out::Int,
         start_in::Int, stride_in::Int, d::Direction) where {T<:Complex} =
    fft_dft!(out, in, N, start_out, stride_out, start_in, stride_in, dft_twiddles(T, N, d))


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
- `d`: Direction of the transform
- `tw`: Twiddle table, see `pow2_twiddles` (omit it to compute the table on the fly)
- `toff`: Offset of the current recursion level in `tw`

"""
function fft_pow2_radix4!(
    out::AbstractVector{T}, in::AbstractVector{U},
    N::Int,
    start_out::Int, stride_out::Int,
    start_in::Int, stride_in::Int,
    d::Direction,
    tw::AbstractVector{T}, toff::Int
) where {T<:Complex, U}
    # If N is 2, compute the size two DFT
    @inbounds if N == 2
        out[start_out]              = in[start_in] + in[start_in + stride_in]
        out[start_out + stride_out] = in[start_in] - in[start_in + stride_in]
        return
    end

    dir = direction_sign(d)

    # If N is 4, compute an unrolled radix-2 FFT and return
    minusi = -dir * im
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
    toff_next = toff + 3m   # the next level's table follows this level's

    fft_pow2_radix4!(out, in, m, start_out                 , stride_out, start_in              , stride_in*4, d, tw, toff_next)
    fft_pow2_radix4!(out, in, m, start_out +   m*stride_out, stride_out, start_in +   stride_in, stride_in*4, d, tw, toff_next)
    fft_pow2_radix4!(out, in, m, start_out + 2*m*stride_out, stride_out, start_in + 2*stride_in, stride_in*4, d, tw, toff_next)
    fft_pow2_radix4!(out, in, m, start_out + 3*m*stride_out, stride_out, start_in + 3*stride_in, stride_in*4, d, tw, toff_next)

    @inbounds for k in 0:m-1
        wkoe = tw[toff + 3k + 1]
        wkeo = tw[toff + 3k + 2]
        wkoo = tw[toff + 3k + 3]
        kee = start_out +  k          * stride_out
        koe = start_out + (k +     m) * stride_out
        keo = start_out + (k + 2 * m) * stride_out
        koo = start_out + (k + 3 * m) * stride_out
        y_kee, y_koe, y_keo, y_koo = out[kee], out[koe], out[keo], out[koo]
        t_koe = y_koe * wkoe
        t_keo = y_keo * wkeo
        t_koo = y_koo * wkoo
        y_kee_p_y_keo = y_kee + t_keo
        y_kee_m_y_keo = y_kee - t_keo
        t_koe_p_t_koo = t_koe + t_koo
        t_koe_m_t_koo = -(t_koe - t_koo) * minusi
        out[kee] = y_kee_p_y_keo + t_koe_p_t_koo
        out[koe] = y_kee_m_y_keo + t_koe_m_t_koo
        out[keo] = y_kee_p_y_keo - t_koe_p_t_koo
        out[koo] = y_kee_m_y_keo - t_koe_m_t_koo
    end
end

fft_pow2_radix4!(out::AbstractVector{T}, in::AbstractVector, N::Int, start_out::Int, stride_out::Int,
                 start_in::Int, stride_in::Int, d::Direction) where {T<:Complex} =
    fft_pow2_radix4!(out, in, N, start_out, stride_out, start_in, stride_in, d, pow2_twiddles(T, N, d), 0)


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
- `minus120`: Depending on direction, perform either ∓120° rotation
- `d`: Direction of the transform
- `tw`: Twiddle table, see `pow3_twiddles` (omit it to compute the table on the fly)
- `toff`: Offset of the current recursion level in `tw`

"""
function fft_pow3!(
    out::AbstractVector{T}, in::AbstractVector{U},
    N::Int,
    start_out::Int, stride_out::Int,
    start_in::Int, stride_in::Int,
    minus120::T,
    d::Direction,
    tw::AbstractVector{T}, toff::Int
) where {T, U}
    plus120 = conj(minus120)
    if N == 3
        @muladd out[start_out + 0]            = in[start_in] + in[start_in + stride_in]          + in[start_in + 2*stride_in]
        @muladd out[start_out +   stride_out] = in[start_in] + in[start_in + stride_in]*plus120  + in[start_in + 2*stride_in]*minus120
        @muladd out[start_out + 2*stride_out] = in[start_in] + in[start_in + stride_in]*minus120 + in[start_in + 2*stride_in]*plus120
        return
    end

    # Size of subproblem
    Nprime = N ÷ 3
    toff_next = toff + 2 * Nprime

    # Dividing into subproblems
    fft_pow3!(out, in, Nprime, start_out,                       stride_out, start_in,               stride_in*3, minus120, d, tw, toff_next)
    fft_pow3!(out, in, Nprime, start_out +   Nprime*stride_out, stride_out, start_in +   stride_in, stride_in*3, minus120, d, tw, toff_next)
    fft_pow3!(out, in, Nprime, start_out + 2*Nprime*stride_out, stride_out, start_in + 2*stride_in, stride_in*3, minus120, d, tw, toff_next)

    @inbounds for k in 0:Nprime-1
        wk1 = tw[toff + 2k + 1]
        wk2 = tw[toff + 2k + 2]
        k0 = start_out + stride_out * k
        k1 = start_out + stride_out * (k + Nprime)
        k2 = start_out + stride_out * (k + 2 * Nprime)
        y_k0, y_k1, y_k2 = out[k0], out[k1], out[k2]
        @muladd out[k0] = y_k0 + y_k1 * wk1            + y_k2 * wk2
        @muladd out[k1] = y_k0 + y_k1 * wk1 * plus120  + y_k2 * wk2 * minus120
        @muladd out[k2] = y_k0 + y_k1 * wk1 * minus120 + y_k2 * wk2 * plus120
    end
end

fft_pow3!(out::AbstractVector{T}, in::AbstractVector, N::Int, start_out::Int, stride_out::Int,
          start_in::Int, stride_in::Int, minus120::T, d::Direction) where {T} =
    fft_pow3!(out, in, N, start_out, stride_out, start_in, stride_in, minus120, d, pow3_twiddles(T, N, d), 0)


"""
$(TYPEDSIGNATURES)
Bluestein's algorithm, still O(N * log(N)) for large primes,
but with a big constant factor.
Zero-pads two sequences derived from the DFT formula to a
3-smooth length ≥ `2N-1` (see `bluestein_pad_length`) and computes their
convolution with FFTs of that length. The chirp, its transform and the work arrays are
precomputed in `scratch` (see `BluesteinScratch`).

# Arguments
- `out`: Output vector
- `in`: Input vector
- `d`: Direction of the transform
- `N`: Size of the transform
- `start_out`: Index of the first element of the output vector
- `stride_out`: Stride of the output vector
- `start_in`: Index of the first element of the input vector
- `stride_in`: Stride of the input vector
- `scratch`: precomputed data and scratch space (omit it to compute it on the fly)

"""
function fft_bluestein!(
    out::AbstractVector{T}, in::AbstractVector{<:Number},
    d::Direction,
    N::Int,
    start_out::Int, stride_out::Int,
    start_in::Int,  stride_in::Int,
    scratch::BluesteinScratch{T}=BluesteinScratch{T}(N, d)
) where T<:Complex
    (; pad_len, chirp, chirp_fft, a, tmp, graph) = scratch
    gt = graph[1].type

    # a_n = x_n · conj(b_n), zero padded
    @inbounds for i in 1:N
        a[i] = in[start_in + (i-1)*stride_in] * conj(chirp[i])
    end
    @inbounds for i in N+1:pad_len
        a[i] = zero(T)
    end

    # Circular convolution of `a` with the periodised chirp via the forward
    # transform only: conv = conj(fft(conj(fft(a) .* fft(b)))) / pad_len, with
    # the 1/pad_len already folded into `chirp_fft`.
    fft!(tmp, a, 1, 1, FFT_FORWARD, gt, graph, 1)
    @inbounds for i in 1:pad_len
        tmp[i] = conj(tmp[i] * chirp_fft[i])
    end
    fft!(a, tmp, 1, 1, FFT_FORWARD, gt, graph, 1)

    # X_k = conj(b_k) · conv_k
    @inbounds for i in 1:N
        out[start_out + (i-1)*stride_out] = conj(chirp[i]) * conj(a[i])
    end
    return nothing
end
