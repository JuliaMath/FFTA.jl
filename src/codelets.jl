# Straight-line "codelets" for the small power-of-two base cases of
# `fft_pow2_radix4!`.
#
# A radix-2 decimation-in-time transform of size `N` is unrolled completely at
# compile time by `@generated` functions keyed on `Val{N}` and the direction:
# every intermediate value lives in its own SSA variable, twiddle factors are
# folded to constants and trivial multiplications (by 1, -1, ±i) are removed.
# This is what FFTW's `genfft` does offline; here Julia's compiler does it on
# first use for each `(N, T, direction)`, at a cost of roughly 0.1-0.6 s each
# (see the `PrecompileTools` workload at the end of the module). The codelets
# are only used for `Complex{Float32}` and `Complex{Float64}`; other element
# types keep the generic recursion, which also serves as the reference
# implementation.

const CODELET_MAX = 64
const CodeletEltype = Union{ComplexF32,ComplexF64}

# Emit the statements computing the length-`length(xs)` DFT of the symbols
# `xs` (radix-2 DIT), returning the symbols holding the outputs.
function _gen_dit!(stmts::Vector{Any}, xs::Vector{Symbol}, ::Type{T}, dir::Int, counter::Ref{Int}) where {T}
    N = length(xs)
    N == 1 && return xs
    newsym() = Symbol(:t, counter[] += 1)
    E = _gen_dit!(stmts, xs[1:2:end], T, dir, counter)
    O = _gen_dit!(stmts, xs[2:2:end], T, dir, counter)
    outs = Vector{Symbol}(undef, N)
    for k in 0:N÷2-1
        o = O[k + 1]
        t = newsym()
        if k == 0
            push!(stmts, :($t = $o))
        elseif 4k == N   # multiply by ∓i
            push!(stmts, dir < 0 ? :($t = Complex{$T}(imag($o), -real($o))) :
                                   :($t = Complex{$T}(-imag($o), real($o))))
        else
            w = cispi(dir * 2 * k / N)
            wr, wi = T(real(w)), T(imag(w))
            # explicit fma (not muladd) so that rounding does not depend on
            # whether LLVM contracts for a particular array type
            push!(stmts, :($t = Complex{$T}(fma($wr, real($o), -$wi * imag($o)),
                                           fma($wr, imag($o), $wi * real($o)))))
        end
        a, b = newsym(), newsym()
        push!(stmts, :($a = $(E[k + 1]) + $t), :($b = $(E[k + 1]) - $t))
        outs[k + 1] = a
        outs[k + 1 + N÷2] = b
    end
    return outs
end

"""
$(TYPEDSIGNATURES)
Fully unrolled length-`N` FFT of the strided input `in[start_in + k*stride_in]`
into `out[start_out + k*stride_out]`, `k = 0..N-1`, in direction `D`
(`-1` forward, `+1` backward).
"""
@generated function fft_pow2_codelet!(
    out::AbstractVector{Complex{T}}, in::AbstractVector{Complex{T}},
    ::Val{N},
    start_out::Int, stride_out::Int,
    start_in::Int, stride_in::Int,
    ::Val{D}
) where {T,N,D}
    counter = Ref(0)
    stmts = Any[]
    xs = [Symbol(:x, i) for i in 1:N]
    for i in 1:N
        push!(stmts, :($(xs[i]) = in[start_in + $(i - 1) * stride_in]))
    end
    outs = _gen_dit!(stmts, xs, T, D, counter)
    for i in 1:N
        push!(stmts, :(out[start_out + $(i - 1) * stride_out] = $(outs[i])))
    end
    body = Expr(:block, stmts...)
    return quote
        @inbounds $body
        return nothing
    end
end

# Dispatch a small power of two to its codelet. Returns `false` when there is
# no codelet for this `N` (the caller then continues with the recursion).
@inline function _pow2_codelet!(
    out::AbstractVector{T}, in::AbstractVector{T}, N::Int,
    start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, d::Direction
) where {T<:CodeletEltype}
    if d === FFT_FORWARD
        N == 64 && (fft_pow2_codelet!(out, in, Val(64), start_out, stride_out, start_in, stride_in, Val(-1)); return true)
        N == 32 && (fft_pow2_codelet!(out, in, Val(32), start_out, stride_out, start_in, stride_in, Val(-1)); return true)
        N == 16 && (fft_pow2_codelet!(out, in, Val(16), start_out, stride_out, start_in, stride_in, Val(-1)); return true)
        N == 8  && (fft_pow2_codelet!(out, in, Val(8),  start_out, stride_out, start_in, stride_in, Val(-1)); return true)
    else
        N == 64 && (fft_pow2_codelet!(out, in, Val(64), start_out, stride_out, start_in, stride_in, Val(1)); return true)
        N == 32 && (fft_pow2_codelet!(out, in, Val(32), start_out, stride_out, start_in, stride_in, Val(1)); return true)
        N == 16 && (fft_pow2_codelet!(out, in, Val(16), start_out, stride_out, start_in, stride_in, Val(1)); return true)
        N == 8  && (fft_pow2_codelet!(out, in, Val(8),  start_out, stride_out, start_in, stride_in, Val(1)); return true)
    end
    return false
end
_pow2_codelet!(out, in, N, start_out, stride_out, start_in, stride_in, d) = false
