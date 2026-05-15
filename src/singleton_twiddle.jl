# Singleton's stable trigonometric recurrence for twiddle factors.
#
# The naïve `w *= step` recurrence gives O(N) relative error growth,
# which is the issue described in #118. Singleton (1967) showed that
# if one precomputes
#     α = 2 sin²(θ/2) = 1 - cos(θ)
#     β = sin(θ)
# then the recurrence
#     c_{k+1} = c_k - (α c_k + β s_k)
#     s_{k+1} = s_k - (α s_k - β c_k)
# produces (cos(kθ), sin(kθ)) with RMS error growing only as √k · ε
# (and worst-case ~ k · ε / something sub-linear in typical cases)
# rather than k · ε. The point is that `α c_k + β s_k` is a *small*
# correction (α ≈ θ²/2, β ≈ θ) subtracted from an ~O(1) base, so the
# update stays in the high-significand bits.
#
# The cost is ~8 flops per step (2× the naïve complex multiply but an
# order of magnitude faster than a fresh `cispi`), and the extra trig
# (`sincospi(θ/2)`) happens once per kernel call.

# Direction lives in `sign(imag(w))`; when `w` is real (N = 2 or any
# degenerate case where `imag` rounds to zero) both directions collapse
# to the same twiddle set so we pick +1.
@inline function twiddle_direction(w::Complex{T}) where {T<:Real}
    s = imag(w)
    copysign(one(T), s)
end

# Recurrence coefficients for stepping by `cispi(freq) = e^(iπ·freq)`.
# Uses `sincospi(hfreq)` so that `α` and `β` are exact-to-ULP even
# for very small frequencies — writing `1 - cos(θ)` directly suffers
# catastrophic cancellation there.
@inline function singleton_params(hfreq::Real)
    s_h, c_h = sincospi(-hfreq)
    α = 2 * s_h * s_h
    β = 2 * s_h * c_h
    Complex(α, β)
end

# Advance `(c, s) = (cos(kθ), sin(kθ))` to `(cos((k+1)θ), sin((k+1)θ))`.
# Computed as `c - (αc + βs)` rather than `(1-α)c - βs` on purpose:
# the correction is small so subtracting it from `c` preserves the
# high-order bits and the recurrence self-heals.
@inline function singleton_step(w::T, z::T) where {T<:Complex}
    # muladd only reduces instructions, doesn't help precision much
    w - @fastmath(z * w)
end
