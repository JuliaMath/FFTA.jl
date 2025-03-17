# FFTA: Fastest Fourier Transform in my Apartment
## A library by Danny Sharp

[![CI](https://github.com/JuliaMath/FFTA.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaMath/FFTA.jl/actions/workflows/ci.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliamath.github.io/FFTA.jl/stable/)
[![codecov](https://codecov.io/gh/JuliaMath/FFTA.jl/graph/badge.svg?token=LEUDt0h3CG)](https://codecov.io/gh/JuliaMath/FFTA.jl)

This is a *pure Julia* implementation of FFTs, with the goal that this could supplant other FFTs for applications that require odd Julia objects. Currently this supports `AbstractArray{T,N}` where `N` in `{1,2}` (i.e. `AbstractVector` and `AbstractMatrix`). If you're looking for more performance, checkout `FFTW.jl`. Regardless of `T`, `one(::Type{T})` must be defined. Additionally, if `T<:Real`, then `convert(::Type{T}, Float64)` has to be defined. Otherwise, `convert(::Type{T},ComplexF64)` must be defined.

Some ideas outside the feature requests in Issues:
- Make the code more readable
- Use `StaticArrays` for the workspace in small cases
- Strictly generate code for certain cases
- Create a SIMD type for Complex numbers
- E-Graphs for the call-graph
- Other performance left on the table....

Interface:
- `fft(x::AbstractVector)`-- Forward FFT
- `fft(x::AbstractMatrix})`-- Forward FFT
- `bfft(x::AbstractVector})`-- Backward FFT (unscaled inverse FFT)
- `bfft(x::AbstractMatrix})`-- Backward FFT (unscaled inverse FFT)

NOTE: Currently, my C++ code is actually faster than this, so "Fastest Fourier Transform in my Apartment" is a bit of a misnomer.


## Why use this?
There's a lot of FFT packages out there, no doubt. Many are great. Some, like mine, are "good enough". Many aren't so great. As far as I know, though, very few are as generic as FFTA. Does that matter? Yes. One of the main draws of Julia is the fact that a lot of functions "just work" with types from other packages. FFTA aims to abide by this philosophy. For example, have you ever wanted to generate what an FFT looks like symbolically? Well, now you can.
```julia
using FFTA, Symbolics
N = 16
@variables x_a[1:N]::Complex
x = collect(x_a)
y = simplify.(fft(x))
```
Now, if you have a signal afterward that you want to substitute in, you can call `map(y_el -> substitute(y_el, Dict(x .=> signal)), y)`. Make no mistake-- it's almost certainly more efficient to just plug your type into `FFTA.fft` than using substitution. But this is an example of how `FFTA` integrates wonderfully and gracefully with the Julia ecosystem. If you want high precision FFTs, use `Complex{BigFloat}`. If you want to use an `SVector` from `StaticArrays` because your data is small, then use that! If you want to use `SizedArray{Complex{BigFloat}}`, be my guest. These are opportunities that won't be provided to you in almost any other package out there.

As of this commit, you can do
```julia
julia> import FFTA, FFTW

julia> N = 64

julia> @btime FFTA.fft(x) setup=(x = @SVector rand(N));
  698.611 ns (8 allocations: 2.11 KiB)

julia> @btime FFTW.fft(x) setup=(x = @SVector rand(N));
  5.433 μs (34 allocations: 4.70 KiB)
```
It's painfully obvious that this example is cherry-picked. Nonetheless, the user can finally take the speedups so much of the Julia community has worked so hard on and propogate them into the FFT.
