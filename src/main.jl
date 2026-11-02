fft(X::AbstractArray{<:Complex}, region::RegionTypes = 1:ndims(X); kwargs...) = plan_fft(X, region; kwargs...) * X
fft(X::AbstractArray, region::RegionTypes = 1:ndims(X); kwargs...) = fft(complex(X), region; kwargs...)

bfft(X::AbstractArray{<:Complex}, region::RegionTypes = 1:ndims(X); kwargs...) = plan_bfft(X, region; kwargs...) * X

ifft(X::AbstractArray{<:Complex}, region::RegionTypes = 1:ndims(X); kwargs...) =
    bfft(X, region; kwargs...) / mapreduce(Base.Fix1(size, X), *, region; init=1)

rfft(X::AbstractArray{<:Real}, region::RegionTypes = 1:ndims(X); kwargs...) = plan_rfft(X, region; kwargs...) * X

brfft(X::AbstractArray{<:Complex}, len::Int, region::RegionTypes = 1:ndims(X); kwargs...) = plan_brfft(X, len, region; kwargs...) * X

function irfft(X::AbstractArray{<:Complex}, len::Int, region::RegionTypes = 1:ndims(X); kwargs...)
    Y = brfft(X, len, region; kwargs...)
    Y ./= mapreduce(Base.Fix1(size, Y), *, region; init=1)
    return Y
end
