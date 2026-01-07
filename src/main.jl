fft(X::AbstractArray{<:Complex}, region::Union{Int,AbstractVector} = 1:ndims(X)) = plan_fft(X, region) * X
fft(X::AbstractArray, region::Union{Int,AbstractVector} = 1:ndims(X)) = fft(complex(X), region)

bfft(X::AbstractArray{<:Complex}, region::Union{Int,AbstractVector} = 1:ndims(X)) = plan_bfft(X, region) * X

ifft(X::AbstractArray{<:Complex}, region::Union{Int,AbstractVector} = 1:ndims(X)) = bfft(X, region) / mapreduce(Base.Fix1(size, X), *, region; init=1)

rfft(X::AbstractArray{<:Real}, region::Union{Int,AbstractVector} = 1:ndims(X)) = plan_rfft(X, region) * X

brfft(X::AbstractArray{<:Complex}, len::Int, region::Union{Int,AbstractVector} = 1:ndims(X)) = plan_brfft(X, len, region) * X

function irfft(X::AbstractArray{<:Complex}, len::Int, region::Union{Int,AbstractVector} = 1:ndims(X))
    Y = brfft(X, len, region)
    Y ./= mapreduce(Base.Fix1(size, Y), *, region; init=1)
    return Y
end
