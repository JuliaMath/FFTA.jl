fft(X::AbstractArray{<:Complex,N}, region::NTuple{D,Int} where D = ntuple(identity, N)) where N = plan_fft(X, region) * X
fft(X::AbstractArray{<:Any,N}, region::NTuple{D,Int} where D = ntuple(identity, N)) where N = fft(complex(X), region)

bfft(X::AbstractArray{<:Complex,N}, region::NTuple{D,Int} where D = ntuple(identity, N)) where N = plan_bfft(X, region) * X

ifft(X::AbstractArray{<:Complex,N}, region::NTuple{D,Int} where D = ntuple(identity, N)) where N = bfft(X, region) / mapreduce(Base.Fix1(size, X), *, region; init=1)

rfft(X::AbstractArray{<:Real,N}, region::NTuple{D,Int} where D = ntuple(identity, N)) where N = plan_rfft(X, region) * X

brfft(X::AbstractArray{<:Complex,N}, len::Int, region::NTuple{D,Int} where D = ntuple(identity, N)) where N = plan_brfft(X, len, region) * X

function irfft(X::AbstractArray{<:Complex,N}, len::Int, region::NTuple{D,Int} where D = ntuple(identity, N)) where N
    Y = brfft(X, len, region)
    Y ./= mapreduce(Base.Fix1(size, Y), *, region; init=1)
    return Y
end
