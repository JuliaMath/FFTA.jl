using NSFFT

for N in [8, 11, 15, 100]
    x = zeros(ComplexF64, N)
    x[1] = 1
    y = NSFFT.fft(x)
    @test y ≈ ones(size(x))
end