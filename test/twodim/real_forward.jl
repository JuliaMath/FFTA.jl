using FFTA, Test
test_nums = [8, 11, 15, 16, 27, 100]
@testset " forward" begin 
    for N in test_nums
        x = ones(Float64, N, N)
        y = rfft(x)
        y_ref = 0*y
        y_ref[1] = length(x)
        @test y ≈ y_ref
    end
end