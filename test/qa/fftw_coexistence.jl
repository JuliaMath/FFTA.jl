using Test

# Loading FFTW.jl makes its methods take over every `fft`/`rfft` call in the
# process, so the coexistence checks run in a separate Julia process to leave
# the rest of the test suite exercising FFTA.
@testset "coexistence with FFTW.jl (subprocess)" begin
    body = joinpath(@__DIR__, "fftw_coexistence_body.jl")
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $body`
    @test success(pipeline(cmd; stdout = stdout, stderr = stderr))
end
