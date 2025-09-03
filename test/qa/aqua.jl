using FFTA, Test
import Aqua

@testset "Aqua" begin
  Aqua.test_all(
    FFTA;
    # This type piracy is caused by the problematic design of AbstractFFTs.jl
    # Ref https://github.com/JuliaMath/AbstractFFTs.jl/issues/32
    piracies = (; treat_as_own = [plan_bfft, plan_brfft, plan_fft, plan_rfft]),
  )
end
