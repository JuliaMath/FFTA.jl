# Thread-scaling probe for the hybrid P/E question.
# WSL2 exposes 22 uniform vCPUs (6P x2 SMT + 8E + 2LP-E on this part) with no
# way to tell them apart or pin to them, so instead of pinning we measure the
# shape of the scaling curve and the run-to-run spread.
#
# A single-worker FFTA plan is re-measured inside every trial round as a
# control: it cannot be affected by thread placement, so its spread is the
# host's baseline noise. Widening beyond that at high thread counts is
# attributable to placement, not to the box being generally noisy.
using FFTA, FFTW, AbstractFFTs, LinearAlgebra, Printf

function tmin(f; reps=20)
    f(); minimum(@elapsed(f()) for _ in 1:reps)
end
spread(v) = 100 * (maximum(v) - minimum(v)) / minimum(v)

const SZ = (1024, 256)
x = randn(ComplexF64, SZ)

# the control: one worker, never varies with nt
p1 = invoke(AbstractFFTs.plan_fft, Tuple{AbstractArray{ComplexF64,2},Any}, x, 1; num_threads = 1)
y1 = similar(x)

println("## Thread scaling on a hybrid P/E part (WSL2, cannot pin)\n")
println("Batched fft along dims=1, $(SZ[1])x$(SZ[2]) ComplexF64. 3 independent trials per point.")
println("`ctrl` = FFTA with 1 worker, re-measured inside every trial round as a noise control.\n")
println("| threads | FFTW (ms) | speedup | FFTA (ms) | speedup | FFTW spread | FFTA spread | ctrl spread |")
println("|--:|--:|--:|--:|--:|--:|--:|--:|")

base_w = nothing; base_a = nothing
allctrl = Float64[]
for nt in (1, 2, 4, 6, 8, 11)
    FFTW.set_num_threads(nt)
    pw = FFTW.plan_fft(x, 1); yw = similar(x)
    pa = invoke(AbstractFFTs.plan_fft, Tuple{AbstractArray{ComplexF64,2},Any}, x, 1; num_threads = nt)
    ya = similar(x)
    tw = Float64[]; ta = Float64[]; tc = Float64[]
    for _ in 1:3
        push!(tw, tmin(() -> mul!(yw, pw, x)))
        push!(ta, tmin(() -> mul!(ya, pa, x)))
        push!(tc, tmin(() -> mul!(y1, p1, x)))     # control, same round
    end
    append!(allctrl, tc)
    w, a = minimum(tw), minimum(ta)
    global base_w, base_a
    if nt == 1; base_w = w; base_a = a; end
    @printf("| %d | %.3f | %.2f× | %.3f | %.2f× | %.1f%% | %.1f%% | %.1f%% |\n",
            nt, w*1e3, base_w/w, a*1e3, base_a/a, spread(tw), spread(ta), spread(tc))
end
@printf("\nPooled control spread over the whole probe (%d samples): **%.1f%%** — this is the\n",
        length(allctrl), spread(allctrl))
println("host's baseline noise floor for a single-worker transform on this box.")
println()
println("A near-linear curve to 6 that flattens after it is the signature of SMT/E-core")
println("placement. FFTA/FFTW spread materially above the control at 8 and 11 threads is")
println("the signature of the scheduler migrating workers between core types.")
