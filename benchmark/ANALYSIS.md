# FFTA.jl performance analysis: implementation gaps vs. structural limits

This document separates the reasons FFTA.jl is slower than FFTW into
**implementation gaps** (fixable inside FFTA's current design) and
**structural limits** (consequences of being a pure-Julia, recursive,
size-generic implementation), and gives a realistic ceiling per size class.
All measurements are from `benchmark/suite.jl` and a set of targeted
experiments on the machine described in [`REPORT.md`](REPORT.md)
(aarch64 Neoverse-N1, Julia 1.12.6, FFTW 3.3.11 single-threaded, `ComplexF64`
unless stated). Ratios are *FFTA time / FFTW time* for planned execution.

## 1. Executive summary

| size class (`ComplexF64` unless noted) | today, vs FFTW `ESTIMATE` | dominant cause | fixable? | realistic ceiling |
|:--|--:|:--|:--|--:|
| 2^k, in cache (n ≤ 2^16) | 2.3–4.4× (7× at 32) | scalar radix-4 recursing to 2/4-point base cases; twiddle recurrence seeded per call | mostly | 1.3–2× |
| 2^k, memory-bound (n ≥ 2^18) | 1.0–1.5× (2.5–3× vs `MEASURE`) | FFTW `ESTIMATE` is itself 2.2–2.8× off FFTW's best here; FFTA has no cache-blocked large-n algorithm | yes | 1.2–1.5× vs `MEASURE` |
| smooth 2^a3^b5^c7^d | 2.6–15× (geomean 7.3×) | O(n²) leaves for 5 and 7 that recompute `sincospi` twiddles per row per execution; no codelets | yes | 1.5–3× |
| primes < 73 | 6–11× | O(n²) DFT with per-execution twiddles | yes | 1.5–3× |
| primes ≥ 73 | 3–9× | Bluestein padded to 2^k, chirp + its FFT recomputed and re-allocated on every call | yes | 1.5–3× |
| prime × small factor | 3.4–11× | as above, inside a composite | yes | 1.5–3× |
| `rfft` 1D | 1.4–7× pow2, 3–25× smooth | already half-size complex; inherits the complex gap + `*`-only API (allocates) | yes | 1.5–2.5× |
| 2D / 3D `fft` | 1.8–15× (worst at small sizes) | per-pencil overhead, copies through per-call buffers, no threading | yes | 1.3–2× (1 thread) |
| 2D `rfft` | 6–25× | full complex transform + copies | yes | 1.5–2.5× |
| batched `dims` | 1.2–5.8× (`fft`), 1.4–10× (`rfft`) | `rfft` along `dims` goes through `mapslices` | yes | 1.3–2× |
| `Float32` | no faster than `Float64` | no SIMD anywhere | partly | see §5.3 |
| vs FFTW with 8 threads | 7–20× | no threading | yes for ND/batched | ~1.5–2× |

The one-line version: **nothing in FFTA's algorithmic approach forces a 10×
gap**; every ≥ 4× class above is an implementation gap. What is structural is
the last 1.3–2× on cache-resident sizes, which is the price of not having a
SIMD codelet library, and a bounded compile-latency cost if FFTA generates
its own.

## 2. What FFTA does today

* **Planning** (`CallGraph{T}(n)`, `src/callgraph.jl`): `n` is recursively
  split. Powers of 2 and 3 become single leaf nodes (`POW2RADIX4_FFT`,
  `POW3_FFT`); primes become `DFT` (n < 73) or `BLUESTEIN` (n ≥ 73) leaves;
  anything else is a `COMPOSITE_FFT` node split into `N1 × N2`, where `N1` is
  the full power of 2 or 3 if present, otherwise the product of prime
  factors closest to √n. Each composite node owns an `n`-element workspace.
  Planning is cheap (microseconds) and never measures anything.
* **Execution** (`src/algos.jl`): `fft!` dispatches on a node-type enum at
  run time. `COMPOSITE_FFT` is textbook Cooley–Tukey: N1 sub-transforms of
  size N2 into the workspace, twiddle multiply, N2 sub-transforms of size N1
  into the output. `POW2RADIX4_FFT` is a recursive decimation-in-time radix-4
  with 2- and 4-point base cases. `POW3_FFT` is recursive radix-3. `DFT` is
  the O(n²) sum. `BLUESTEIN` pads to the next power of two ≥ 2n−1 and does
  three `POW2RADIX4` transforms.
* **Twiddles**: nowhere stored. Every kernel seeds Singleton's recurrence
  with `singleton_params(...)`, which calls `sincospi`, then steps it. The
  `DFT` leaf does this once per output row; the composite node once per
  `j1`; the radix-4 and radix-3 kernels three/two times per recursion level.
* **Real transforms**: 1D even-length `rfft` uses the standard
  half-length complex trick (two real sequences packed into one complex
  transform, then a butterfly). Odd lengths and 2D real transforms run a
  full complex transform and discard half. `rfft`/`irfft` plans implement
  only `*`, not `mul!`; along `dims` of an N-d array they use `mapslices`.
* **Multidimensional**: transforms along each dimension by copying each
  pencil into a contiguous buffer, transforming, and copying back. The two
  buffers are allocated per call.
* **No threading, no SIMD, no in-place (`plan_fft!`) plans.**

## 3. What FFTW does that matters here

Looking at the plans FFTW's `ESTIMATE` planner actually chose on this
machine (`FFTW.plan_fft(x)` prints them):

* Sizes ≤ 64 (and many "twiddle" steps) are single **codelets**:
  `n1fv_5_neon`, `n1fv_13_neon`, `t3fv_32_neon`, `n1fv_128_neon`, …. These
  are straight-line C functions generated offline by `genfft`, with all
  twiddles folded to constants, common sub-expressions eliminated, and
  operations scheduled for the target's register file. The `v` variants
  operate on two transforms (or two halves) at once so that a 128-bit NEON
  register holds a full `ComplexF64`.
* Power-of-two sizes are one or two radix-32 twiddle passes on top of a
  64- or 128-point codelet: `4096 = t3fv_32 ∘ n1fv_128`. There is no
  recursion below 64.
* Composite sizes use the same machinery with mixed radices:
  `1000 = t3fv_25 ∘ t3fv_5 ∘ n2fv_8`, `720 = t3fv_20 ∘ t1fv_6 ∘ n2fv_6`.
* Primes: **Rader** for 61, 73 and 65537 (n−1 is smooth), **Bluestein** for
  127, 1009, 4099 — padded to a *smooth* length (1009 → 2025 = 3⁴·5²,
  4099 → 8640), not a power of two, and with the chirp's transform stored in
  the plan.
* `rfft` uses dedicated real-input codelets (`r2cf_*`, `hc2cfdftv_*`), i.e.
  it does not go through a complex transform at all.
* With `MEASURE`, the planner additionally *times* candidate decompositions.
  On this machine that buys FFTW 1.1–1.8× over `ESTIMATE` for cache-resident
  powers of two and **2.2–2.8× for n ≥ 2^19** (see `REPORT.md`,
  "FFTW.MEASURE" table): the measured plans pick decompositions that make
  fewer passes over DRAM. Below the cache size most of FFTW's advantage is
  the codelet library and the heuristics, not the measuring; above it, the
  measuring matters.

## 4. Implementation gaps (fixable)

Each item quantifies the cost with a targeted experiment and says what the
fix looks like. Items are ordered by expected impact on the benchmark suite.

### 4.1 Twiddle factors are recomputed on every execution

`singleton_params` costs a `sincospi` (tens of ns). The `DFT` leaf calls it
n−1 times per execution; for n = 5 that is the entire cost of the transform:

| n | FFTA today | same loop with a precomputed table | FFTW codelet |
|--:|--:|--:|--:|
| 5 | 172 ns | 37 ns | 23 ns |
| 7 | 292 ns | 74 ns | 28 ns |
| 11 | 612 ns | 191 ns | 42 ns |
| 13 | 807 ns | 270 ns | 50 ns |
| 31 | 3.7 µs | 1.6 µs | 0.39 µs |
| 61 | 13.3 µs | 6.4 µs | 1.4 µs |

The composite kernel likewise seeds a recurrence per `j1` and the radix-4/3
kernels per recursion level. For n = 1000 (= 8 · 5 · 5 · 5) roughly 1 800
`sincospi` evaluations happen per execution, on a transform FFTW finishes in
10 µs. The recurrence also costs accuracy: in `Float32`, FFTA's error against a
`Float64` reference grows from 10 ulp at 2^16 to ~1000 ulp (1.3e-4
relative) at 2^22, while FFTW stays at ~1.5 ulp (the upstream accuracy test
stops at 2^18). `Float64` stays within 3e-14. **Fix:** store per-node
twiddle tables in the `CallGraph` at plan time (this is what a plan is
for), computed directly with `sincospi` so they are correctly rounded.
Memory cost is ≤ n complex numbers per node, i.e. comparable to the
workspace FFTA already allocates. Expected gain: 2–5× on every size with a
prime factor ≥ 5, ~10–20% on powers of 2/3, and `Float32` accuracy back to
a few ulp.

### 4.2 Small transforms are O(n²) loops, not codelets

Even with tables, the O(n²) leaf is 1.6× (n = 5) to 4.6× (n = 61) slower
than FFTW's codelet, and every composite size in the smooth class bottoms
out in such leaves. The experiment in §5.2 shows that Julia can generate
straight-line codelets that match FFTW within 1.0–1.7× for n ≤ 256.
**Fix:** `@generated` codelets keyed on `Val{n}` for the small prime leaves
(5, 7, 11, 13) and for the power-of-two base cases (8/16/32/64), selected by
the planner. Radix-5 and radix-7 *butterflies* (the analogue of
`fft_pow3!`) would remove the composite-node overhead for 5^a and 7^b
factors entirely (upstream issue #105 discusses this).

### 4.3 The power-of-two kernel recurses to 2- and 4-point base cases

FFTA is 2.3–4.4× slower than FFTW for cache-resident powers of two (7× at
n = 32), and the ratio is roughly flat across 2^4…2^14, i.e. it is
per-butterfly cost, not memory traffic. Odd powers of two (whose radix-4
recursion ends in a radix-2 step) are consistently worse than even ones
(4.2–4.4× vs 3.1×). Causes, in decreasing order: (i) the recursion bottoms out at
n = 2 or 4, so half the work is function-call and index arithmetic;
(ii) three recurrence steps per butterfly instead of loads from a table;
(iii) no SIMD (see §5.3). At n ≥ 2^18 both libraries are memory-bound and
the ratio drops to 1.0–1.5× against `ESTIMATE` plans — but see §5.4. **Fix:** stop the recursion at a 16/32/64-point
straight-line base case (§5.2 measures 1.0–1.1× of FFTW at 16 and 64),
load twiddles from the plan, and consider a radix-8 pass to reduce passes
over memory for the large sizes.

### 4.4 Bluestein allocates and recomputes on every call

For n = 1009 (207 µs vs FFTW's 54 µs):

| component | time |
|:--|--:|
| `prealloc_blue` (3 × 2048-element allocations + chirp via `cispi`) | 23 µs |
| three 2048-point radix-4 transforms | 3 × 56 µs |
| everything else | ~15 µs |

One of the three transforms is the chirp's, which depends only on n and
direction. **Fix:** precompute chirp, its transform and the scratch
buffers in the plan (saves ≈ 40%); pad to the smallest *smooth* length
≥ 2n−1 instead of a power of two once the composite path is fast (FFTW pads
1009 to 2025 where FFTA uses 2048 — similar — but pads 4099 to 8640 where
FFTA uses 16384, and 65537 to a Rader plan where FFTA uses 262144). Expected: 2–3× on primes ≥ 73 and
on the "awkward" class. In the composite path (`fft_composite!`) the
scratch *is* hoisted, but re-allocated per outer call, so batched/ND
transforms with a Bluestein factor still allocate ~n × 48 bytes per pencil.

### 4.5 The Bluestein cutoff (73) is too high, and Rader is absent

The O(n²) leaf at n = 61 costs 13 µs (6.4 µs with tables); Bluestein at 73
costs 16 µs. FFTW does both in 1.4–1.7 µs using Rader's algorithm, which
maps a prime-length DFT onto a (n−1)-length convolution — n−1 is even and
usually smooth. **Fix:** after 4.1/4.2, re-tune the cutoff (likely ~20–30),
and add Rader for primes whose n−1 is smooth. Rader is not required for
the ceiling in §6 but is the difference between 2–3× and ~1.5× on primes.

### 4.6 Real transforms

The 1D even-`rfft` path is structurally right (half-length complex +
butterfly): at n = 4096 the complex half-transform is 85% of the time and
the butterfly plus output allocation 15%. So `rfft` mostly inherits the
complex gap, plus: (i) plans implement `*` only, so there is no
zero-allocation path and `mul!`-based consumers (DSP.jl's `fftfilt`,
periodograms) cannot use a preallocated output; (ii) odd lengths and 2D run
a full complex transform (2× the work); (iii) along `dims` of a matrix the
`mapslices` path allocates per column and is 10× slower than FFTW, versus
3–6× for the complex `dims` path. **Fix:** `mul!` for real plans, a
strided/batched real path that reuses the complex `fft_along_dim!` loop,
and the odd-length / 2D cases via the same half-length trick. One thing the
`mapslices` path did right, discovered only after replacing it: its copy of
each pencil into a contiguous vector was an unlabelled copy-in
optimisation for strided pencils (a `dims=2` pencil of a 64-row matrix
strides a cache line per element), worth 1.2–1.3× on x86-64; the
replacement keeps a copy-in for non-unit-stride pencils, from plan-owned
buffers.

### 4.7 Multidimensional execution

2D 256×256: 2.2× FFTW, of which ~80% is the 512 pencil transforms (i.e. the
1D gap) and ~20% is the copy-in/copy-out through `ibuf`/`obuf`, which are
also allocated per call (8 KiB for 256², up to 2 × n × 16 bytes). At small
sizes the per-pencil overhead dominates: 8×8 is 13×, 8×8×8 is 15×, 32×32
is 11× FFTW. Composite 2D sizes inherit the 1D composite gap (1000×1000:
10×, 720×480: 8×), and 2D `rfft` — a full complex transform plus copies —
is 6–25×. Sizes with a Bluestein factor allocate 6–10 MiB per call
(1009×64, 127×257). FFTW
transforms strided pencils in place (its "vrank" plans) and batches
several columns per codelet call. **Fix:** allocate buffers in the plan,
transform contiguous dimension-1 pencils directly (the kernels already take
strides; `fft_along_dim!` copies even when the pencil is contiguous), and
thread across pencils (§5.6).

### 4.8 API and dispatch issues found while benchmarking

* Loading FFTW.jl and FFTA.jl together makes `plan_rfft(::Vector{Float64},
  ::Int)` — and hence `rfft(x)` — a **method ambiguity error**, because
  FFTA annotates `region::RegionTypes` while FFTW annotates the array as
  `StridedArray`. Neither package is "more specific". This breaks any
  environment that has both loaded (e.g. a DSP.jl user trying FFTA).
* No `mul!` for real plans (only `*`), which is what DSP.jl uses at 10 of
  its 13 plan-execution sites; `mul!` into a `SubArray` output fails even
  for complex plans.
* No in-place plans (`plan_fft!`, `plan_bfft!`, hence `ifft!`).
* `inv(p)` throws a `TypeError` rather than working or giving a clean
  `MethodError`: FFTA defines no `AbstractFFTs.plan_inv` method, and the
  dummy `pinv::FFTAInvPlan` field makes `AbstractFFTs.inv`'s
  `pinv_type(p)` resolve to `Union{}`. Consequently `p \ x` and `ldiv!`
  are unavailable too.
* `plan_rfft`/`plan_brfft` for 3D arrays throw (DSP.jl's `conv` is tested
  for N = 3).
* The region argument is not inferred (`plan_fft(x, 2)` returns a plan
  whose type depends on `region`'s run-time type; upstream #78/#91).
* `Float32` transforms are exactly as slow as `Float64` (no SIMD), while
  FFTW is 1.75× faster in single precision.

## 5. Structural limits

### 5.1 No measuring planner

FFTA chooses a factorization by a fixed rule. FFTW's `ESTIMATE` mode also
does not measure; on this machine `MEASURE` improves `ESTIMATE` by 0–25% for
powers of two. FFTA could adopt FFTW's heuristics (radix-32/16 first,
codelet sizes at the bottom) without timing anything. A timing planner is
implementable in Julia (it is just a loop over candidate call graphs at plan
time), but its value is small and its plan-time cost large (FFTW `MEASURE`
takes seconds at 2^22). **Verdict: not a real limit.** Plan creation is
actually an FFTA advantage: FFTA plans in 0.1–30 µs, FFTW `ESTIMATE` in
2 µs–190 ms (large primes are slow: 10 ms at n = 120 779, 190 ms at
1 727 797) and `MEASURE` in milliseconds to seconds; consumers that plan
per call (DSP.jl does this in several places) benefit, and one-shot
`fft(x)` on large primes is within 1.5–2× of FFTW today.

### 5.2 Codelets and Julia's compilation model

FFTW ships ~150 pre-generated codelets per SIMD flavour. In Julia the
equivalent is an `@generated` function keyed on `Val{n}` that emits
straight-line code. A 40-line generator (radix-2 DIT, every intermediate in
its own SSA variable, twiddles folded to constants) gives, on this machine:

| n | statements | compile (first call) | codelet | FFTA today | FFTW | codelet / FFTW |
|--:|--:|--:|--:|--:|--:|--:|
| 8 | 52 | 0.02 s | 16 ns | 114 ns | 26 ns | 0.61 |
| 16 | 128 | 0.05 s | 43 ns | 148 ns | 41 ns | 1.05 |
| 32 | 304 | 0.17 s | 128 ns | 584 ns | 80 ns | 1.61 |
| 64 | 704 | 0.62 s | 341 ns | 817 ns | 309 ns | 1.10 |
| 128 | 1 600 | 2.5 s | 875 ns | 2.75 µs | 644 ns | 1.36 |
| 256 | 3 584 | 6.8 s | 2.27 µs | 4.05 µs | 1.31 µs | 1.73 |
| 5 (O(n²)) | — | 0.03 s | 31 ns | 179 ns | 22 ns | 1.39 |
| 7 (O(n²)) | — | 0.05 s | 68 ns | 301 ns | 26 ns | 2.60 |
| 13 (O(n²)) | — | 0.18 s | 264 ns | 819 ns | 48 ns | 5.5 |

So: **codelets up to 64 points are within 1.0–1.6× of FFTW and compile in
under a second**; LLVM's compile time grows super-linearly beyond that
(a naive O(n²) 64-point unroll took 39 s), so the codelet set must be small
and fixed — **2–64 for powers of two**, plus 3, 5, 7, 9, 11, 13, 25 — with
the generic recursion above it. The upper end is set by register pressure,
not compile time: on x86-64 (AVX2, 16 architectural vector registers
against NEON's 32) the same generator reproduces the aarch64 ratios up to
64 points (0.51 / 1.11 / 1.54 / 1.13 at 8/16/32/64) but spills earlier and
about twice as often beyond, giving 1.70× at 128 and **3.62× at 256**
(aarch64: 1.36× and 1.73×) although it compiles them *faster* (4.9 s vs
6.8 s at 256). 128- and 256-point codelets are therefore pessimisations on
x86-64 while looking tolerable on aarch64; #135 stops at 64. The compile cost is paid once per (n, T) per
Julia session unless the package precompiles them with `PrecompileTools`
(then it is paid once at package install, ~10–30 s for the set above, and
the package's precompile cache grows by a few MB). This is the structural
cost of "planner-like behaviour" in Julia: **a fixed codelet set with
precompilation, not per-size specialisation** — an arbitrary size like 1009
must still be built from the fixed set at run time, exactly as FFTW does.
The O(n²) codelets for 11 and 13 remain 5× off because FFTW's are
Winograd/Rader-derived with far fewer multiplications; those algorithms can
be generated too but are more work.

### 5.3 SIMD on interleaved complex data

On this CPU a NEON register is 128 bits = one `ComplexF64`, so SIMD gains
for double precision must come from operating on *two independent
butterflies* per instruction (FFTW's `*v_*` codelets) or from split
real/imaginary storage. LLVM does not do this transformation for scalar
straight-line code: the `ComplexF32` codelet in §5.2 is no faster than the
`ComplexF64` one, while FFTW's is 1.6× faster. Reaching it in Julia means
explicit `SIMD.jl`/`VectorizationBase` vectors in the codelet generator, or
vectorising across pencils in ND/batched transforms (each lane a different
column — the easy and big win for the downstream use case). On AVX2/AVX-512
x86-64 the same argument applies with 2–4 `ComplexF64` per register, so the
FFTA/FFTW gap on x86-64 was expected to be *larger* than the aarch64 numbers
here until this is done. **Verdict: the last ~1.3–2× on cache-resident
sizes is structural until FFTA has an explicitly vectorised kernel
generator; it is not blocked by the language.**

**What the x86-64 companion run says (see §5.3a).** The gap *is* larger on
AVX2, by 1.22× overall — but the vector-width argument explains only the
power-of-two part of it (1.21×). The bulk of the x86 penalty sits on
composite sizes with factors of 5 and 7 (up to 3.8×), i.e. on the kernels
that are *not* the radix-4 one. So the first x86-specific lever is not the
vectorised generator but radix-5/7 butterflies (work item 3), and the
ceiling table below needs a per-class x86 adjustment rather than a uniform
one.

### 5.3a The x86-64 penalty is mostly a radix-5/7 penalty

The same suite on an AVX2 machine (Intel Core Ultra 7 165H, no AVX-512;
`benchmark/x86-64/` in the report) gives, per size class, the ratio of the
x86-64 FFTA/FFTW gap to the aarch64 one (geometric mean over the four
element types, FFTA 0.3.1):

| class | x86-64 gap / aarch64 gap |
|:--|--:|
| smooth composites | **1.60×** |
| batched `dims=1` | 1.23× |
| powers of two | 1.21× |
| primes / prime × small | 1.07× |
| 2D | 1.06× |
| 3D | 0.96× |
| batched `dims=2` | 0.78× |

Within the smooth class the penalty tracks the factor content
(`ComplexF64`): 98304 = 2^15·3 → 1.21×; 54 = 2·3^3 → 1.33×; 25 = 5^2 →
2.10×; 46305 = 3^3·5·7^3 → 2.39×; 10^6 = 2^6·5^6 → 2.58×;
441000 = 2^3·3^2·5^3·7^2 → **3.77×**. Factors of 2 and 3 cost ~1.2–1.3×
more on x86-64 than on aarch64; factors of 5 and 7 cost 2.1–3.8× more. The
Bluestein padding calibration (§4.4, `benchmark/x86-64/calib_*.md`) found
the same thing independently — a 3-smooth length costs 2.0–3.1× a power of
two per element on x86-64 against 1.3–2.3× on aarch64. Two measurements,
one conclusion: FFTA's non-radix-4 kernels (the O(n²) DFT leaves for 5 and
7, the composite step, the radix-3 kernel) degrade disproportionately on
x86-64, presumably because the radix-4 kernel is the one LLVM vectorises
best and the wider registers make the difference larger. The batched
`dims=2` class is the one place x86-64 is relatively *faster*; it is also
the shape where the merged optimisations gain least on aarch64 and regress
on x86-64 (§6), i.e. that shape behaves qualitatively differently between
the two architectures in both directions.

### 5.4 Memory-bound regime

For n ≥ 2^18 (4 MiB of `ComplexF64`) FFTA is within 1.0–1.5× of FFTW's
`ESTIMATE` plans, which is misleading: FFTW's `MEASURE` plans are 2.2–2.8×
faster than its `ESTIMATE` plans in this regime, so the gap to the best
FFTW plan is ~2.5–3×. What `MEASURE` finds is a decomposition with fewer
passes over DRAM (large radices, in-place transposes, buffered
sub-transforms — the "four-step"/cache-blocked family). FFTA's depth-first
radix-4 recursion is cache-oblivious to a degree but its twiddle pass at
each level touches the whole array. Closing this needs a cache-blocked
large-n path (transform as an n₁×n₂ matrix with a transpose step and
contiguous sub-transforms) — algorithmically standard, and independent of
SIMD. This is fixable but is a new code path rather than a tweak.

### 5.5 Primes

Asymptotically both approaches are O(n log n) via Bluestein; FFTW's
constant is better because of smooth padding and a stored chirp transform,
and Rader wins when n−1 is smooth. After §4.4/4.5, primes should sit at
2–3× the cost of a same-size smooth transform in both libraries, i.e. a
1.5–3× ratio. **No structural limit** beyond §5.3.

### 5.6 Threading

FFTW threads *inside* one transform (splitting radix passes). FFTA has no
threading, and its `CallGraph` owns a single workspace, so a plan is not
safe to share between threads. Threading across independent pencils in ND
and batched transforms is straightforward (one workspace per thread) and
gives near-linear speedups for the dominant signal-processing pattern
(`fft(X, 1)` over many columns). Threading within a single 1D transform
needs a breadth-first pass structure, which is a bigger restructuring.
FFTW with 8 threads is 5–12× faster than with one on n ≥ 2^16, 2D ≥ 256²
and batched matrices (see `REPORT.md`, threading section), leaving FFTA
7–20× behind in those configurations. For downstream signal-processing
users the batched and ND cases are the ones that matter, and those are the
easy ones to thread.

### 5.7 Generic element types

FFTA's selling point — `Complex{BigFloat}`, dual numbers, symbolic
elements — is preserved by every fix above as long as the codelet path is
gated on `isbitstype` / `T <: Union{Float32,Float64}` and the generic
recursion remains the fallback. This constrains the design (two code paths)
but not the ceiling.

## 6. Realistic ceiling

Combining §4 and §5 (single-threaded, this machine; x86-64 with AVX2/512
will show larger ratios until §5.3 is addressed):

| class | today (vs `ESTIMATE`) | after §4.1–4.7 (no SIMD) | with vectorised codelets |
|:--|--:|--:|--:|
| 2^k in cache | 2.3–4.4× | 1.5–2× | 1.1–1.5× |
| 2^k memory-bound (vs `MEASURE`) | 2.5–3× | 1.5–2× (needs the blocked path of §5.4) | 1.2–1.5× |
| smooth composites | 2.6–15× | 1.5–3× | 1.2–2× |
| primes / awkward | 3–11× | 2–3× | 1.5–2.5× |
| `rfft` | 1.4–25× | 1.5–2.5× | 1.2–2× |
| 2D/3D/batched, 1 thread | 1.2–25× | 1.3–2× | 1.1–1.5× |
| 2D/3D/batched vs FFTW 8 threads | 7–20× | ~1.5–2× with pencil threading | ~1.2–1.5× |

"1.5–3× of FFTW across the sweep" is achievable without SIMD work; getting
under 1.5× everywhere needs the vectorised generator. A 10× gap is
**nowhere unavoidable**; the sizes where it exists today (5- and 7-smooth
composites, primes) are exactly the ones fixed by storing twiddles and
codelets, i.e. by making the plan carry the work a plan is supposed to
carry.

**x86-64 adjustment (from §5.3a).** The table above is for aarch64. On
AVX2 x86-64 multiply the power-of-two, prime and ND rows by ~1.1–1.2×, but
the smooth-composite row by 1.6× on average and by 2–4× for sizes rich in
factors of 5 and 7 — until radix-5/7 butterflies replace the O(n²) leaves
(work item 3). On AVX-512 parts both factors are expected to be larger
still (not measured).

## 7. Recommendation for downstream projects today

Until the fixes land: FFTA is a reasonable substitute for 1D power-of-two
transforms (within 2.3–4.4× in cache, 1.0–1.5× of FFTW-`ESTIMATE` above
2^18); it is a poor substitute (5–15×) for sizes with factors of 5 or 7,
for primes, for small 2D/3D arrays and for 2D `rfft`, and its `rfft` along
`dims` should be avoided. Since signal-processing consumers overwhelmingly use
power-of-two `nfft` (Welch, periodograms, `fftfilt`) the practical cost is
the 2–5× on `rfft`, plus the ambiguity in §4.8 which must be fixed before
FFTA and FFTW can coexist in one environment at all.

## 8. Work list (one PR each)

1. Store twiddle tables and Bluestein chirp/scratch in the plan (§4.1, §4.4).
2. `mul!` for real plans; zero-allocation `rfft`/`irfft`; `dims` path for
   real transforms without `mapslices` (§4.6).
3. Straight-line base cases for the power-of-two kernel (16/32/64) and
   codelets for 3/5/7 leaves; radix-5/7 butterflies (§4.2, §4.3). **On
   x86-64 the radix-5/7 butterflies are the single largest
   architecture-specific gap in the sweep (§5.3a) and are worth more than
   item 7 there.**
4. Plan-owned ND buffers, contiguous-pencil fast path, threading across
   pencils (§4.7, §5.6).
5. Smooth-length Bluestein padding; retune the cutoff; Rader (§4.4, §4.5).
6. Method-ambiguity fix with FFTW, `plan_fft!`, 3D `rfft` (§4.8).
7. Vectorised codelet generator (§5.3) — larger, separate design discussion.

## 9. Where things stand after the work list (measured)

Items 1–6 of §8 have been implemented (twiddle tables and Bluestein data in
the plan; real-plan `mul!`, zero-allocation and N-d real paths; 8–64-point
straight-line base cases; plan-owned per-thread workers with threading
across pencils; a re-calibrated Bluestein cutoff with a dormant 3-smooth
pad chooser; inverse and in-place plans; the FFTW.jl ambiguity fix). With
all of them merged, the same suite on the same aarch64 machine gives, over
the 510 cases comparable with the baseline, a geometric-mean FFTA speedup of
**2.57× with no case slower by more than 5%**; on an AVX2 x86-64 machine
**3.24×** with three x86-specific regressions (wide `64×N` `rfft` along
`dims=2`, 1.16–1.34×, alongside the removal of 16–66 MiB of allocation per
call). FFTA/FFTW after (aarch64, single thread, FFTW `ESTIMATE`):

| type | pow2 | smooth | prime | awkward | 2D | 3D | batched dim=1 | batched dim=2 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| ComplexF64 fft | 1.19× (0.8–1.7) | 2.79× (1.4–4.5) | 1.97× (1.1–3.9) | 2.29× (1.0–5.9) | 1.90× (1.0–4.1) | 2.96× (2.0–4.9) | 1.61× (1.1–2.3) | 1.42× (1.0–2.7) |
| ComplexF32 fft | 1.68× (0.6–2.8) | 4.28× (2.0–8.7) | 2.37× (1.4–6.5) | 2.75× (1.3–6.7) | 2.48× (1.2–4.9) | 4.31× (3.1–5.9) | 2.43× (1.8–3.6) | 1.66× (1.2–4.4) |
| Float64 rfft | 1.74× (0.8–3.5) | 3.49× (1.6–6.9) | 2.26× (1.0–4.5) | 2.43× (1.0–7.3) | 2.33× (1.1–7.7) | 3.55× (1.8–8.0) | 1.69× (1.3–3.0) | 1.51× (1.0–3.3) |
| Float32 rfft | 2.30× (1.6–3.3) | 4.52× (2.1–8.2) | 2.13× (1.1–4.2) | 2.85× (1.4–6.7) | 2.92× (1.6–6.8) | 4.33× (2.8–7.6) | 2.40× (2.1–2.8) | 1.64× (1.0–3.1) |

*Measurement notes.* The aarch64 numbers were taken on a shared host and
the x86-64 ones on a laptop part under WSL2; on both, deltas below ~1.3×
between two suite runs an hour apart are not evidence of anything — four
apparent regressions in this work turned out to be run-to-run drift, and
the one real one (§4.6's strided real pencils, found by the x86-64 run and
fixed in the plan) was confirmed only by measuring all candidate branches
back-to-back in one session. Every per-PR claim above rests on such a
same-session comparison or on a matched-cases geometric mean over hundreds
of sizes, never on a single cross-run cell.

Against the ceiling in §6: every class landed inside the "after §4.1–4.7
(no SIMD)" column or better, and the composite/prime classes went from
5–10× to 2–3× as predicted. Two qualifications matter for anyone quoting
these numbers. First, the ratios are against FFTW's default `ESTIMATE`
plans; FFTA now matches or beats those for mid-size powers of two, but
FFTW `MEASURE` plans are 1.6–2.4× faster than `ESTIMATE` there, so against
tuned FFTW the state is 1.5–1.8× behind at 2^14–2^19 and 2.7–3.2× behind
at n ≥ 2^20 (§5.4). Second, the structural items are untouched: `ComplexF32`
is still no faster than `ComplexF64` in FFTA (§5.3), the non-radix-4 kernels
are where x86-64 loses most (§5.3a), Rader is absent (§5.5), and threading
within a single 1D transform does not exist (§5.6). Those are the remaining
work list: radix-5/7 butterflies (first on x86-64), a vectorised codelet
generator, a cache-blocked large-n path, Rader.
