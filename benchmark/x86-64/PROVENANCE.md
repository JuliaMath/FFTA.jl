# Commit provenance for the x86-64 runs

Measured on: Intel Core Ultra 7 165H (alderlake), Julia 1.12.6, FFTW 3.3.11 (provider `fftw`), -t 8.
All runs use the *same* `suite.jl`/`report.jl`/`compare.jl` copied from
pankgeorg/FFTA.jl `bench/comprehensive-suite` @ 73b190c, so before/after pairs are
methodologically identical. (The bench branch later moved to d0064b1, which adds a
`"complete": true` flag and line-flushing but changes nothing that is measured;
deliberately not adopted mid-sequence.)

| run | branch | commit measured | note |
|:--|:--|:--|:--|
| baseline | bench/comprehensive-suite | 73b190c | FFTA `src/` byte-identical to JuliaMath/FFTA.jl main @ 7aeb327 (FFTA 0.3.1); verified with `git diff origin/main...HEAD -- src/ Project.toml` (empty) |
| integration_all | integration/all | b06c41c | branch head later moved to 45696ac; `git diff b06c41c 45696ac -- src/ Project.toml` is empty (test-only changes), so these numbers are valid for 45696ac |
| A_plan_twiddles | feat/plan-twiddles | 8411861 | |
| C_pow2_codelets | feat/pow2-codelets | 5b02c9d | `--only 1d --kinds fft` |
| B_real_plan_mul | feat/real-plan-mul | c9d0a99 | `--kinds rfft` |
| D_nd_buffers | feat/nd-buffers-threads | 17ed942 | `--only nd,batched,threads` |
| E_bluestein | feat/bluestein-padding | 919628e | `--only 1d --kinds fft` |

## Final run set

| run | branch | commit | scope |
|:--|:--|:--|:--|
| baseline | bench/comprehensive-suite | 73b190c (FFTA src == main 7aeb327) | full + band |
| integration_all | integration/all | 5c472f8 (src == dd5148f, 9331e2f) | full + band |
| integration_all_prefix | integration/all | 5f17d1c (before the #132 fix) | full |
| A_plan_twiddles | feat/plan-twiddles | 8411861 | full + band |
| C_pow2_codelets | feat/pow2-codelets | 5b02c9d | 1d fft + band |
| B_real_plan_mul | feat/real-plan-mul | 600d633 (with the fix) | rfft + band |
| D_nd_buffers | feat/nd-buffers-threads | 17ed942 | nd,batched,threads |
| E_bluestein | feat/bluestein-padding | cc6a69b | 1d fft + band |

The prime band (23, 29, 31, 37, 43) was measured with the newer suite's
`--sizes` filter and merged into each run above; the diff from 73b190c touches
only the size list, that filter, a `complete` flag and a stdout flush — `timeit`
and `bench_case!` are byte-identical, which is what makes the merge valid.
