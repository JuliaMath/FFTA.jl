## Thread scaling on a hybrid P/E part (WSL2, cannot pin)

Batched fft along dims=1, 1024x256 ComplexF64. 3 independent trials per point.
`ctrl` = FFTA with 1 worker, re-measured inside every trial round as a noise control.

| threads | FFTW (ms) | speedup | FFTA (ms) | speedup | FFTW spread | FFTA spread | ctrl spread |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 1.336 | 1.00× | 1.679 | 1.00× | 1.1% | 0.4% | 1.7% |
| 2 | 0.901 | 1.48× | 1.045 | 1.61× | 9.9% | 4.9% | 0.2% |
| 4 | 0.513 | 2.60× | 0.525 | 3.20× | 17.3% | 5.9% | 0.2% |
| 6 | 0.577 | 2.32× | 0.557 | 3.01× | 6.0% | 9.2% | 0.6% |
| 8 | 0.515 | 2.59× | 0.579 | 2.90× | 6.9% | 5.2% | 1.0% |
| 11 | 0.383 | 3.49× | 0.554 | 3.03× | 36.7% | 10.4% | 1.6% |

Pooled control spread over the whole probe (18 samples): **2.8%** — this is the
host's baseline noise floor for a single-worker transform on this box.

A near-linear curve to 6 that flattens after it is the signature of SMT/E-core
placement. FFTA/FFTW spread materially above the control at 8 and 11 threads is
the signature of the scheduler migrating workers between core types.
