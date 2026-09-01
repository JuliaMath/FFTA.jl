```@meta
CurrentModule = FFTA
```

# Development Documentation

Here is the documentation for key parts of the development side of the package.

```@docs
CallGraphNode
CallGraph
CallGraphNode!
fft!
fft_composite!
fft_dft!
fft_pow2_radix4!
fft_pow3!
fft_bluestein!
```

## Twiddle factors

Twiddle factors are computed once, at plan time, and stored per node of the
call graph in the layout each kernel reads sequentially.

```@docs
twiddle
unit_roots
dft_twiddles
composite_twiddles
pow2_twiddles
pow3_twiddles
node_twiddles
BluesteinScratch
bluestein_pad_length
```
