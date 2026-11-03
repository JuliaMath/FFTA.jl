-- genfft-style straight-line radix-2 codelet (type-stable call): compile vs runtime vs FFTW vs FFTA
N=   4  stmts=    20  compile=   0.32 s  codelet=   2.5 ns  ffta=  18.7 ns  fftw= 19.0 ns   codelet/fftw=0.13  ffta/fftw=1.0
N=   8  stmts=    52  compile=   0.03 s  codelet=  10.4 ns  ffta= 110.1 ns  fftw= 20.1 ns   codelet/fftw=0.51  ffta/fftw=5.5
N=  16  stmts=   128  compile=   0.06 s  codelet=  31.3 ns  ffta= 137.0 ns  fftw= 28.1 ns   codelet/fftw=1.11  ffta/fftw=4.9
N=  32  stmts=   304  compile=   0.15 s  codelet=  89.5 ns  ffta= 553.3 ns  fftw= 58.0 ns   codelet/fftw=1.54  ffta/fftw=9.5
N=  64  stmts=   704  compile=   0.42 s  codelet= 234.4 ns  ffta= 770.8 ns  fftw=207.3 ns   codelet/fftw=1.13  ffta/fftw=3.7
N= 128  stmts=  1600  compile=   2.11 s  codelet= 765.7 ns  ffta=2576.6 ns  fftw=451.2 ns   codelet/fftw=1.70  ffta/fftw=5.7
N= 256  stmts=  3584  compile=   4.85 s  codelet=3719.4 ns  ffta=3761.2 ns  fftw=1027.3 ns   codelet/fftw=3.62  ffta/fftw=3.7

-- odd leaves (O(N^2) straight-line) ComplexF64
N=   3  compile=   0.03 s  codelet=   8.5 ns  ffta=  21.5 ns  fftw= 19.7 ns   codelet/fftw=0.43
N=   5  compile=   0.03 s  codelet=  29.4 ns  ffta= 249.7 ns  fftw= 20.0 ns   codelet/fftw=1.47
N=   7  compile=   0.06 s  codelet=  57.4 ns  ffta= 423.5 ns  fftw= 22.0 ns   codelet/fftw=2.61
N=  11  compile=   0.15 s  codelet= 141.7 ns  ffta= 913.2 ns  fftw= 30.8 ns   codelet/fftw=4.60
N=  13  compile=   0.18 s  codelet= 227.0 ns  ffta=1231.0 ns  fftw= 36.9 ns   codelet/fftw=6.15

-- ComplexF32 pow2 codelet
N=  16  compile=   0.22 s  codelet=  34.5 ns  fftw= 62.5 ns   codelet/fftw=0.55
N=  64  compile=   0.64 s  codelet= 291.3 ns  fftw=113.2 ns   codelet/fftw=2.57
N= 256  compile=   4.40 s  codelet=3246.8 ns  fftw=488.1 ns   codelet/fftw=6.65
