-- codelet vs A's radix-4 recursion (ComplexF64, contiguous), and FFTW
N=  8 compile=0.30s  codelet=  11.8 ns  A-recursion=  37.5 ns  fftw=  20.3 ns  codelet/fftw=0.58  speedup vs A=3.2x
N= 16 compile=0.29s  codelet=  36.1 ns  A-recursion=  56.4 ns  fftw=  28.1 ns  codelet/fftw=1.28  speedup vs A=1.6x
N= 32 compile=0.17s  codelet=  99.4 ns  A-recursion= 201.0 ns  fftw=  58.5 ns  codelet/fftw=1.70  speedup vs A=2.0x
N= 64 compile=0.42s  codelet= 266.1 ns  A-recursion= 307.9 ns  fftw= 207.5 ns  codelet/fftw=1.28  speedup vs A=1.2x

-- full transforms with codelet base cases (BASE=64) vs A vs FFTW
N=     256  base64= 1414.5 ns  base16= 1271.5 ns  A= 1434.2 ns  fftw= 1125.8 ns   base64/fftw=1.26  A/fftw=1.27
N=    4096  base64=36323.0 ns  base16=34057.0 ns  A=35785.0 ns  fftw=32423.0 ns   base64/fftw=1.12  A/fftw=1.10
N=   65536  base64=818460.0 ns  base16=798171.0 ns  A=815769.0 ns  fftw=823978.0 ns   base64/fftw=0.99  A/fftw=0.99
N= 1048576  base64=53424055.0 ns  base16=66597136.0 ns  A=57165204.0 ns  fftw=33230754.0 ns   base64/fftw=1.61  A/fftw=1.72

64-pt codelet with input stride 64: 276.5 ns  (contiguous: 266.3 ns)
