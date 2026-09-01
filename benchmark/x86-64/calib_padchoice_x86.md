## Bluestein pad choice, measured end-to-end on x86-64

`A` = `feat/plan-twiddles` @ 8411861 (pad is always `nextpow2(2N-1)`).
`E` = `feat/bluestein-padding` @ 071712b (3-smooth chooser, factor 1.9, 2048 floor).
ComplexF64, planned execution, 3 trials of min-of-300.

| n | pad on A | pad on E | A (µs) | E (µs) | E/A |
|--:|--:|--:|--:|--:|--:|
| 73 | 256 | 256 | 4.13 | 4.19 | 1.01× |
| 79 | 256 | 256 | 7.49 | 8.37 | — (same pad) |
| 4099 | 16384 | 8748 (2^2·3^7) | 365.9 | 434.9 | **1.19×** |
| 8443 | 32768 | 17496 (2^2·3^7) | 935.8 | 1109.8 | **1.19×** |
| 65537 | 262144 | 139968 | 9879.0 | 9983.0 | 1.01× |

Repeat trials (µs), showing the spread:

```
A  n=  4099 pad=  16384  365.9, 370.6, 366.8
A  n=  8443 pad=  32768  936.2, 935.8, 938.9
A  n= 65537 pad= 262144  9884.6, 9879.0, 10075.7
E  n=  4099 pad=   8748  434.8, 507.3, 506.9
E  n=  8443 pad=  17496  1109.8, 1117.7, 1111.3
E  n= 65537 pad= 139968  10104.1, 10145.7, 9983.0
```

The 2048 floor behaves correctly on both architectures (73 and 79 keep the 256 pad).
At n = 8443 the choice measured as a 1.58x win on aarch64 is a 1.19x loss here.

Factor at which each case flips from accept to reject: 4099 -> 2.00, 8443 -> 1.99,
65537 -> 1.97. The measured 3-smooth cost factor is 2.01-3.11 on x86-64 and
1.28-2.26 on aarch64, so the threshold falls between the two machines.
