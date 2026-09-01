## DFT vs Bluestein crossover — x86-64 (ComplexF64, planned exec, min of 200)

| n | DFT leaf (µs) | Bluestein (µs) | Bluestein/DFT | cheaper |
|--:|--:|--:|--:|:--|
| 13 | 0.34 | 0.85 | 2.50× | DFT |
| 17 | 0.57 | 0.83 | 1.45× | DFT |
| 19 | 0.71 | 0.82 | 1.14× | DFT |
| 23 | 1.19 | 0.91 | 0.76× | **Bluestein** |
| 29 | 1.87 | 1.70 | 0.91× | **Bluestein** |
| 31 | 2.19 | 0.90 | 0.41× | **Bluestein** |
| 37 | 3.15 | 2.28 | 0.72× | **Bluestein** |
| 41 | 3.88 | 2.29 | 0.59× | **Bluestein** |
| 43 | 4.24 | 2.30 | 0.54× | **Bluestein** |
| 47 | 5.13 | 2.31 | 0.45× | **Bluestein** |
| 53 | 6.53 | 2.33 | 0.36× | **Bluestein** |
| 59 | 8.11 | 2.35 | 0.29× | **Bluestein** |
| 61 | 8.70 | 2.34 | 0.27× | **Bluestein** |
| 67 | 10.54 | 3.68 | 0.35× | **Bluestein** |
| 71 | 11.82 | 3.71 | 0.31× | **Bluestein** |
| 73 | 12.54 | 3.71 | 0.30× | **Bluestein** |
| 79 | 14.71 | 3.72 | 0.25× | **Bluestein** |
| 83 | 16.32 | 3.74 | 0.23× | **Bluestein** |
| 89 | 18.85 | 3.77 | 0.20× | **Bluestein** |
| 97 | 22.36 | 3.80 | 0.17× | **Bluestein** |

**x86-64 crossover: Bluestein first wins at n = 23** (aarch64 value: 47)
