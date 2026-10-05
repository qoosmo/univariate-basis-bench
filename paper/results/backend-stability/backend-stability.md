# Backend stability at N = 2^20

Independent Criterion invocations: **5**

## Per-run central estimates

| Run | Kernel special (us) | Kernel Montgomery (us) | Monomial special (us) | Monomial Montgomery (us) | Kernel M/S | Monomial M/S |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 883.610 | 1172.400 | 1453.500 | 1531.000 | 1.3268x | 1.0533x |
| 2 | 989.030 | 1227.700 | 1559.400 | 1654.200 | 1.2413x | 1.0608x |
| 3 | 1080.900 | 1159.900 | 1537.000 | 1901.100 | 1.0731x | 1.2369x |
| 4 | 1036.000 | 1323.700 | 2208.300 | 1666.000 | 1.2777x | 0.7544x |
| 5 | 952.150 | 1277.600 | 1607.800 | 1702.800 | 1.3418x | 1.0591x |

## Across-run summary

| Benchmark | Median (us) | Min (us) | Max (us) | CV |
|---|---:|---:|---:|---:|
| kernel_goldilocks_special | 989.030 | 883.610 | 1080.900 | 7.69% |
| kernel_montgomery | 1227.700 | 1159.900 | 1323.700 | 5.63% |
| monomial_goldilocks_special | 1559.400 | 1453.500 | 2208.300 | 18.19% |
| monomial_montgomery | 1666.000 | 1531.000 | 1901.100 | 7.93% |

## Within-run ratios

- Median Montgomery / specialized kernel ratio: **1.2777x**
- Median Montgomery / specialized monomial ratio: **1.0591x**

Interpretation rule for the paper:
- Use within-run ratios for backend comparisons.
- Use the median of independent runs as the reproducibility summary.
- Do not mix absolute timings from this experiment with the earlier cross-basis sweep.
- If the monomial M/S ratios straddle 1 materially across runs, report no clear backend separation for monomial evaluation.
- If the kernel M/S ratio remains consistently above 1, report the kernel-specific specialization effect.
