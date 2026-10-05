# Final cross-basis benchmark

All representations encode the same deterministic polynomial.
Lagrange values are generated outside the timed loop using a radix-2 NTT.

| N | Method | Time (ms) | Relative to kernel |
|---:|---|---:|---:|
| 1024 | kernel_optimized | 0.004343 | 1.000x |
| 1024 | monomial_parallel_4096 | 0.008613 | 1.983x |
| 1024 | monomial_horner | 0.008708 | 2.005x |
| 1024 | lagrange_parallel | 0.109710 | 25.260x |
| 1024 | lagrange_serial | 0.021612 | 4.976x |
| 4096 | kernel_optimized | 0.016975 | 1.000x |
| 4096 | monomial_parallel_4096 | 0.034599 | 2.038x |
| 4096 | monomial_horner | 0.035491 | 2.091x |
| 4096 | lagrange_parallel | 0.250220 | 14.741x |
| 4096 | lagrange_serial | 0.084269 | 4.964x |
| 16384 | kernel_optimized | 0.059153 | 1.000x |
| 16384 | monomial_parallel_4096 | 0.069174 | 1.169x |
| 16384 | monomial_horner | 0.146720 | 2.480x |
| 16384 | lagrange_parallel | 0.439930 | 7.437x |
| 16384 | lagrange_serial | 0.334240 | 5.650x |
| 65536 | kernel_optimized | 0.091116 | 1.000x |
| 65536 | monomial_parallel_4096 | 0.139480 | 1.531x |
| 65536 | monomial_horner | 0.558910 | 6.134x |
| 65536 | lagrange_parallel | 1.250300 | 13.722x |
| 65536 | lagrange_serial | 1.339100 | 14.697x |
| 262144 | kernel_optimized | 0.263230 | 1.000x |
| 262144 | monomial_parallel_4096 | 0.380910 | 1.447x |
| 262144 | monomial_horner | 2.244000 | 8.525x |
| 262144 | lagrange_parallel | 4.561200 | 17.328x |
| 262144 | lagrange_serial | 5.367900 | 20.392x |
| 1048576 | kernel_optimized | 1.056000 | 1.000x |
| 1048576 | monomial_parallel_4096 | 1.411600 | 1.337x |
| 1048576 | monomial_horner | 9.040200 | 8.561x |
| 1048576 | lagrange_parallel | 17.199000 | 16.287x |
| 1048576 | lagrange_serial | 21.389000 | 20.255x |
