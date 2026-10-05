# Goldilocks special reduction vs Montgomery

| Benchmark | Central estimate |
|---|---:|
| mul_backend/goldilocks_special_mul | 1.910700 ns |
| mul_backend/montgomery_mul | 1.916200 ns |
| mul_backend/goldilocks_kernel_fold | 3.029200 ns |
| mul_backend/montgomery_kernel_fold | 3.795200 ns |
| backend_full_eval_n1048576/kernel_goldilocks_special/1048576 | 918.660000 us |
| backend_full_eval_n1048576/kernel_montgomery/1048576 | 1.191100 ms |
| backend_full_eval_n1048576/monomial_goldilocks_special/1048576 | 1.643100 ms |
| backend_full_eval_n1048576/monomial_montgomery/1048576 | 1.563600 ms |

## Ratios

- Montgomery / Goldilocks multiplication: 1.0029x
- Montgomery / Goldilocks kernel fold: 1.2529x
- Montgomery / Goldilocks full kernel: 1.2966x
- Montgomery / Goldilocks full monomial: 0.9516x
