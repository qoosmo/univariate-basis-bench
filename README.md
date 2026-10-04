# Univariate Basis Bench

Rust benchmarks for **univariate polynomial representations relevant to ZK provers**.

The project intentionally keeps the scope narrow:

1. monomial basis,
2. Lagrange/evaluation basis on a power-of-two roots-of-unity domain,
3. Boolean-kernel basis

\[
K_y(X)=\prod_{i=0}^{m-1}(X^{2^i}+y_i),\qquad y_i\in\{0,1\}.
\]

No multilinear-polynomial benchmarks are included.

## Benchmark 01 — single-point evaluation

Given one degree-`< N` univariate polynomial and a point `r`, measure the wall-clock cost of computing `P(r)` from each representation.

Implementations:

- `monomial_horner`: ordinary Horner evaluation.
- `monomial_parallel_blocks`: fair intra-polynomial parallel baseline obtained by splitting the coefficient vector into independent blocks.
- `lagrange_barycentric`: roots-of-unity Lagrange evaluation using one batch inversion.
- `lagrange_barycentric_parallel`: parallel denominator construction and parallel final accumulation.
- `kernel_serial`: tree folding in the Boolean-kernel basis.
- `kernel_parallel`: level-by-level Rayon parallel tree folding.

The kernel fold at level `i` is

\[
(a,b)\mapsto t(a+b)+b,\qquad t=r^{2^i}.
\]

All independent pairs at one level can execute concurrently.

## Why Goldilocks first?

Benchmark 01 needs a genuine power-of-two roots-of-unity Lagrange domain, so the initial common field is Goldilocks

\[
p=2^{64}-2^{32}+1.
\]

This avoids giving Lagrange an artificial non-FFT domain. The basis algorithms are field-generic conceptually. A production-field backend (including M31 where the target protocol requires it) can be added after the operation-level methodology is fixed.

## Run

```bash
cargo run --release > results.csv
cargo bench --bench eval
```

For thread scaling:

```bash
for t in 1 2 4 8 16; do
  RAYON_NUM_THREADS=$t cargo run --release > results-${t}t.csv
done
```

## Methodology rules

- Same field implementation for every basis in a given comparison.
- Same `N` and same benchmark machine.
- Report wall-clock time and ns/element.
- Benchmark serial and parallel algorithms separately.
- Do not intentionally serialize competing bases when a legitimate parallel algorithm exists.
- Keep setup/precomputation out of a timed loop only when that setup is reusable in the corresponding prover workflow; document every such exclusion.
- Validate algorithms against a common polynomial before treating timings as publishable.

## Roadmap

We add operations one at a time:

- [x] 01 single-point evaluation
- [x] correctness tests for monomial ↔ kernel evaluation and Lagrange reconstruction
- [ ] 02 full-domain evaluation / encoding
- [ ] 03 low-degree extension
- [ ] 04 basis conversion
- [ ] 05 polynomial folding/reduction used by FRI-style provers
- [ ] thread-scaling plots and memory-bandwidth counters

