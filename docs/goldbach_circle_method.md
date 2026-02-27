# Goldbach / Circle-Method Cutoff Conventions

This document clarifies the **N vs 2N cutoff** used in Hardy–Littlewood
circle-method style sums that may appear in experiment scripts within this
repository.

---

## Background

The **Hardy–Littlewood circle method** represents the number of ways `2N` can
be written as a sum of two primes via the integral

```
r(2N) = ∫₀¹ S(α)² e(-2N α) dα
```

where the **exponential sum** (major-arc approximant) is

```
S(α) = Σ_{p ≤ cutoff_max} e(p α)
```

and `e(x) = exp(2πi x)`.

---

## N vs 2N: What changes?

| Setting | Effect on S(α) | Effect on r(2N) |
|---|---|---|
| `cutoff_max = 2N` | Sums over all primes up to `2N` | Captures every pair `(p, q)` with `p + q = 2N`, `p, q ≤ 2N`. This is the **correct support** for the convolution. |
| `cutoff_max = N` | Sums over primes up to `N` only | Only pairs with both `p ≤ N` and `q ≤ N` are counted. Pairs such as `(3, 2N-3)` where `2N-3 > N` are **silently dropped**, under-counting `r(2N)`. |

### Why `2N` is the correct default

Because `p + q = 2N` and both `p, q ≥ 2`, the largest prime that can appear
in any valid pair is `2N - 2`.  Setting `cutoff_max = 2N` guarantees that
no valid pair is excluded from the sum.

Setting `cutoff_max = N` is sometimes used as a **symmetry shortcut** in
analytic proofs (where both primes are assumed ≤ N by symmetry and the count
is then doubled), but this requires an explicit factor-of-2 correction and
can silently produce wrong counts when that correction is omitted.

---

## Normalization

When `cutoff_max = 2N` (the default):

```
r(2N) = #{(p, q) : p prime, q prime, p + q = 2N, p ≤ 2N, q ≤ 2N}
```

When `cutoff_max = N` (symmetric shortcut):

```
r_sym(2N) ≈ ½ · r(2N)      (approximate, valid for large N)
```

so a corrected estimate is `r(2N) ≈ 2 · r_sym(2N)`.  The approximation
becomes exact when `2N` has no representation with `p = N` (i.e. `N` is
not prime, or `N` is prime but the `(N, N)` pair is counted only once).

**Recommended practice:** always use `cutoff_max = 2N` in experiment scripts
unless you are deliberately applying the symmetric shortcut *and* applying
the factor-of-2 correction explicitly.

---

## Recommended Defaults for Experiment Scripts

```python
# Correct default – covers full convolution support
cutoff_max = 2 * N

# Symmetric shortcut – remember to multiply result by 2
cutoff_max = N   # WARNING: under-counts r(2N) without ×2 correction
```

See [`../experiments/goldbach_experiment.py`](../experiments/goldbach_experiment.py)
for a reference implementation that accepts `--cutoff-max` and defaults to `2N`.

---

## Quick Reference

- **Target sum:** `2N`
- **Correct cutoff:** `cutoff_max = 2N`
- **Symmetric shortcut:** `cutoff_max = N` + multiply count by 2
- **Formula:** `r(2N) = ∫₀¹ S(α)² e(-2Nα) dα`  with  `S(α) = Σ_{p ≤ cutoff_max} e(pα)`
