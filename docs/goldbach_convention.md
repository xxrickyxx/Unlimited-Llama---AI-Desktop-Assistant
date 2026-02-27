# Goldbach Circle-Method Convention: N as the Target Even Integer

This document describes the cutoff and DFT conventions used in this repository's
Goldbach-experiment tooling.  Throughout, **N denotes the even integer we wish to
represent as a sum of two primes** (the "target").

---

## 1. Exponential Sum and Representation Count

Define the von-Mangoldt-weighted exponential sum truncated at **N**:

```
S(α) = ∑_{n=1}^{N} Λ(n) · e(αn)       where e(x) = exp(2πi x)
```

The corresponding weighted representation count is

```
r_Λ(N) = ∫₀¹ S(α)² · e(−Nα) dα
```

which counts all ordered pairs `(m, n)` with `m + n = N`, weighted by `Λ(m)·Λ(n)`.

Goldbach's conjecture (in this weighted form) asserts that `r_Λ(N) > 0` for every
even integer `N ≥ 4`.

---

## 2. Discrete Approximation (DFT)

For numerical experiments the integral is approximated by a length-`N` DFT.
Choosing the DFT grid size `W = N` (or any `W ≥ N`), the relevant coefficient
satisfies

```
k ≡ N  (mod W)
```

so that the phase `e(−Nα_k)` aligns with the target.  The experiment script
defaults to `W = N` so that `k = N mod N = 0` selects the zero-frequency bin,
which carries `r_Λ(N)`.

In code (see `scripts/goldbach_experiment.py`):

```python
W = N                       # DFT grid size equals the target even integer
cutoff_max = N              # sum runs from n=1 to N (inclusive)
k = N % W                   # DFT coefficient index  →  k = 0 for W = N
```

---

## 3. Variable-naming Convention

| Symbol | Meaning |
|--------|---------|
| `N` / `target_even` | The even integer being represented (e.g. 10, 28, 100) |
| `cutoff_max` | Upper limit of the von-Mangoldt sum; equals `N` in this convention |
| `W` | DFT grid length; defaults to `N` |
| `k` | DFT frequency bin for the target; `k = N % W` |
| `Λ(n)` | Von Mangoldt function (`log p` if `n = p^r`, else `0`) |

The variable `target_even` is used in the script to make the role of `N`
self-documenting:

```python
target_even = 28            # the even N we want to write as p + q
```

---

## 4. Alternative Convention: target = 2N

Some references use a shifted convention where the argument passed to the sum is
**half** the target, i.e. they write the target as `2M` and build the sum up to `M`:

```
S_alt(α) = ∑_{n=1}^{M} Λ(n) · e(αn)

r_Λ(2M) = ∫₀¹ S_alt(α)² · e(−2Mα) dα
```

### How formulas change

| Quantity | This repo (target = N, even) | Alternative (target = 2M) |
|----------|------------------------------|---------------------------|
| Sum upper limit | `cutoff_max = N` | `cutoff_max = M = N/2` |
| Phase exponent | `e(−Nα)` | `e(−2Mα)` |
| DFT grid size | `W = N` | `W = 2M` (same integer!) |
| DFT bin | `k = N % W = 0` | `k = 2M % W = 0` |
| Printed target | `N` (e.g. 28) | `2M` (e.g. 28, M=14) |

To switch the experiment script to the alternative convention, set
`use_half_convention = True` (see `scripts/goldbach_experiment.py`).  The script
will then use `cutoff_max = N // 2` and label output with `M = N // 2`.

> **Note:** Both conventions produce the same `r_Λ` value for the same even integer
> because `S_alt(α)² e(−2Mα)` evaluated on the full unit interval is identical to
> the full-range integral—the sum just runs to a smaller cutoff.  However, the
> **truncation error** differs: using `cutoff_max = N` retains more Λ-weight and
> gives a more accurate numerical estimate for fixed grid size.

---

## 5. Quick Reference

```
# Primary convention (this repo):
N          = 28             # target even integer
cutoff_max = N              # = 28
k          = N % N          # = 0
S(α)       = Σ_{n≤28} Λ(n) e(αn)
r_Λ(28)    = DFT[S·S][k=0]

# Alternative convention:
M          = N // 2         # = 14
cutoff_max = M              # = 14
k          = (2*M) % (2*M)  # = 0
S_alt(α)   = Σ_{n≤14} Λ(n) e(αn)
r_Λ(28)    = DFT[S_alt·S_alt][k=0]
```
