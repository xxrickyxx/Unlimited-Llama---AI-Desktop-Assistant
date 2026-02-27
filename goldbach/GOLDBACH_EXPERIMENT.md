# Goldbach Experiment: Major and Minor Arcs

## Overview

This experiment studies the Hardy–Littlewood circle method applied to Goldbach's conjecture.
For a target even integer **N**, we analyse the exponential sum

```
Φ(α) = Σ_{p ≤ N} log(p) · e(α p)
```

where `e(x) = exp(2πi x)` and the sum is over primes `p ≤ N`.

The integral `∫₀¹ Φ(α)² e(−α N) dα` counts (with log-weight) the number of ways to write
`N = p + q` with primes `p, q ≤ N`.

---

## Major / Minor Arcs — Case 1

### Parameter

```
P = (log N)^A
```

where **A > 0** is a user-chosen exponent (default `A = 2`).

### Definition

The **major arc** `𝔐(a, q)` around the rational point `a/q` is the set of real `α` satisfying

```
|α − a/q| ≤ P / (N · q)
```

The **major arcs** `𝔐` are the union of all such arcs over integers `a`, `q` with

```
1 ≤ q ≤ P,   0 ≤ a < q,   gcd(a, q) = 1
```

The **minor arcs** `𝔪 = [0, 1) \ 𝔐` are the complement.

### Width of each arc

Each major arc `𝔐(a, q)` has width `2P/(N·q)`.  The total measure of `𝔐` is at most
`O(P² / N)`, which is small for `P = (log N)^A`.

---

## W-trick and Fiber DFT

When using a smoothing modulus **W** (e.g., `W = primorial(B)` for a small bound `B`),
the circle `[0, 1)` is partitioned into **W** cosets.  The **fiber DFT** of the
residue-`b` contribution uses the Fourier coefficient at index

```
k ≡ N (mod W),   k ∈ {0, 1, …, W−1}
```

because we are looking for representations `N = p + q` and the shift by `N` selects the
correct fiber.

### Aliasing count

For a fixed `α` and modulus `W`, the **aliasing count** is the number of shifts
`t ∈ {0, 1, …, W−1}` such that `α + t/W` lies inside the major arcs `𝔐`.

---

## Files

| File | Description |
|------|-------------|
| `goldbach_arcs.py` | Helper functions and CLI for major/minor arc classification |

---

## CLI Quick Reference

```
python goldbach_arcs.py --N 1000000 --A 2.0 --alpha 0.00031 --W 30
```

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--N` | int | *required* | Target even integer N |
| `--A` | float | 2.0 | Exponent in P = (log N)^A |
| `--alpha` | float | *required* | Value of α to classify |
| `--W` | int | 1 | Smoothing modulus W |

Example output:

```
N = 1000000, A = 2.0, P = 190.87 (= (log N)^A)
alpha = 0.00010
Classification: MAJOR  (closest rational: 0/1, |alpha - a/q| = 1.00e-04 <= P/(Nq) = 1.91e-04)
Fiber DFT coefficient index k = N mod W = 10
Aliasing count (t in [0,W) s.t. alpha+t/W in major arcs): ...
Phi(alpha) = ...
```
