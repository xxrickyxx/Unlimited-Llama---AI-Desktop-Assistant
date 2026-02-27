# Goldbach Conjecture – Circle Method Research Plan

## Overview

**Goldbach's Conjecture** states that every even integer `N ≥ 4` is the sum of two primes.
This document outlines a research roadmap towards a numerical/analytic verification
strategy using the Hardy–Littlewood **circle method**, combined with the
fiber/W-trick, dispersion estimates, and Kuznetsov/Kloosterman summation.

---

## 1. The Circle Method Setup

### 1.1 Exponential Sum Decomposition

Define the von Mangoldt exponential sum (major/minor arc):

```
S(alpha) = sum_{n <= N} Lambda(n) * e(n * alpha)
```

where `e(x) = exp(2*pi*i*x)` and `Lambda` is the von Mangoldt function.

The Goldbach counting function is:

```
R(N) = integral_{0}^{1} S(alpha)^2 * e(-N * alpha) d alpha
```

We split the unit circle into **major arcs** M and **minor arcs** m:

```
R(N) = R_major(N) + R_minor(N)
```

Goal: show `R(N) > 0` for all even `N >= 4`.

---

## 2. The W-Trick / Fiber Reduction

### 2.1 Motivation

Direct treatment of `S(alpha)` near rational approximations `a/q` is hard.
The **W-trick** (Green–Tao style) introduces a modular filter:

```
W = prod_{p <= w(N)} p          (primorial up to w(N) = log log N)
```

and replaces primes by **W-smooth residue classes**:

```
Lambda_W(n) = (phi(W)/W) * Lambda(n) * 1[gcd(n, W) = 1]
```

### 2.2 Fiber Decomposition

Fix a residue `b mod W` with `gcd(b, W) = 1`.  Define the fiber sum:

```
S_b(alpha) = sum_{n <= N, n ≡ b (W)} Lambda(n) * e(n * alpha)
```

The full sum decomposes as:

```
S(alpha) = sum_{b : gcd(b,W)=1} S_b(alpha)
```

### 2.3 Shifted Alpha

For the fiber analysis it is natural to write:

```
alpha = a/q + beta,   and separately   alpha = t/W + gamma
```

This **aliasing** picture is the core of the fiber experiment (see §5).

---

## 3. Major Arc Analysis

### 3.1 Standard Major Arc Estimate

On the major arc around `a/q` (with `q <= Q = (log N)^A`):

```
S(a/q + beta) ≈ mu(q)/phi(q) * V(beta)
```

where `V(beta) = sum_{n<=N} e(n*beta)` is the standard Fejér-type kernel.

### 3.2 Singular Series

The major arc contribution gives the **singular series**:

```
S(N) = prod_p [ 1 - (-1)/((p-1)^2) ] * prod_{p | N, p>2} (p-1)/(p-2)
```

**Lemma 3.1 (Lower bound):** For all even `N`, `S(N) >= c > 0`.

_Proof sketch:_ Standard Euler product bound; see Hardy–Littlewood (1923)._

---

## 4. Minor Arc Bounds via Dispersion and Kloosterman

### 4.1 Dispersion Method

To bound `integral_m |S(alpha)|^2 d alpha` we use the **dispersion method**:

```
integral |S|^2 d alpha = sum_{n} r(n) * Lambda(n)
```

where `r(n)` counts representations `n = a - b` with `a,b <= N`.

Vaughan's identity decomposes `Lambda` into:

```
Lambda = Lambda_1 - Lambda_2 * d - Lambda_3 * Lambda_4
```

allowing bilinear forms amenable to Cauchy-Schwarz.

### 4.2 Kuznetsov / Kloosterman Summation

For the bilinear sum in the type-II terms:

```
Sigma_II = sum_{m,n} a_m b_n e(mn/q)
```

we apply the **Kuznetsov trace formula** to convert the sum to spectral data
(Maass forms, holomorphic forms, Eisenstein series), then use the
**Weil bound** on Kloosterman sums:

```
|Kl(a, b; q)| <= d(q) * sqrt(q)
```

**Lemma 4.1 (Vaughan type-II bound):** Under GRH (or unconditionally for `q <= N^{1/2}`):

```
|Sigma_II| << N^{1/2 + epsilon} * (sum |a_m|^2)^{1/2} * (sum |b_n|^2)^{1/2}
```

---

## 5. The Fiber Aliasing Experiment

### 5.1 Discrete Grid Picture

Fix `N`, `W`, and evaluate `Phi(alpha)` on the fine grid:

```
alpha_t = t / W,   t = 0, 1, ..., W-1
```

For each residue class `b mod W` the fiber sum `S_b(alpha_t)` exhibits
**aliasing**: contributions at `alpha_t` and `alpha_{t+W}` are identified.

### 5.2 Aliasing Count

Define the aliasing multiplicity at frequency `t`:

```
A(t) = #{(b, k) : b*k ≡ t (mod W), gcd(b,W) = 1, k >= 1, b*k <= N}
```

Large `A(t)` signals constructive interference and dominant major arc behavior.

### 5.3 Numerical Goals

For moderate `N` (say `N = 200`, `1000`, `10000`):

1. Compute `|S(t/W)|` for all `t` and identify peaks.
2. Compare peak positions to the major arc rational approximations `a/q`.
3. Estimate `integral_m |S|^2 d alpha` numerically and compare to the
   theoretical upper bound from §4.

See `tools/goldbach_fiber_experiments.py` for the implementation.

---

## 6. Roadmap / Next Steps

| Step | Status | Description |
|------|--------|-------------|
| S1 | ✅ done | Set up W-trick decomposition and fiber sums |
| S2 | ✅ done | Major arc singular series computation |
| S3 | 🔄 in progress | Numerical aliasing grid experiment |
| S4 | ⬜ todo | Vaughan type-II bilinear form bound (numerical) |
| S5 | ⬜ todo | Kuznetsov spectral expansion (symbolic/numerical) |
| S6 | ⬜ todo | Full minor arc integral upper bound vs `R_major` lower bound |
| S7 | ⬜ todo | Verification loop: check `R(N) > 0` for `N <= 10^6` |

---

## 7. Key References

- Hardy, G. H. and Littlewood, J. E. (1923). *Some problems of Partitio Numerorum III*.
- Vinogradov, I. M. (1937). *Representation of an odd number as a sum of three primes*.
- Vaughan, R. C. (1997). *The Hardy–Littlewood Method*, 2nd ed.
- Green, B. and Tao, T. (2008). *The primes contain arbitrarily long arithmetic progressions*.
- Kuznetsov, N. V. (1980). *The Petersson conjecture for cusp forms of weight zero*.
- Iwaniec, H. and Kowalski, E. (2004). *Analytic Number Theory*.
