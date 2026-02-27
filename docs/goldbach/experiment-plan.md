# Goldbach Fiber-Sum Experiments — Plan (S1.2 / Option 2)

## Overview

This document describes the numerical experiments that probe the
**Goldbach fiber-sum** approach (S1.2, option 2).  The central objects
are the exponential sum

```
Phi(α) = Σ_{p ≤ N}  log(p) · e(α·p),   e(x) = exp(2πi·x)
```

and its squared convolution

```
F(t)  = Phi(α + t/W)²,   t ∈ {0, 1, …, W−1}
```

whose discrete Fourier transform

```
F̂(k) = Σ_{t=0}^{W−1}  F(t) · e(−kt/W)
```

controls the number of representations of `n ≈ 2N` as a sum of two
primes via

```
R(n) ≈ (1/W) Σ_{k}  F̂(k) · e(kn/W).
```

---

## Parameters

| Symbol | Meaning | Suggested range |
|--------|---------|-----------------|
| `N`    | Prime-sum upper bound | 10³ – 10⁶ |
| `W`    | Modulus (primorial or product of small primes) | 30, 210, 2310, … |
| `a/q`  | Base major-arc centre | 0/1, 1/2, 1/3, … |
| `Q`    | Major-arc half-width threshold | 1/N … (log N)²/N |
| `S`    | Number of minor-arc α samples | 100 – 10 000 |

---

## Metrics and Theoretical Lemma Mappings

### M1 — Mean and variance of |Phi(α + t/W)|

**What we measure:**  
For a fixed minor-arc base `α`, collect `|Phi(α + t/W)|` for
`t = 0, …, W−1`.  Compute mean `μ`, standard deviation `σ`, and the
flatness ratio `σ/μ`.

**Theoretical lemma:**  
The Hardy–Littlewood circle method (Lemma 2.1 in the standard
exposition) predicts that for α in the *minor arc*,
`|Phi(α)| = O(N^{1/2+ε})`.  A flat distribution with `σ/μ` small
indicates good cancellation consistent with this bound.  Spikes near
`σ/μ = 1` signal proximity to an undiscovered major arc.

---

### M2 — Aliasing count

**What we measure:**  
Given base `α` and threshold `Q`, count how many of the shifted points
`α + t/W` (mod 1) fall within `Q` of a rational `a/q` with
`q ≤ (log N)^B` (major-arc condition).

**Theoretical lemma:**  
Lemma 3.4 (minor-arc separation): if `α` is δ-far from every major arc
centre then at most `O(W·Q·(log N)^B)` shifts land inside major arcs.
Measuring this count validates that our minor-arc samples are
well-separated.

---

### M3 — Spectrum of F̂(k) for k ≈ 2N mod W

**What we measure:**  
Compute the DFT `F̂(k)` and report the magnitude at `k* = 2N mod W`.
Also report `max_{k≠0} |F̂(k)| / |F̂(k*)|` (off-diagonal leakage).

**Theoretical lemma:**  
By the fiber-sum identity (Lemma S1.2-A), the dominant contribution to
R(2N) comes from the major-arc integral, which corresponds to a
near-delta spike at `k*`.  Large off-diagonal ratios indicate that minor
arcs contribute non-negligibly, which would conflict with the
Goldbach conjecture's prediction.

---

### M4 — Character correlations on ℤ/Wℤ

**What we measure:**  
For each Dirichlet character χ (mod W), compute

```
C(χ) = (1/W) Σ_{t=0}^{W−1}  F(t) · conj(χ(t))
```

and report `|C(χ)|` normalised by `|F̂(0)|`.

**Theoretical lemma:**  
Lemma 4.2 (character orthogonality): if `F(t)` were a pure character
multiple, `C(χ)` would equal 1 for exactly one χ and 0 for all others.
Broadly spread `|C(χ)|` values indicate that F(t) is approximately
equidistributed over characters, consistent with minor-arc
cancellation.

---

### M5 — Convergence with N

**What we measure:**  
Run M1–M4 for a geometric sequence `N = N0, 2N0, 4N0, …` and fit power
laws to `μ`, `σ`, and `|F̂(k*)|` as functions of N.

**Theoretical lemma:**  
If `Phi(α) ~ N^{1/2}` on average (square-root cancellation, consistent
with GRH), then `|F̂(k*)| ~ N`.  Deviation from this scaling would
suggest either major-arc contamination or a breakdown of the
cancellation hypothesis.

---

## Workflow

```
generate primes ≤ N          (sieve_primes)
          │
          ▼
evaluate Phi(α + t/W)        (eval_phi, vectorised FFT-based)
  for t = 0, …, W−1
          │
          ▼
compute F(t) = Phi(...)²     (compute_F)
          │
      ┌───┴────────────┐
      ▼                ▼
DFT of F(t)      character correlations
(numpy.fft)      (character_correlations)
      │
      ▼
report M1–M5 statistics      (run_experiment, print_report)
```

---

## CLI Usage

See [`tools/goldbach_fiber_experiments.py`](../../tools/goldbach_fiber_experiments.py)
for the implementation.

```bash
# Basic experiment: N=10000, W=30, base α=0 (major arc at 0/1)
python tools/goldbach_fiber_experiments.py --N 10000 --W 30

# Minor-arc sample: α near 3/7, 200 random samples, W=210
python tools/goldbach_fiber_experiments.py --N 50000 --W 210 --a 3 --q 7 --samples 200

# Self-check / doctest mode
python tools/goldbach_fiber_experiments.py --selfcheck

# Full parameter sweep (slow)
python tools/goldbach_fiber_experiments.py --N 100000 --W 2310 --a 1 --q 11 --samples 500 --Q 0.01
```

---

## Interpretation Guide

| Observation | Likely meaning |
|-------------|---------------|
| `σ/μ < 0.3` on minor arcs | Good cancellation; consistent with GRH |
| Aliasing count > `W/q` | α is too close to a major arc; resample |
| Off-diagonal leakage < 1 % | Major arc dominates R(2N) as expected |
| `|C(χ)| ≈ 1/√φ(W)` for all χ | F(t) is equidistributed over characters |
| `|F̂(k*)| ~ N` scaling | Square-root cancellation confirmed numerically |
