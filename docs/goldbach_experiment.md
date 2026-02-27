# Goldbach Conjecture – Fiber-Sum / Discrete Fourier-Shift Experiment

## Overview

This document describes the experimental protocol implemented in
[`goldbach_experiment.py`](../goldbach_experiment.py) for the
**fiber-sum / discrete Fourier-shift** approach to numerical experiments
related to the Goldbach conjecture.

The conjecture states that every even integer ≥ 4 is the sum of two primes.
The circle-method framework reformulates this as a question about the magnitude
of a Fourier coefficient of a product of exponential sums over primes.

---

## Mathematical Background

### Exponential sum Φ(θ)

Let *N* be a positive integer and let **P** = {p prime : p ≤ N}.  Define the
Λ-weighted (von Mangoldt) exponential sum

```
Φ(θ) = Σ_{p ∈ P}  log(p) · e(θ · p),    e(x) := exp(2πi · x)
```

An unweighted variant replaces log(p) with 1 (controlled via `--no-weighted`).

The circle-method identity for the Goldbach problem relates the number of
representations r₂(2N) = #{(p₁,p₂) ∈ P² : p₁+p₂ = 2N} to the integral

```
r₂(2N) ~ ∫₀¹ Φ(θ)² e(-2N·θ) dθ
```

### Fiber lattice and the fiber-square

Fix a rational approximation **α = a/q + β** to a major-arc centre.  Sample
Φ at *W* equally-spaced shifts around α:

```
θ_t  =  α + t/W,    t = 0, 1, …, W-1
```

Form the **fiber-square**

```
F(t)  =  Φ(θ_t)²
```

### Discrete Fourier Transform of F

Compute the DFT (via FFT):

```
F̂(k)  =  Σ_{t=0}^{W-1}  F(t) · e(-k·t/W),    k = 0,…,W-1
```

The coefficient **k₀ = 2N mod W** is the *Goldbach target frequency*: by
periodicity, F̂(k₀) accumulates the contribution of all Goldbach pairs whose
sum equals 2N modulo W.

### Spike / flatness metrics

| Metric | Meaning |
|--------|---------|
| `|F̂(k₀)|` | Signal strength at the Goldbach target frequency |
| Top-*k* `|F̂(k)|` | Whether k₀ dominates the spectrum (spike) or not (flat) |
| `|Φ|` distribution | Concentration/flatness of the exponential sum along the fiber |
| Aliasing count | Number of near-neighbours of k₀ that also spike above 2×median |

A large `|F̂(k₀)|` relative to neighbours is consistent with a major-arc
contribution dominating; flatness (no spike) suggests minor-arc dominance or
cancellation.

### Major-arc aliasing

For each k in the window [k₀ − w, k₀ + w] (excluding k₀ itself) check
whether `|F̂(k)| > threshold` where threshold = 2 × median(|F̂|).  A non-zero
aliasing count means nearby frequencies are also elevated, which can indicate
either genuine arithmetic structure or numerical artefacts from the finite
lattice.

---

## Experimental Protocol

1. **Choose N** (prime bound) and **W** (fiber width, power of 2 recommended).
2. **Choose α = a/q + β** where a/q is a Farey fraction (major-arc centre)
   and β ∈ [0, 1/q) is a small perturbation.
3. Compute Φ at all W fiber points (cached on disk of the process).
4. Form F = Φ², compute F̂ via FFT, inspect top-k magnitudes and |F̂(k₀)|.
5. Summarise the |Φ| distribution.
6. Optionally sweep many random β values (`--samples`) to obtain statistics.

---

## Running the Tool

### Prerequisites

```bash
pip install numpy
```

### Single experiment

```bash
python goldbach_experiment.py --N 100000 --W 256 --a 1 --q 6 --beta 0.0
```

### Vary the fiber width

```bash
python goldbach_experiment.py --N 500000 --W 1024 --a 1 --q 6 --beta 0.0 --topk 10
```

### Random-beta sweep (50 samples)

```bash
python goldbach_experiment.py --N 100000 --W 256 --a 1 --q 6 --samples 50
```

### Unit-weight sums (unweighted Φ)

```bash
python goldbach_experiment.py --N 100000 --W 256 --a 1 --q 6 --no-weighted
```

### Self-check / smoke tests

```bash
python goldbach_experiment.py --selfcheck
```

### Full CLI reference

```
usage: goldbach_experiment.py [-h] [--N N] [--W W] [--a A] [--q Q]
                               [--beta BETA] [--weighted | --no-weighted]
                               [--topk TOPK] [--arc_width ARC_WIDTH]
                               [--samples SAMPLES] [--seed SEED] [--selfcheck]

options:
  --N N                    Upper bound for prime sieve (default: 100000)
  --W W                    Fiber lattice width (default: 256)
  --a A                    Numerator of rational part of α (default: 1)
  --q Q                    Denominator of rational part of α (default: 6)
  --beta BETA              Irrational perturbation of α (default: 0.0)
  --weighted/--no-weighted Use log(p) weights for Φ (default: weighted)
  --topk TOPK              Number of top |F̂(k)| values to display (default: 5)
  --arc_width ARC_WIDTH    Half-width of major-arc window (default: 2)
  --samples SAMPLES        Number of random-beta samples (default: 0)
  --seed SEED              RNG seed for random sweep (default: 42)
  --selfcheck              Run smoke tests and exit
```

---

## Performance Notes

- The prime sieve runs in O(N log log N) time and O(N) memory.
- Φ evaluation is O(W · π(N)) where π(N) ≈ N/ln N is the prime-counting
  function.  For N = 10⁶ and W = 4096 this is ~280 M operations; expect a few
  seconds on a modern CPU.
- The FFT step is O(W log W) — negligible.
- Results are cached (per process) so re-running with the same (N, W, a, q,
  beta, weighted) tuple is instantaneous.

---

## Interpreting Output

```
============================================================
Goldbach Fiber-Sum Experiment
  N=100000, W=256, α=1/6 + 0
  Weighted (log p): True
  Primes ≤ N: 9592
============================================================

Top-5 |F̂(k)| magnitudes:
  # 1  k=   0  |F̂(k)|=2.7432e+11
  # 2  k= 128  |F̂(k)|=1.1201e+09
  ...

Goldbach target  k₀ = 200  (2N mod W)
  |F̂(k₀)| = 3.4817e+08

|Φ(α + t/W)| distribution (t = 0 … 255):
  min=1.2e+04  max=9.8e+05  mean=2.1e+05  median=1.9e+05  std=1.5e+05

Aliasing count (arc_width neighbours of k₀ above 2×median): 0
```

- **k=0 dominates** – expected; the DC component captures the total weight.
- **|F̂(k₀)|** – compare across different α, N, W to look for growth trends.
- **Aliasing count = 0** – clean spectrum around k₀.

---

## Unit Tests

```bash
python -m pytest tests/test_goldbach_experiment.py -v
```

See [`tests/test_goldbach_experiment.py`](../tests/test_goldbach_experiment.py).
