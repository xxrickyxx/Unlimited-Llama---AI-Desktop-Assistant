# Goldbach Fiber-Sum Experiments

This document describes the mathematical framework, conventions, and practical
usage of `tools/goldbach_fiber_experiment.py`.

---

## 1. Definitions and Conventions

| Symbol | Meaning |
|--------|---------|
| **N** | Target even integer whose Goldbach representation is studied |
| **p** | A prime number ≤ N |
| **log p** | Natural logarithm of *p*, used as the von-Mangoldt-style weight |
| **α** | Real frequency parameter, sampled uniformly in [0, 1) |
| **W** | Fiber width and DFT modulus (must be a positive integer; default 1024) |
| **A** | Major-arc exponent (default 6) |
| **P** | Major-arc threshold P = (log N)^A |

### 1.1 Exponential sum Phi

```
Phi(α) = Σ_{p ≤ N}  log(p) · exp(2πi α p)
```

### 1.2 Major and minor arcs

A frequency α ∈ [0, 1) is on a **major arc** if there exist integers *q*, *a*
with

```
1 ≤ q ≤ P,   gcd(a, q) = 1,   |α − a/q| ≤ P / (N · q)
```

Otherwise α is on a **minor arc**.  The experiment samples α uniformly from
minor arcs using rejection sampling.

### 1.3 Fiber vector

For a fixed minor-arc α the **fiber vector** is

```
F_α(t) = Phi(α + t/W)^2,    t = 0, 1, …, W−1
```

### 1.4 Fiber DFT coefficient

The DFT of F_α at frequency index *k* is

```
F̂_α(k) = Σ_{t=0}^{W−1}  F_α(t) · exp(−2πi k t / W)
```

The **key index** is **k\* = N mod W**.  By a standard character-sum identity,
F̂_α(k\*) concentrates the contribution of prime pairs (p₁, p₂) satisfying
p₁ + p₂ ≡ N (mod W).

---

## 2. Output Metrics

For each sampled minor-arc α the tool reports:

| Metric | Description |
|--------|-------------|
| `B(α)` | **Aliasing count** — number of fiber shifts α + t/W (t = 0…W−1) that land on a major arc |
| `\|F̂(k*)\|` | Magnitude of the fiber DFT at k\* = N mod W |
| `top-K` | The K largest `\|F̂(k)\|` values and their indices |
| `\|Phi\|` mean / std / max | Distribution statistics of `\|Phi(α + t/W)\|` across the fiber |

---

## 3. Running the Script

### Prerequisites

```bash
pip install numpy
```

### Default run (N = 1 000 000, W = 1024, 50 samples)

```bash
python tools/goldbach_fiber_experiment.py
```

### Override parameters

```bash
python tools/goldbach_fiber_experiment.py --N 100000 --W 512 --A 5 --samples 20 --seed 42
```

| Flag | Default | Description |
|------|---------|-------------|
| `--N` | 1000000 | Target even integer |
| `--W` | 1024 | Fiber width / DFT modulus |
| `--A` | 6 | Major-arc exponent |
| `--samples` | 50 | Number of minor-arc α samples |
| `--seed` | *none* | Random seed for reproducibility |
| `--top-k` | 5 | Number of top \|F̂(k)\| values to print |
| `--csv` | *none* | Write results to the given CSV file |
| `--self-check` | *off* | Run identity-verification checks then exit |

> **Note on the A parameter.**  With the natural-logarithm convention the
> major-arc threshold is P = (ln N)^A.  For N = 10⁶ the condition P < N
> (needed for any minor arcs to exist) requires A < ln N / ln(ln N)
> (which equals approximately 5.3 for N = 10⁶).
> The default A = 6 is intentionally generous (all arcs are major) and is
> suited to studying the near-complete-coverage regime; for typical minor-arc
> sampling use A = 2 or A = 3.

### Save results to CSV

```bash
python tools/goldbach_fiber_experiment.py --N 100000 --W 256 --samples 10 --csv results.csv
```

### Self-check mode (validates correctness for small N/W)

```bash
python tools/goldbach_fiber_experiment.py --N 100 --W 16 --self-check
```

The self-check verifies two identities:

1. **FFT vs direct DFT** — confirms that `numpy.fft.fft` and the explicit
   DFT summation agree to machine precision.
2. **F̂(k\*) vs prime-pair sum** — confirms that the fiber DFT coefficient
   matches the weighted count of prime pairs (p₁, p₂) with p₁ + p₂ ≡ N (mod W).

---

## 4. Interpreting Results

### Spikes vs flatness in |F̂(k)|

* A **spike at k\*** (i.e., `|F̂(k*)|` significantly larger than the median
  `|F̂(k)|`) indicates that prime pairs congruent to N mod W are abundant for
  this α — consistent with the Goldbach conjecture holding modulo W.
* A **flat spectrum** (no dominant spike) is expected for most minor-arc α
  values; it means the exponential sum behaves like random noise, which is
  the desired "cancellation" property.

### Aliasing count B(α)

* A small `B(α)` (close to 0) means that almost none of the fiber shifts
  land on a major arc.  This is the typical minor-arc regime.
* A large `B(α)` suggests that the fiber passes near rational approximations,
  which can inflate `|F̂(k*)|` artificially (aliasing).

### |Phi| distribution

* Large mean / max of `|Phi(α + t/W)|` on a minor arc is unusual and may
  indicate a near-major-arc sample that slipped through the rejection filter
  (tighten `--A` or reduce `--W` to diagnose).
* Roughly uniform `|Phi|` values across t confirm that α is well into the
  minor-arc regime.

---

## 5. Performance Notes

* Runtime scales approximately as O(π(N) · W) per sample, where π(N) is the
  number of primes up to N.
* For N = 10⁶ and W = 1024 expect roughly an order-of-magnitude of seconds
  per sample depending on hardware (NumPy vectorised inner loop).
* For W > 4096 or N > 10⁶ consider reducing `--samples`.
