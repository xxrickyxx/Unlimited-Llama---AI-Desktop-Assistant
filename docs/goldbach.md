# Goldbach Experiment — Parameter Reference & Theory Guide

This document explains the mathematical conventions, CLI parameters, and
example configurations for the `goldbach_experiment.py` tooling included in
this repository.

---

## Background

**Goldbach's conjecture** states that every even integer greater than 2 can be
expressed as the sum of two prime numbers.  For example:

```
4  = 2 + 2
6  = 3 + 3
100 = 3 + 97 = 11 + 89 = 29 + 71 = 41 + 59 = 47 + 53
```

This tool approaches the conjecture experimentally using **wheel factorisation**
and a **Discrete Fourier Transform (DFT) on the reduced residue fiber** of a
wheel modulus W to produce spectral fingerprints for each even candidate N.

---

## Parameter Conventions

### N — Target Even Integer

`N` is the even integer you want to study.  It must satisfy:

- `N ≥ 4`
- `N % 2 == 0`

`N` drives every downstream computation:

| Derived quantity | Formula | Meaning |
|---|---|---|
| `k` | `k = N mod W` | Residue class of N in the wheel |
| Goldbach pairs | `{(p, N-p) : p prime, N-p prime}` | All prime-pair decompositions |
| DFT coefficient | `C(N,W,A)` (see below) | Spectral fingerprint of N on the fiber |

**CLI flag:** `--N <integer>`  
**Default:** `100`

```bash
python goldbach_experiment.py --N 200
```

---

### W — Wheel Modulus

`W` is the product of the first few primes (**primorial**).  It defines the
"wheel" used in the fiber DFT: only residues coprime to W (the *reduced
residue system* R(W)) are included in the fiber.

Common values:

| W | Factorisation | φ(W) |
|---|---|---|
| 6 | 2·3 | 2 |
| 30 | 2·3·5 | 8 |
| 210 | 2·3·5·7 | 48 |
| 2310 | 2·3·5·7·11 | 480 |

A larger W gives a finer-grained fiber and a sharper spectral resolution, at
the cost of more computation.

**CLI flag:** `--W <integer>`  
**Default:** `30`

```bash
python goldbach_experiment.py --N 100 --W 210
```

---

### A — Amplitude (Scaling Coefficient)

`A` is a real-valued scaling factor applied to every term in the fiber DFT
sum.  Its default value of `1.0` gives un-scaled coefficients.  Use `A > 1`
to amplify the signal or `A < 1` to attenuate it when comparing experiments
across different wheel sizes.

**CLI flag:** `--A <float>`  
**Default:** `1.0`

```bash
python goldbach_experiment.py --N 100 --W 30 --A 2.5
```

---

## The Fiber DFT and the k ≡ N (mod W) Coefficient

### Reduced Residue Fiber R(W)

The *fiber* for wheel modulus W is the set of all integers in [1, W) that are
coprime to W:

```
R(W) = { r : 1 ≤ r < W,  gcd(r, W) = 1 }
```

Its size is Euler's totient φ(W).  For W=30:

```
R(30) = {1, 7, 11, 13, 17, 19, 23, 29}   →   φ(30) = 8
```

### DFT Coefficient Formula

For a given N, W, and A, the Goldbach fiber DFT coefficient is:

```
C(N, W, A) = A · Σ_{r ∈ R(W)}  exp( 2πi · r · k / φ(W) )
```

where **k = N mod W** is the residue class of N modulo the wheel.

The role of k:

- `k` selects which "frequency bin" of the DFT spectrum the target N falls
  into.
- When k = 0 (N is a multiple of W), all additive characters align and the
  coefficient reaches its maximum absolute value A · φ(W).
- For other values of k, constructive and destructive interference produce a
  smaller |C|.
- The argument of C encodes which half of the fiber carries the Goldbach
  prime pairs.

### Example: N=100, W=30

```
k  = 100 mod 30 = 10
R(30) = {1, 7, 11, 13, 17, 19, 23, 29},   φ(30) = 8

C = 1.0 · Σ_r  exp(2πi · r · 10 / 8)
  = Σ_r  exp(2.5πi · r)
```

Run it yourself:

```bash
python goldbach_experiment.py --N 100 --W 30 --A 1.0 --spectrum
```

---

## CLI Reference

```
usage: goldbach_experiment [-h]
       [--N N | --N-start N_START] [--N-end N_END]
       [--W W] [--A A]
       [--output FILE] [--spectrum] [--pairs | --no-pairs]
```

| Flag | Type | Default | Description |
|---|---|---|---|
| `--N` | int | `100` | Single target even integer |
| `--N-start` | int | — | Start of batch range (even, ≥ 4) |
| `--N-end` | int | — | End of batch range (even, inclusive) |
| `--W` | int | `30` | Wheel modulus |
| `--A` | float | `1.0` | Amplitude / scaling coefficient |
| `--output` | path | — | Write CSV results to this file |
| `--spectrum` | flag | off | Print full fiber DFT spectrum |
| `--pairs` / `--no-pairs` | flag | on | Print/suppress Goldbach prime pairs |

---

## Example Configurations

### 1. Single target (defaults)

```bash
python goldbach_experiment.py
# N=100, W=30, A=1.0
```

### 2. Custom N and wheel

```bash
python goldbach_experiment.py --N 1000 --W 210 --A 1.0
```

### 3. Full spectral output for one N

```bash
python goldbach_experiment.py --N 60 --W 30 --A 1.0 --spectrum
```

### 4. Batch scan and CSV export

```bash
python goldbach_experiment.py --N-start 4 --N-end 200 --W 30 --A 1.0 \
    --output goldbach_4_200.csv
```

### 5. Larger wheel with amplified signal

```bash
python goldbach_experiment.py --N 2310 --W 2310 --A 0.1 --spectrum
```

### 6. Suppress pair listing (fast, DFT only)

```bash
python goldbach_experiment.py --N 500 --W 30 --no-pairs
```

---

## Output Fields (CSV / Summary)

| Field | Description |
|---|---|
| `N` | Target even integer |
| `W` | Wheel modulus |
| `A` | Amplitude |
| `k (= N mod W)` | Residue class of N modulo W |
| `phi(W)` | Euler's totient of W (fiber size) |
| `DFT_real` | Real part of the DFT coefficient |
| `DFT_imag` | Imaginary part |
| `DFT_abs` | Absolute value \|C(N,W,A)\| |
| `num_pairs` | Number of Goldbach prime pairs found |

---

## Further Reading

- Hardy, G. H. & Littlewood, J. E. (1923). "Some Problems of 'Partitio Numerorum'; III: On the Expression of a Number as a Sum of Primes." *Acta Mathematica* 44, 1–70.
- Goldbach, C. (1742). Letter to Euler, 7 June 1742.
- Apostol, T. M. (1976). *Introduction to Analytic Number Theory*. Springer.
