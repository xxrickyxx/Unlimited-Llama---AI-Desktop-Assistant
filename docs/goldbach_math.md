# Goldbach Fiber-Sum Experiment — Mathematical Background

## Overview

This experiment probes **Goldbach's conjecture** (every even integer > 2 is the
sum of two primes) using the **Hardy–Littlewood circle method**.  The key object
is the exponential sum

```
Phi(alpha) = sum_{p <= N, p prime} log(p) * exp(2*pi*i * alpha * p)
```

where `log` denotes the natural logarithm and `N` is the prime upper bound.

The weighted count of representations of an even integer `n` as `p1 + p2`
(with `p1, p2 <= N`) is related to the Fourier coefficient

```
r(n) = integral_0^1  Phi(alpha)^2 * exp(-2*pi*i * alpha * n)  d(alpha).
```

The circle method partitions `[0, 1)` into **major arcs** (neighbourhoods of
rational points) and **minor arcs** (the complement) to estimate this integral.

---

## Parameters

| Symbol | Name           | Default     | Description                                      |
|--------|----------------|-------------|--------------------------------------------------|
| N      | prime bound    | 1 000 000   | Upper bound for the prime sum; must be even      |
| W      | window size    | 1024        | Number of frequency bins; power-of-2 recommended |
| A      | arc exponent   | 6           | Controls P = (ln N)^A                            |
| P      | arc threshold  | (ln N)^A    | Denominator bound that defines major arcs        |

---

## Major and Minor Arcs

**Major arc** M(a, q): the set of alpha in [0, 1) satisfying

```
|alpha - a/q| <= P / (N * q)
```

for integers `a`, `q` with `1 <= q <= P` and `gcd(a, q) = 1`.

**Minor arcs**: the complement `[0, 1) \ union_{a, q} M(a, q)`.

On major arcs, `Phi(alpha)` is well-approximated by explicit Gauss sums (the
"singular series"); on minor arcs it is small in absolute value (quantified by
Weyl/Vinogradov estimates), making the minor-arc integral negligible in the
asymptotic expansion of `r(n)`.

### Implementation note

For a given `alpha`, the major-arc test is equivalent to asking whether there
exists `q` in `{1, ..., P}` such that

```
dist(q * alpha, Z) <= P * q / N
```

where `dist(x, Z)` is the distance from `x` to the nearest integer.  This form
avoids enumerating all fractions `a/q` and allows a vectorised implementation:
for each `q`, the nearest integer `a = round(q * alpha)` gives the minimising
fraction.  Fractions with `gcd(a, q) > 1` are covered by smaller denominators
already checked, so no explicit gcd test is needed.

---

## The Fiber-Sum and FFT Analysis

Fix an alpha (sampled from the minor arcs).  Define the **fiber**

```
F(t) = Phi(alpha + t/W)^2,   t = 0, 1, ..., W-1.
```

The discrete Fourier transform

```
Fhat(k) = sum_{t=0}^{W-1}  F(t) * exp(-2*pi*i * k * t / W)
```

satisfies

```
Fhat(k)  ~  sum_{p1, p2 <= N,  p1+p2 ≡ k (mod W)}
               log(p1) * log(p2) * exp(2*pi*i * alpha * (p1 + p2)).
```

Thus `|Fhat(k)|` measures the weighted count of prime pairs whose sum falls in
residue class `k mod W`, modulated by the exponential twist `exp(2*pi*i*alpha*)`.

The bin `k = N mod W` is singled out because it corresponds to the Goldbach
target `N`; the value `|Fhat(k_N)|` is the main observable of interest.

---

## Efficient Computation via the g-Vector

Direct evaluation of `Phi(alpha + t/W)` for each `t` costs O(W * pi(N))
operations (pi(N) ≈ N / ln N primes).  Instead, observe:

```
Phi(alpha + t/W) = sum_{k=0}^{W-1}  g[k] * exp(2*pi*i * k * t / W)
```

where

```
g[k] = sum_{p <= N,  p ≡ k (mod W)}  log(p) * exp(2*pi*i * alpha * p).
```

Hence the full array `{Phi(alpha + t/W) : t = 0,...,W-1}` is obtained by a
single **Inverse FFT** of the W-vector `g` (scaled by W), reducing the cost to
`O(pi(N) + W log W)`.  For N = 10^6 and W = 1024 this is roughly 400x faster
than the naive approach.

### Code path

```python
# 1. Per-prime complex weights
c[p] = log(p) * exp(2*pi*i * alpha * p)

# 2. Accumulate into residue-class bins
g[k] = sum of c[p] for all p with p % W == k

# 3. Single IFFT gives all Phi values
phi_vals = ifft(g) * W          # shape (W,)

# 4. Square and FFT
F    = phi_vals ** 2
Fhat = fft(F)                   # shape (W,)
```

---

## Aliasing Count B(alpha)

```
B(alpha) = #{t in {0,...,W-1} : alpha + t/W  mod 1  lies in major arcs}
```

A large `B(alpha)` means many of the `W` fiber points fall near rational
approximations, so major-arc contributions contaminate the minor-arc signal.
For a well-chosen minor-arc alpha with `P << N`, `B(alpha)` should be small
(often 0 or a handful of points out of W = 1024).

---

## Interpreting the Output

| Output field        | Meaning                                                      |
|---------------------|--------------------------------------------------------------|
| `alpha`             | Sampled point (minor arc if feasible, otherwise random)      |
| `P`                 | Major-arc threshold (ln N)^A                                 |
| `B(alpha)`          | Number (and fraction) of fiber points in major arcs          |
| `k_N = N mod W`     | Target Goldbach residue bin                                  |
| `|Fhat[k_N]|`       | Fiber-sum magnitude at the Goldbach bin                      |
| Top-K `|Fhat[k]|`   | Dominant residue classes of prime-pair sums                  |
| `mean / std` of `|Fhat|` | Overall spectral energy and spread                      |

A **Goldbach-consistent** experiment will show `|Fhat[k_N]|` clearly elevated
relative to the mean, consistent with the abundance of prime pairs summing to N.

---

## Notes on Default Parameters

With the defaults `N = 10^6`, `A = 6`:

- `P = (ln 10^6)^6 ≈ 6.96 × 10^6 > N`.
- Every alpha lies in a major arc; genuine minor-arc sampling is not possible.
- The script falls back to a random alpha and prints a warning.

**To run a genuine minor-arc experiment** with N = 1 000 000:

```bash
python goldbach_experiment.py --A 2    # P ≈ 191  (narrow major arcs)
python goldbach_experiment.py --A 3    # P ≈ 2638 (moderate major arcs)
```

Alternatively, keep A = 6 and use a much larger N.  For N = 10^12 the threshold
`(ln 10^12)^6 ≈ 2.7 × 10^9 << 10^12`, but sieving and storing all primes up to
10^12 requires significant RAM (a segmented sieve or prime-counting approach
would be needed).

The default A = 6 is chosen to match the standard Siegel–Walfisz exponent used
in the literature; it illustrates the full API and all computation steps without
requiring a very large N.

---

## References

- Hardy, G. H. & Littlewood, J. E. (1923). *Some problems of 'partitio numerorum' III:
  On the expression of a number as a sum of primes.* Acta Mathematica, 44, 1–70.
- Vaughan, R. C. (1997). *The Hardy–Littlewood Method* (2nd ed.). Cambridge University Press.
- Tao, T. (2012). *Every odd number greater than 1 is the sum of at most five primes.*
  Mathematics of Computation, 83(286), 997–1038.
