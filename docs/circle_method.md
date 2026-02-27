# Hardy-Littlewood Circle Method: Definitions and Conventions

This document defines the exponential sums and related quantities used in the
circle-method analysis of Goldbach-type problems, establishes the convention for
N, specifies how the fiber/DFT decomposition uses the congruence k ≡ N (mod W),
and explains how to translate between the Λ-weighted sum and the prime-only
log-weighted sum Φ(α).

---

## 1. Notation

| Symbol | Meaning |
|--------|---------|
| N | A fixed **even** positive integer — the target for the Goldbach representation |
| p, q | Prime numbers |
| Λ | Von Mangoldt function: Λ(n) = log p if n = p^k (k ≥ 1), else 0 |
| e(x) | e(x) = exp(2πix) |
| W | Wheel modulus: product of all primes up to a chosen bound B, W = ∏_{p ≤ B} p |
| k | A residue class mod W satisfying k ≡ N (mod W) |

---

## 2. Definition of Φ(α) — Option B (Prime-Only Log-Weighted Sum)

**Definition (Option B).**

```
Φ(α) = ∑_{p ≤ N} (log p) · e(α p)
```

where:
- the sum runs over all **primes** p with p ≤ N,
- N is the target even integer (see Section 1),
- e(αp) = exp(2πi α p).

This is a prime-only, log-weighted exponential sum.  Only genuine primes
contribute; prime-power terms (p², p³, …) are excluded.

---

## 3. Λ-Weighted Exponential Sum S(α)

The von Mangoldt sum, standard in the circle method, is:

```
S(α) = ∑_{n ≤ N} Λ(n) · e(α n)
```

Because Λ(n) = log p when n = p^k and zero otherwise, S(α) includes both
primes and higher prime powers.

---

## 4. Relation Between Φ(α) and S(α)

Decomposing S(α) by the prime-power order k gives:

```
S(α) = Φ(α) + R(α)
```

where the **prime-power remainder** is:

```
R(α) = ∑_{k ≥ 2} ∑_{p : p^k ≤ N} (log p) · e(α p^k)
```

### Size of R(α)

The number of prime powers p^k with k ≥ 2 and p^k ≤ N is O(N^{1/2}), so:

```
|R(α)| ≤ ∑_{k ≥ 2} ∑_{p^k ≤ N} log p  =  O(N^{1/2})
```

In contrast, Φ(α) and S(α) are both O(N / log N) on average, so R(α) is a
lower-order error term of relative size O(N^{-1/2} log N).

---

## 5. Goldbach Representation Counts

### 5.1 Λ-Weighted Count r_Λ(N)

Define the **Λ-weighted Goldbach sum**:

```
r_Λ(N) = ∑_{n=1}^{N-1} Λ(n) · Λ(N − n)
```

This counts representations N = m + n weighted by Λ(m)Λ(n), so it receives
its main contribution from prime pairs (p, q) with p + q = N, each contributing
(log p)(log q).

### 5.2 Prime-Pair Count r(N)

The ordinary Goldbach count is:

```
r(N) = #{(p, q) : p + q = N, p prime, q prime}
```

(Ordered pairs; if N = 2p then the pair (p, p) is counted once with multiplicity 1
under the convention that unordered pairs are used, or twice if ordered.)

### 5.3 Approximate Relation

Because Λ(p) = log p and all primes near N satisfy log p ≈ log N:

```
r_Λ(N) ≈ (log N)² · r(N)          [rough approximation]
```

More precisely, the **log-weighted prime-pair sum** is:

```
∫₀¹ Φ(α)² e(−α N) dα  =  ∑_{p + q = N} (log p)(log q)
```

and this equals exactly (log N)² · r(N) only when all contributing primes have
exactly log p = log N, which is the leading-order approximation.  The precise
relation is:

```
r(N) = ∑_{p + q = N} 1
     = ∑_{p + q = N} (log p)(log q) / (log p)(log q)

r_Λ(N) − ∑_{p + q = N} (log p)(log q)  =  O(r(N) · N^{1/2})   [prime-power cross terms]
```

---

## 6. Circle Method Identities

### 6.1 Exact Identity for r_Λ(N)

By the orthogonality of additive characters:

```
r_Λ(N)  =  ∫₀¹ S(α)² · e(−α N) dα
```

This is the **fundamental circle method identity**.

### 6.2 Identity for the Log-Weighted Prime-Pair Sum

Using Φ(α) in place of S(α):

```
∫₀¹ Φ(α)² · e(−α N) dα  =  ∑_{p + q = N} (log p)(log q)
```

### 6.3 Approximation for r(N)

Combining the two identities with the expansion S = Φ + R:

```
r_Λ(N)  =  ∫₀¹ Φ(α)² e(−α N) dα
          + 2 ∫₀¹ Φ(α) R(α) e(−α N) dα
          + ∫₀¹ R(α)² e(−α N) dα
```

The last two integrals are lower-order error terms (O(N^{1/2}) each in absolute
value after applying Parseval-type bounds on R).  Hence:

```
r_Λ(N)  =  ∑_{p + q = N} (log p)(log q)  +  O(N^{1/2} log² N)
```

Dividing both sides by (log N)²:

```
r(N)  ≈  r_Λ(N) / (log N)²   [up to lower-order corrections]
```

---

## 7. Fiber / DFT Decomposition with k ≡ N (mod W)

### 7.1 Wheel Sieve and Residue Classes

Choose W = ∏_{p ≤ B} p (e.g., W = 2·3·5·7 = 210 for B = 7).  Every prime
p > B satisfies gcd(p, W) = 1, so p lies in one of the φ(W) reduced residue
classes mod W.

### 7.2 Fiber Constraint for Goldbach Pairs

For a Goldbach pair (p, q) with p + q = N, write:

```
p ≡ a (mod W),    q ≡ N − a (mod W)
```

for some a with gcd(a, W) = 1 and gcd(N − a, W) = 1.

Since p + q = N, the residues (a, N − a mod W) always satisfy:

```
a + (N − a)  ≡  N  (mod W)
```

**Convention:** we label a fiber by the residue k = a of the *first* prime p,
so k ≡ a (mod W).  Because we need N − k ≡ N − a (mod W) to also be coprime
to W, the relevant fibers satisfy:

```
k ≡ N (mod W)   is NOT required individually on k,
```

but the **pair** (k, N − k mod W) must consist of reduced residues.  However,
a common simplification used in DFT-based computations is to restrict to fibers
with k ≡ N/2 (mod W) when N/2 is an integer — or more generally to fix k such
that the DFT bin index k satisfies:

```
k  ≡  N  (mod W)        [DFT bin / fiber index convention]
```

This means the DFT is evaluated at frequency bin k = N mod W, which is the
unique bin (within 0 ≤ k < W) for which the exponential e(−αN) "lines up"
with the fiber after the wheel-sieve reduction.

### 7.3 DFT Formula on a Single Fiber

Let P_a = {p prime : p ≤ N, p ≡ a (mod W)}.  Define the fiber exponential sum:

```
Φ_a(α) = ∑_{p ∈ P_a} (log p) · e(α p)
```

Then:

```
Φ(α) = ∑_{a : gcd(a,W)=1} Φ_a(α)        [sum over φ(W) fibers]
```

To compute Φ_a(α) efficiently via DFT, sample α at points α_j = j / (W·M)
for j = 0, …, W·M − 1 (M controls resolution).  The bin corresponding to N is
the bin index:

```
j_N  =  N mod (W · M)
```

which satisfies j_N ≡ N (mod W) after reducing j_N modulo W.  This is the
**k ≡ N (mod W)** condition for the DFT output bin.

---

## 8. Translating Between Λ-Weighted and Prime-Only Log-Weighted Sums

The following steps convert between the two formulations:

### Step 1 — From S(α) to Φ(α)

Given S(α) computed over all n ≤ N with Λ(n):

```
Φ(α) = S(α) − R(α)
```

where R(α) = ∑_{k≥2, p^k ≤ N} (log p) e(α p^k).  In practice R(α) is computed
by iterating over prime squares, cubes, etc., up to N (there are only O(√N) such
terms).

### Step 2 — From r_Λ(N) to the Log-Weighted Prime-Pair Sum

```
∑_{p+q=N} (log p)(log q)  =  r_Λ(N)
                             − 2 ∑_{k≥2} ∑_{p^k ≤ N} (log p) · Λ(N − p^k)
                             − ∑_{j,k≥2} ∑_{p^j + q^k = N} (log p)(log q)
```

The second and third lines are the **cross-term** and **prime-power-pair**
corrections, both O(N^{1/2} log² N).

### Step 3 — From the Log-Weighted Sum to r(N)

```
r(N) = ∑_{p+q=N} (log p)(log q) / (log p · log q)
     = ∑_{p+q=N} 1
```

There is no direct algebraic formula expressing r(N) purely in terms of r_Λ(N)
without knowing individual prime values.  The approximation

```
r(N)  ≈  r_Λ(N) / (log N)²
```

is accurate up to a relative error of O(log log N / log N), arising from the
spread of log p around log N for primes p ≈ N/2.

A sharper estimate uses partial summation: if π(x; a, W) denotes primes in
the arithmetic progression a (mod W), then for each fiber:

```
∑_{p+q=N, p≡a(W)} (log p)(log q)  ≈  (log N)² · #{p ≡ a (W) : p+q=N, q prime}
```

and summing over a gives the improved approximation.

---

## 9. Summary of Key Formulas

| Quantity | Formula |
|----------|---------|
| Φ(α) | ∑_{p ≤ N} (log p) e(αp) — prime-only log-weighted sum |
| S(α) | ∑_{n ≤ N} Λ(n) e(αn) — von Mangoldt sum (includes prime powers) |
| R(α) | S(α) − Φ(α) — prime-power remainder, O(N^{1/2}) |
| r_Λ(N) | ∫₀¹ S(α)² e(−αN) dα = ∑_n Λ(n) Λ(N−n) |
| log-wt sum | ∫₀¹ Φ(α)² e(−αN) dα = ∑_{p+q=N} (log p)(log q) |
| r(N) | ≈ r_Λ(N) / (log N)² — ordinary prime-pair count |
| DFT bin | k ≡ N (mod W) — fiber index for the Goldbach frequency bin |

---

## References

- G.H. Hardy and J.E. Littlewood, "Some problems of 'Partitio Numerorum' III:
  On the expression of a number as a sum of primes," *Acta Math.* 44 (1923),
  1–70.
- I.M. Vinogradov, *The Method of Trigonometrical Sums in the Theory of Numbers*,
  Dover, 2004.
- H.L. Montgomery and R.C. Vaughan, *Multiplicative Number Theory I*, Cambridge
  University Press, 2006.
