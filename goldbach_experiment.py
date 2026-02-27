#!/usr/bin/env python3
"""
Goldbach Conjecture – Fiber-Sum / Discrete Fourier-Shift Experiment
====================================================================
Numerical tooling for the fiber-sum approach to Goldbach-type problems.

Method overview
---------------
Given N, the width W, and a rational approximation α = a/q + β:

1. Sieve primes p ≤ N with the Sieve of Eratosthenes.
2. Build the exponential sum
       Φ(θ) = Σ_{p≤N} log(p) · e(θ·p),   e(x) = exp(2πi·x)
   evaluated at the W lattice points θ_t = α + t/W,  t = 0,…,W-1.
3. Form the fiber-square  F(t) = Φ(θ_t)².
4. Compute F̂(k) = FFT(F) and report top-|k| magnitudes.
5. Report |F̂(k₀)| where k₀ ≡ 2N (mod W) – the Goldbach target frequency.
6. Summarise the |Φ| distribution and count "aliasing" events where
   |F̂(k)| > threshold for k in [k₀-w, k₀+w] excluding k₀.

Caching
-------
Φ evaluations are cached so that repeated calls with identical (N, W, a, q,
beta) parameters reuse previously computed arrays without re-sieving.

Usage
-----
    python goldbach_experiment.py --N 100000 --W 256 --a 1 --q 6 --beta 0.0
    python goldbach_experiment.py --selfcheck
    python goldbach_experiment.py --N 500000 --W 1024 --a 1 --q 6 \\
        --samples 50 --topk 10 --arc_width 3

Run  python goldbach_experiment.py --help  for all options.
"""

from __future__ import annotations

import argparse
import sys
from functools import lru_cache
from typing import Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Prime sieve
# ---------------------------------------------------------------------------

def sieve_primes(N: int) -> np.ndarray:
    """Return all primes p with 2 ≤ p ≤ N as a NumPy integer array."""
    if N < 2:
        return np.array([], dtype=np.int64)
    is_prime = np.ones(N + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(N**0.5) + 1):
        if is_prime[i]:
            is_prime[i * i :: i] = False
    return np.nonzero(is_prime)[0].astype(np.int64)


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_phi(primes: np.ndarray, thetas: np.ndarray, weighted: bool = True) -> np.ndarray:
    """
    Compute Φ(θ_t) for each θ_t in *thetas*.

    If *weighted* is True (default), use Λ-weights log(p) per prime.
    Otherwise use unit weights (unweighted exponential sum).

    Parameters
    ----------
    primes  : 1-D int array of primes ≤ N
    thetas  : 1-D float array of evaluation points (length W)
    weighted: bool

    Returns
    -------
    phi : complex128 array of length W
    """
    weights = np.log(primes) if weighted else np.ones(len(primes), dtype=np.float64)
    # phi[t] = sum_p  weights[p] * exp(2πi * theta_t * p)
    # Vectorised: outer product (W x P) then sum over primes axis
    phases = 2.0 * np.pi * np.outer(thetas, primes)   # shape (W, P)
    phi = (weights * np.exp(1j * phases)).sum(axis=1)  # shape (W,)
    return phi


# Cached wrapper keyed on immutable parameters
@lru_cache(maxsize=128)
def _cached_phi_key(N: int, W: int, a: int, q: int, beta: float, weighted: bool
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return (primes, thetas, phi) for the given parameters.
    The result is cached so identical calls skip re-computation.
    """
    primes = sieve_primes(N)
    alpha = a / q + beta
    thetas = alpha + np.arange(W, dtype=np.float64) / W
    phi = compute_phi(primes, thetas, weighted=weighted)
    return primes, thetas, phi


def run_experiment(
    N: int,
    W: int,
    a: int,
    q: int,
    beta: float,
    weighted: bool = True,
    topk: int = 5,
    arc_width: int = 2,
    verbose: bool = True,
) -> dict:
    """
    Run a single fiber-sum experiment and return a results dict.

    Parameters
    ----------
    N         : upper bound for primes
    W         : width of the fiber lattice (number of evaluation points)
    a, q      : rational part of α = a/q + beta
    beta      : irrational perturbation
    weighted  : use log(p) weights for Φ
    topk      : number of top |F̂(k)| values to report
    arc_width : half-width of the "major arc" window around k ≡ 2N (mod W)
    verbose   : print results to stdout

    Returns
    -------
    dict with keys: primes, thetas, phi, F, F_hat, top_k_indices,
                    top_k_magnitudes, goldbach_k, goldbach_magnitude,
                    aliasing_count, phi_stats
    """
    primes, thetas, phi = _cached_phi_key(N, W, a, q, beta, weighted)

    # Fiber-square
    F = phi ** 2                         # shape (W,)

    # FFT-based DFT
    F_hat = np.fft.fft(F)               # shape (W,)
    magnitudes = np.abs(F_hat)

    # Top-k frequencies
    top_idx = np.argsort(magnitudes)[::-1][:topk]
    top_mag  = magnitudes[top_idx]

    # Goldbach target frequency: k₀ ≡ 2N (mod W)
    k0 = (2 * N) % W
    goldbach_mag = magnitudes[k0]

    # Aliasing count: how many k in the arc window (excluding k₀) exceed
    # the median magnitude by a factor of 2
    threshold = 2.0 * np.median(magnitudes)
    arc_ks = [(k0 + dk) % W for dk in range(-arc_width, arc_width + 1) if dk != 0]
    aliasing_count = sum(1 for k in arc_ks if magnitudes[k] > threshold)

    # |Φ| distribution statistics
    phi_abs = np.abs(phi)
    phi_stats = {
        "min":    float(phi_abs.min()),
        "max":    float(phi_abs.max()),
        "mean":   float(phi_abs.mean()),
        "median": float(np.median(phi_abs)),
        "std":    float(phi_abs.std()),
    }

    results = {
        "N": N, "W": W, "a": a, "q": q, "beta": beta, "weighted": weighted,
        "primes": primes,
        "thetas": thetas,
        "phi": phi,
        "F": F,
        "F_hat": F_hat,
        "top_k_indices": top_idx,
        "top_k_magnitudes": top_mag,
        "goldbach_k": k0,
        "goldbach_magnitude": float(goldbach_mag),
        "aliasing_count": aliasing_count,
        "phi_stats": phi_stats,
    }

    if verbose:
        _print_results(results)

    return results


def _print_results(r: dict) -> None:
    n_primes = len(r["primes"])
    print(f"\n{'='*60}")
    print(f"Goldbach Fiber-Sum Experiment")
    print(f"  N={r['N']}, W={r['W']}, α={r['a']}/{r['q']} + {r['beta']:.6g}")
    print(f"  Weighted (log p): {r['weighted']}")
    print(f"  Primes ≤ N: {n_primes}")
    print(f"{'='*60}")

    print(f"\nTop-{len(r['top_k_indices'])} |F̂(k)| magnitudes:")
    for rank, (k, m) in enumerate(zip(r["top_k_indices"], r["top_k_magnitudes"]), 1):
        marker = " ← Goldbach k₀" if k == r["goldbach_k"] else ""
        print(f"  #{rank:2d}  k={k:5d}  |F̂(k)|={m:.4e}{marker}")

    print(f"\nGoldbach target  k₀ = {r['goldbach_k']}  (2N mod W)")
    print(f"  |F̂(k₀)| = {r['goldbach_magnitude']:.4e}")

    s = r["phi_stats"]
    print(f"\n|Φ(α + t/W)| distribution (t = 0 … {r['W']-1}):")
    print(f"  min={s['min']:.4e}  max={s['max']:.4e}")
    print(f"  mean={s['mean']:.4e}  median={s['median']:.4e}  std={s['std']:.4e}")

    print(f"\nAliasing count (arc_width neighbours of k₀ above 2×median): {r['aliasing_count']}")
    print()


# ---------------------------------------------------------------------------
# Random-sample sweep
# ---------------------------------------------------------------------------

def random_sample_sweep(
    N: int,
    W: int,
    a: int,
    q: int,
    n_samples: int,
    topk: int = 5,
    arc_width: int = 2,
    seed: int = 42,
    verbose: bool = True,
) -> list[dict]:
    """
    Run *n_samples* experiments with beta drawn uniformly from [0, 1/q).

    Returns a list of result dicts (verbose=False for each individual run;
    a summary is printed if verbose=True).
    """
    rng = np.random.default_rng(seed)
    betas = rng.uniform(0.0, 1.0 / max(q, 1), size=n_samples)

    results = []
    for i, beta in enumerate(betas):
        r = run_experiment(N, W, a, q, float(beta),
                           topk=topk, arc_width=arc_width, verbose=False)
        results.append(r)

    if verbose:
        goldbach_mags = [r["goldbach_magnitude"] for r in results]
        aliasing_counts = [r["aliasing_count"] for r in results]
        print(f"\n{'='*60}")
        print(f"Random-sample sweep: {n_samples} samples")
        print(f"  N={N}, W={W}, a/q={a}/{q}")
        print(f"  |F̂(k₀)| – mean={np.mean(goldbach_mags):.4e} "
              f"std={np.std(goldbach_mags):.4e} "
              f"max={np.max(goldbach_mags):.4e}")
        print(f"  Aliasing count – mean={np.mean(aliasing_counts):.2f} "
              f"max={np.max(aliasing_counts)}")
        print(f"{'='*60}\n")

    return results


# ---------------------------------------------------------------------------
# Self-check / smoke test
# ---------------------------------------------------------------------------

def self_check() -> bool:
    """
    Run a small self-check to verify correctness of core routines.

    Returns True if all checks pass, False otherwise.
    Prints a summary of each check.
    """
    ok = True

    # --- 1. Prime sieve ---
    expected_primes_30 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    got = sieve_primes(30).tolist()
    if got == expected_primes_30:
        print("[PASS] sieve_primes(30)")
    else:
        print(f"[FAIL] sieve_primes(30): expected {expected_primes_30}, got {got}")
        ok = False

    # --- 2. Phi at theta=0 should equal sum of log(p) ---
    primes = sieve_primes(30)
    phi0 = compute_phi(primes, np.array([0.0]))[0]
    expected = float(np.log(primes).sum())
    if abs(phi0.real - expected) < 1e-9 and abs(phi0.imag) < 1e-9:
        print("[PASS] Φ(0) == Σ log(p)")
    else:
        print(f"[FAIL] Φ(0): expected {expected}, got {phi0}")
        ok = False

    # --- 3. FFT round-trip ---
    x = np.array([1.0, 2.0, 3.0, 4.0])
    if np.allclose(np.fft.ifft(np.fft.fft(x)), x):
        print("[PASS] FFT round-trip")
    else:
        print("[FAIL] FFT round-trip")
        ok = False

    # --- 4. Small end-to-end experiment ---
    r = run_experiment(50, 8, 1, 6, 0.0, verbose=False)
    if isinstance(r["goldbach_magnitude"], float) and r["goldbach_magnitude"] >= 0:
        print("[PASS] end-to-end run_experiment (N=50, W=8)")
    else:
        print(f"[FAIL] end-to-end: {r['goldbach_magnitude']}")
        ok = False

    # --- 5. Goldbach k ≡ 2N (mod W) ---
    N, W = 50, 8
    expected_k = (2 * N) % W
    if r["goldbach_k"] == expected_k:
        print(f"[PASS] goldbach_k == (2N mod W) = {expected_k}")
    else:
        print(f"[FAIL] goldbach_k: expected {expected_k}, got {r['goldbach_k']}")
        ok = False

    # --- 6. Parseval-type check: sum |F̂(k)|² == W * sum |F(t)|² ---
    F = r["F"]
    F_hat = r["F_hat"]
    lhs = np.sum(np.abs(F_hat) ** 2)
    rhs = W * np.sum(np.abs(F) ** 2)
    if np.isclose(lhs, rhs, rtol=1e-9):
        print("[PASS] Parseval's identity for F̂")
    else:
        print(f"[FAIL] Parseval: |F̂|²={lhs:.6e}, W·|F|²={rhs:.6e}")
        ok = False

    print()
    if ok:
        print("All self-checks passed.")
    else:
        print("One or more self-checks FAILED.")
    return ok


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Goldbach Fiber-Sum / Discrete Fourier-Shift Experiment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--N", type=int, default=100_000,
                   help="Upper bound for prime sieve")
    p.add_argument("--W", type=int, default=256,
                   help="Fiber lattice width (number of evaluation points)")
    p.add_argument("--a", type=int, default=1,
                   help="Numerator of rational part of α = a/q + beta")
    p.add_argument("--q", type=int, default=6,
                   help="Denominator of rational part of α")
    p.add_argument("--beta", type=float, default=0.0,
                   help="Irrational perturbation of α")
    p.add_argument("--weighted", action=argparse.BooleanOptionalAction, default=True,
                   help="Use log(p) weights for Φ (--no-weighted for unit weights)")
    p.add_argument("--topk", type=int, default=5,
                   help="Number of top |F̂(k)| values to display")
    p.add_argument("--arc_width", type=int, default=2,
                   help="Half-width of major-arc window for aliasing count")
    p.add_argument("--samples", type=int, default=0,
                   help="If > 0, run this many random-beta samples instead of a single run")
    p.add_argument("--seed", type=int, default=42,
                   help="RNG seed for random-sample sweep")
    p.add_argument("--selfcheck", action="store_true",
                   help="Run self-check / smoke tests and exit")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.selfcheck:
        passed = self_check()
        return 0 if passed else 1

    if args.samples > 0:
        random_sample_sweep(
            N=args.N,
            W=args.W,
            a=args.a,
            q=args.q,
            n_samples=args.samples,
            topk=args.topk,
            arc_width=args.arc_width,
            seed=args.seed,
            verbose=True,
        )
    else:
        run_experiment(
            N=args.N,
            W=args.W,
            a=args.a,
            q=args.q,
            beta=args.beta,
            weighted=args.weighted,
            topk=args.topk,
            arc_width=args.arc_width,
            verbose=True,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
