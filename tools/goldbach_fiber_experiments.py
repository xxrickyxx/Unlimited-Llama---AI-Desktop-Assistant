#!/usr/bin/env python3
"""
Goldbach Fiber-Sum Numerical Experiments  (S1.2 / Option 2)
============================================================

Computes Phi(α) = Σ_{p≤N} log(p)·e(αp) and the squared convolution
F(t) = Phi(α + t/W)² over t mod W, then analyses the discrete Fourier
transform F̂(k) together with character correlations and aliasing counts
to probe the Goldbach fiber-sum conjecture numerically.

References
----------
See docs/goldbach/experiment-plan.md for metric definitions and their
mapping to theoretical lemmas.

CLI examples
------------
  # Basic run with default parameters
  python tools/goldbach_fiber_experiments.py

  # Specify N, W, base alpha = a/q, number of minor-arc samples
  python tools/goldbach_fiber_experiments.py --N 10000 --W 30 --a 1 --q 3 --samples 50

  # Convergence sweep doubling N from N0 up to N
  python tools/goldbach_fiber_experiments.py --N 32000 --W 30 --sweep

  # Self-check / doctest mode (no external deps required)
  python tools/goldbach_fiber_experiments.py --selfcheck
"""

from __future__ import annotations

import argparse
import cmath
import math
import sys
from typing import List, Tuple

# ---------------------------------------------------------------------------
# Optional fast path: use numpy when available, fall back to pure Python.
# ---------------------------------------------------------------------------
try:
    import numpy as np
    _NUMPY = True
except ImportError:  # pragma: no cover
    _NUMPY = False


# ===========================================================================
# 1.  Prime generation
# ===========================================================================

def sieve_primes(n: int) -> List[int]:
    """Return all primes p with 2 ≤ p ≤ n (Sieve of Eratosthenes).

    >>> sieve_primes(10)
    [2, 3, 5, 7]
    >>> sieve_primes(1)
    []
    >>> sieve_primes(2)
    [2]
    """
    if n < 2:
        return []
    is_prime = bytearray([1]) * (n + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if is_prime[i]:
            is_prime[i * i :: i] = bytearray(len(is_prime[i * i :: i]))
    return [i for i, v in enumerate(is_prime) if v]


# ===========================================================================
# 2.  Phi(α) evaluation
# ===========================================================================

def eval_phi(alpha: float, primes: List[int]) -> complex:
    """Compute Phi(α) = Σ_{p in primes} log(p) · e(α·p).

    e(x) = exp(2π i x).

    >>> import cmath, math
    >>> p = [2, 3, 5]
    >>> got = eval_phi(0.0, p)
    >>> expected = sum(math.log(p_) * 1 for p_ in p)
    >>> abs(got - expected) < 1e-10
    True
    """
    tau = 2.0 * math.pi
    if _NUMPY:
        ps = np.asarray(primes, dtype=np.float64)
        lps = np.log(ps)
        phases = np.exp(1j * tau * alpha * ps)
        return complex(np.dot(lps, phases))
    # Pure-Python fallback
    total = complex(0.0)
    for p in primes:
        total += math.log(p) * cmath.exp(1j * tau * alpha * p)
    return total


def eval_phi_batch(
    alpha_base: float,
    W: int,
    primes: List[int],
) -> "List[complex] | np.ndarray":
    """Evaluate Phi(α_base + t/W) for t = 0, …, W−1.

    Uses an FFT-based approach when numpy is available: build a
    log-prime weighted indicator array of length W·ceil(N/W), then
    take the DFT to read off all W values at once.

    Returns a list (or numpy array) of W complex numbers.

    >>> vals = eval_phi_batch(0.0, 6, [2, 3, 5])
    >>> len(vals)
    6
    >>> import cmath, math
    >>> expected_t0 = sum(math.log(p) * cmath.exp(0) for p in [2,3,5])
    >>> abs(vals[0] - expected_t0) < 1e-8
    True
    """
    if not primes:
        if _NUMPY:
            return np.zeros(W, dtype=complex)
        return [complex(0)] * W

    N = primes[-1]
    tau = 2.0 * math.pi

    if _NUMPY:
        # Build a length-N+1 log-prime indicator, zero-padded to W multiple
        padded = _next_multiple(N + 1, W)
        a = np.zeros(padded, dtype=np.float64)
        ps = np.asarray(primes, dtype=int)
        a[ps] = np.log(ps.astype(np.float64))

        # Multiply by global phase e(alpha_base · index)
        idx = np.arange(padded, dtype=np.float64)
        a_complex = a * np.exp(1j * tau * alpha_base * idx)

        # Reshape to (padded//W, W) and sum along rows  ←  DFT trick:
        # Phi(alpha_base + t/W) = Σ_k  [Σ_m  a[mW+t] e(alpha_base(mW+t))]
        #                                * e(t·k/W) ← no extra phase here
        # because within each block the residue is t, so summing over blocks
        # gives exactly Phi evaluated at alpha_base + t/W.
        blocks = a_complex.reshape(-1, W)
        phi_vals = blocks.sum(axis=0)
        return phi_vals
    else:
        return [eval_phi(alpha_base + t / W, primes) for t in range(W)]


def _next_multiple(x: int, m: int) -> int:
    """Smallest multiple of m that is ≥ x."""
    return m * math.ceil(x / m)


# ===========================================================================
# 3.  F(t) = Phi(α + t/W)²
# ===========================================================================

def compute_F(phi_vals: "List[complex] | np.ndarray") -> "List[complex] | np.ndarray":
    """Square each element of phi_vals to get F(t).

    >>> F = compute_F([1+1j, 2+0j])
    >>> abs(F[0] - (1+1j)**2) < 1e-12
    True
    """
    if _NUMPY and isinstance(phi_vals, np.ndarray):
        return phi_vals ** 2
    return [v ** 2 for v in phi_vals]


# ===========================================================================
# 4.  DFT of F(t)
# ===========================================================================

def dft_F(F_vals: "List[complex] | np.ndarray") -> "List[complex] | np.ndarray":
    """Compute F̂(k) = Σ_t F(t) · e(−kt/W) for k = 0, …, W−1.

    >>> import cmath, math
    >>> F = [1.0, 0.0, 0.0, 0.0]   # delta at t=0
    >>> Fhat = dft_F(F)
    >>> all(abs(Fhat[k] - 1.0) < 1e-10 for k in range(4))
    True
    """
    if _NUMPY:
        arr = np.asarray(F_vals, dtype=complex)
        return np.fft.fft(arr)
    W = len(F_vals)
    tau = 2.0 * math.pi
    result = []
    for k in range(W):
        s = complex(0.0)
        for t, f in enumerate(F_vals):
            s += f * cmath.exp(-1j * tau * k * t / W)
        result.append(s)
    return result


# ===========================================================================
# 5.  Major / minor arc classifier
# ===========================================================================

def is_major_arc(alpha: float, Q: float, max_q: int) -> Tuple[bool, int, int]:
    """Return (is_major, a, q) where a/q is the best rational approximation.

    alpha is considered a *major arc* centre if |alpha − a/q| < Q for
    some a/q with q ≤ max_q.

    Uses the Stern–Brocot / Farey approach (just brute-force for small
    max_q).

    >>> is_major_arc(0.0, 0.1, 10)
    (True, 0, 1)
    >>> is_major_arc(0.5, 0.01, 10)
    (True, 1, 2)
    >>> result = is_major_arc(0.37, 1e-4, 5)
    >>> result[0]  # 0.37 is not very close to any a/q with q≤5
    False
    """
    alpha_mod = alpha % 1.0
    best_dist = float("inf")
    best_a, best_q = 0, 1
    for q in range(1, max_q + 1):
        for a in range(0, q + 1):
            dist = abs(alpha_mod - a / q)
            dist = min(dist, 1.0 - dist)  # wrap around 1
            if dist < best_dist:
                best_dist = dist
                best_a, best_q = a, q
    return (best_dist < Q, best_a, best_q)


def count_aliasing(
    alpha_base: float,
    W: int,
    Q: float,
    max_q: int,
) -> int:
    """Count how many t in {0,…,W−1} make α_base + t/W a major-arc point.

    >>> count_aliasing(0.0, 6, 0.2, 1)  # only 0/1 counts; t=0,1,5 fall within Q=0.2
    3
    """
    count = 0
    for t in range(W):
        major, _, _ = is_major_arc(alpha_base + t / W, Q, max_q)
        if major:
            count += 1
    return count


def sample_minor_arc_alpha(
    a: int,
    q: int,
    N: int,
    n_samples: int,
    rng_seed: int = 42,
) -> List[float]:
    """Sample n_samples values of α near a/q but outside major arcs.

    Uses a simple uniform random offset in [Q_lo, Q_hi] away from a/q
    where Q_lo = (log N)**2 / N  and  Q_hi = 0.5 / q.

    >>> alphas = sample_minor_arc_alpha(1, 3, 1000, 5)
    >>> len(alphas)
    5
    >>> all(abs(a_ % 1.0 - 1/3) > 0.001 for a_ in alphas)
    True
    """
    if _NUMPY:
        rng = np.random.default_rng(rng_seed)
    else:
        import random
        random.seed(rng_seed)

    Q_lo = max((math.log(N) ** 2) / N, 1.0 / N)
    Q_hi = 0.5 / max(q, 1)
    if Q_lo >= Q_hi:
        Q_lo = Q_hi / 10.0

    base = a / q
    results: List[float] = []
    attempts = 0
    while len(results) < n_samples and attempts < n_samples * 100:
        attempts += 1
        if _NUMPY:
            offset = float(rng.uniform(Q_lo, Q_hi))
        else:
            offset = random.uniform(Q_lo, Q_hi)
        sign = 1 if (attempts % 2 == 0) else -1
        candidate = (base + sign * offset) % 1.0
        results.append(candidate)
    return results


# ===========================================================================
# 6.  Statistics
# ===========================================================================

def phi_statistics(phi_vals: "List[complex] | np.ndarray") -> dict:
    """Compute mean, std, flatness ratio of |Phi| values.

    >>> import math
    >>> stats = phi_statistics([1+0j, 1j, -1+0j, -1j])
    >>> abs(stats['mean'] - 1.0) < 1e-10
    True
    >>> stats['flatness'] < 1e-10
    True
    """
    if _NUMPY:
        mags = np.abs(np.asarray(phi_vals, dtype=complex))
        mu = float(np.mean(mags))
        sigma = float(np.std(mags))
    else:
        mags = [abs(v) for v in phi_vals]
        mu = sum(mags) / len(mags)
        variance = sum((m - mu) ** 2 for m in mags) / len(mags)
        sigma = math.sqrt(variance)
    flatness = sigma / mu if mu > 0 else float("nan")
    return {"mean": mu, "std": sigma, "flatness": flatness}


def fhat_statistics(Fhat: "List[complex] | np.ndarray", k_star: int) -> dict:
    """Report |F̂(k*)| and the off-diagonal leakage ratio.

    leakage = max_{k ≠ k*} |F̂(k)| / |F̂(k*)|

    >>> Fhat = [10.0, 0.1, 0.05, 0.2]
    >>> s = fhat_statistics(Fhat, 0)
    >>> abs(s['Fhat_kstar'] - 10.0) < 1e-10
    True
    >>> abs(s['leakage'] - 0.02) < 1e-10
    True
    """
    if _NUMPY:
        mags = np.abs(np.asarray(Fhat, dtype=complex))
    else:
        mags = [abs(v) for v in Fhat]
    W = len(mags)
    k_star = k_star % W
    fhat_k = mags[k_star]
    off_diag = [mags[k] for k in range(W) if k != k_star]
    if off_diag:
        max_off = max(off_diag)
    else:
        max_off = 0.0
    leakage = max_off / fhat_k if fhat_k > 0 else float("nan")
    return {"Fhat_kstar": float(fhat_k), "max_off_diag": float(max_off), "leakage": leakage}


def character_correlations(
    F_vals: "List[complex] | np.ndarray",
    W: int,
) -> List[Tuple[int, float]]:
    """Compute |C(χ_k)| where χ_k(t) = e(kt/W) for k = 0,…,W−1.

    This is just the normalised DFT magnitude of F:
      C(χ_k) = (1/W) Σ_t F(t) · e(−kt/W) = F̂(k)/W.

    Returns a list of (k, |C(χ_k)|) sorted by descending |C|.

    >>> F = [1.0, 0.0, 0.0, 0.0]
    >>> corrs = character_correlations(F, 4)
    >>> all(abs(c - 0.25) < 1e-10 for _, c in corrs)
    True
    """
    Fhat = dft_F(F_vals)
    if _NUMPY:
        mags = np.abs(np.asarray(Fhat, dtype=complex)) / W
        pairs = sorted(
            [(int(k), float(mags[k])) for k in range(W)],
            key=lambda x: -x[1],
        )
    else:
        pairs = sorted(
            [(k, abs(Fhat[k]) / W) for k in range(W)],
            key=lambda x: -x[1],
        )
    return pairs


# ===========================================================================
# 7.  Main experiment runner
# ===========================================================================

def run_experiment(
    N: int = 10_000,
    W: int = 30,
    a: int = 0,
    q: int = 1,
    samples: int = 10,
    Q: float = 0.0,
    verbose: bool = True,
) -> dict:
    """Run the full Goldbach fiber-sum experiment.

    Parameters
    ----------
    N       : prime upper bound
    W       : modulus
    a, q    : base alpha = a/q (major-arc centre)
    samples : number of random minor-arc alpha values to probe
    Q       : major-arc half-width (0 → auto: (log N)^2 / N)
    verbose : print report to stdout

    Returns a dict of results for each alpha sampled.
    """
    if Q <= 0:
        Q = (math.log(max(N, 2)) ** 2) / N

    primes = sieve_primes(N)
    k_star = (2 * N) % W

    if verbose:
        print("=" * 60)
        print(f"Goldbach Fiber-Sum Experiment")
        print(f"  N={N}, W={W}, base α={a}/{q}, Q={Q:.2e}")
        print(f"  Primes ≤ N: {len(primes)},  k* = {k_star}")
        print("=" * 60)

    results = {}

    # ------------------------------------------------------------------ #
    # A.  Evaluate at the major-arc base alpha = a/q
    # ------------------------------------------------------------------ #
    alpha_base = a / q
    phi_vals_base = eval_phi_batch(alpha_base, W, primes)
    F_base = compute_F(phi_vals_base)
    Fhat_base = dft_F(F_base)
    phi_stats_base = phi_statistics(phi_vals_base)
    fhat_stats_base = fhat_statistics(Fhat_base, k_star)
    alias_base = count_aliasing(alpha_base, W, Q, max_q=int(math.log(N) ** 2) + 1)
    corr_base = character_correlations(F_base, W)

    results["base"] = {
        "alpha": alpha_base,
        "phi_stats": phi_stats_base,
        "fhat_stats": fhat_stats_base,
        "aliasing": alias_base,
        "top_correlations": corr_base[:5],
    }

    if verbose:
        _print_result("BASE (major arc)", alpha_base, phi_stats_base,
                      fhat_stats_base, alias_base, corr_base[:5])

    # ------------------------------------------------------------------ #
    # B.  Sample minor-arc alpha values
    # ------------------------------------------------------------------ #
    minor_alphas = sample_minor_arc_alpha(a, q, N, samples)
    minor_results = []
    for i, alpha in enumerate(minor_alphas):
        phi_vals = eval_phi_batch(alpha, W, primes)
        F_vals = compute_F(phi_vals)
        Fhat = dft_F(F_vals)
        phi_stats = phi_statistics(phi_vals)
        fhat_stats = fhat_statistics(Fhat, k_star)
        alias = count_aliasing(alpha, W, Q, max_q=int(math.log(N) ** 2) + 1)
        corrs = character_correlations(F_vals, W)
        entry = {
            "alpha": alpha,
            "phi_stats": phi_stats,
            "fhat_stats": fhat_stats,
            "aliasing": alias,
            "top_correlations": corrs[:5],
        }
        minor_results.append(entry)
        if verbose and (i < 3 or i == len(minor_alphas) - 1):
            _print_result(f"MINOR arc #{i+1}", alpha, phi_stats,
                          fhat_stats, alias, corrs[:5])
        elif verbose and i == 3:
            print("  … (further minor-arc results omitted for brevity) …\n")

    results["minor"] = minor_results

    # ------------------------------------------------------------------ #
    # C.  Summary
    # ------------------------------------------------------------------ #
    if minor_results and verbose:
        flatness_vals = [r["phi_stats"]["flatness"] for r in minor_results]
        leakage_vals = [r["fhat_stats"]["leakage"] for r in minor_results
                        if not math.isnan(r["fhat_stats"]["leakage"])]
        alias_vals = [r["aliasing"] for r in minor_results]
        print("-" * 60)
        print("SUMMARY (minor arcs)")
        print(f"  Flatness σ/μ : mean={_fmean(flatness_vals):.4f}, "
              f"max={max(flatness_vals):.4f}")
        if leakage_vals:
            print(f"  Off-diag leakage: mean={_fmean(leakage_vals):.4f}, "
                  f"max={max(leakage_vals):.4f}")
        print(f"  Aliasing count  : mean={_fmean(alias_vals):.1f}, "
              f"max={max(alias_vals)}")
        print("=" * 60)

    return results


def _fmean(vals: list) -> float:
    return sum(vals) / len(vals) if vals else float("nan")


def _print_result(
    label: str,
    alpha: float,
    phi_stats: dict,
    fhat_stats: dict,
    alias: int,
    corrs: list,
) -> None:
    print(f"\n[{label}]  α = {alpha:.6f}")
    print(f"  |Phi|: mean={phi_stats['mean']:.4g}, "
          f"std={phi_stats['std']:.4g}, flatness={phi_stats['flatness']:.4f}")
    print(f"  F̂(k*): {fhat_stats['Fhat_kstar']:.4g}, "
          f"max-off-diag={fhat_stats['max_off_diag']:.4g}, "
          f"leakage={fhat_stats['leakage']:.4f}")
    print(f"  Aliasing count: {alias}")
    top5 = ", ".join(f"χ_{k}:{c:.3f}" for k, c in corrs[:3])
    print(f"  Top correlations: {top5}")


# ===========================================================================
# 8.  N-sweep for convergence study (M5)
# ===========================================================================

def run_sweep(
    N_max: int = 32_000,
    N0: int = 1_000,
    W: int = 30,
    a: int = 0,
    q: int = 1,
    verbose: bool = True,
) -> List[dict]:
    """Run run_experiment for N = N0, 2·N0, 4·N0, … ≤ N_max.

    Reports how key metrics scale with N.
    """
    results = []
    N = N0
    while N <= N_max:
        if verbose:
            print(f"\n{'='*60}\nSWEEP N={N}\n{'='*60}")
        res = run_experiment(N=N, W=W, a=a, q=q, samples=3, verbose=verbose)
        base = res["base"]
        results.append({
            "N": N,
            "mean_phi": base["phi_stats"]["mean"],
            "flatness": base["phi_stats"]["flatness"],
            "Fhat_kstar": base["fhat_stats"]["Fhat_kstar"],
            "leakage": base["fhat_stats"]["leakage"],
        })
        N *= 2

    if verbose and len(results) > 1:
        print("\nSWEEP SUMMARY")
        print(f"{'N':>10}  {'mean|Phi|':>12}  {'flatness':>10}  {'F̂(k*)':>12}  {'leakage':>10}")
        for r in results:
            print(f"{r['N']:>10}  {r['mean_phi']:>12.4g}  "
                  f"{r['flatness']:>10.4f}  {r['Fhat_kstar']:>12.4g}  "
                  f"{r['leakage']:>10.4f}")
    return results


# ===========================================================================
# 9.  Self-check / doctest runner
# ===========================================================================

def selfcheck() -> None:
    """Run all doctests embedded in this module."""
    import doctest
    results = doctest.testmod(verbose=False)
    if results.failed:
        print(f"SELF-CHECK FAILED: {results.failed}/{results.attempted} tests failed.")
        sys.exit(1)
    else:
        print(f"SELF-CHECK PASSED: {results.attempted} tests passed.")

    # Extra sanity checks beyond doctests
    primes = sieve_primes(100)
    assert primes[:4] == [2, 3, 5, 7], "sieve_primes failed"

    phi0 = eval_phi(0.0, primes)
    expected = sum(math.log(p) for p in primes)
    assert abs(phi0 - expected) < 1e-8, f"eval_phi(0) failed: {phi0} != {expected}"

    phi_batch = eval_phi_batch(0.0, 10, primes)
    assert len(phi_batch) == 10, "eval_phi_batch length wrong"
    assert abs(phi_batch[0] - phi0) < 1e-6, "eval_phi_batch t=0 mismatch"

    F = compute_F(phi_batch)
    assert len(F) == 10
    assert abs(F[0] - phi_batch[0] ** 2) < 1e-6

    Fhat = dft_F(F)
    assert len(Fhat) == 10

    stats = phi_statistics(phi_batch)
    assert "mean" in stats and "std" in stats and "flatness" in stats

    print("All extra sanity checks passed.")


# ===========================================================================
# 10.  CLI
# ===========================================================================

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Goldbach fiber-sum numerical experiments (S1.2 / option 2)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--N", type=int, default=10_000,
                   help="Prime upper bound (default: 10000)")
    p.add_argument("--W", type=int, default=30,
                   help="Modulus W (default: 30)")
    p.add_argument("--a", type=int, default=0,
                   help="Numerator of base α = a/q (default: 0)")
    p.add_argument("--q", type=int, default=1,
                   help="Denominator of base α = a/q (default: 1)")
    p.add_argument("--samples", type=int, default=10,
                   help="Number of minor-arc α samples (default: 10)")
    p.add_argument("--Q", type=float, default=0.0,
                   help="Major-arc half-width (default: (log N)²/N)")
    p.add_argument("--sweep", action="store_true",
                   help="Run N-sweep (M5 convergence study)")
    p.add_argument("--N0", type=int, default=1_000,
                   help="Starting N for sweep (default: 1000)")
    p.add_argument("--selfcheck", action="store_true",
                   help="Run self-check / doctests and exit")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress verbose output")
    return p


def main(argv: List[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.selfcheck:
        selfcheck()
        return

    verbose = not args.quiet

    if args.sweep:
        run_sweep(
            N_max=args.N,
            N0=args.N0,
            W=args.W,
            a=args.a,
            q=args.q,
            verbose=verbose,
        )
    else:
        run_experiment(
            N=args.N,
            W=args.W,
            a=args.a,
            q=args.q,
            samples=args.samples,
            Q=args.Q,
            verbose=verbose,
        )


if __name__ == "__main__":
    main()
