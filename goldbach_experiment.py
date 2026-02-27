#!/usr/bin/env python3
"""goldbach_experiment.py — Goldbach Fiber-Sum Experiment
==========================================================

Implements the Hardy–Littlewood circle-method experiment for probing
Goldbach representations via exponential sums over primes.

Default parameters (all overridable via CLI):
  N = 1_000_000  (even upper bound for primes)
  W = 1024       (window / frequency-bin count)
  A = 6          (major-arc exponent; P = (ln N)^A)

Usage
-----
    python goldbach_experiment.py [--N N] [--W W] [--A A] [options]

Run with --help for the full option list.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys

import numpy as np


# ---------------------------------------------------------------------------
# Prime sieve
# ---------------------------------------------------------------------------

def sieve_primes(n: int) -> np.ndarray:
    """Return all primes <= *n* as an int64 numpy array (Sieve of Eratosthenes)."""
    if n < 2:
        return np.array([], dtype=np.int64)
    sieve = np.ones(n + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, math.isqrt(n) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.where(sieve)[0].astype(np.int64)


# ---------------------------------------------------------------------------
# Major-arc detection
# ---------------------------------------------------------------------------

def major_arc_threshold(N: int, A: float) -> float:
    """Return P = (ln N)^A, the major-arc denominator threshold."""
    return math.log(N) ** A


def _major_arc_flags(alphas: np.ndarray, P: float, N: int) -> np.ndarray:
    """Return a boolean array: True where each alpha lies in a major arc.

    Definition: alpha in a major arc iff there exist integers q with
    1 <= q <= P and a such that |alpha - a/q| <= P / (N * q).

    The minimum of |alpha - a/q| over all integers a equals ||q*alpha||/q,
    where ||x|| = dist(x, Z) is the distance to the nearest integer.  This
    allows the check to be fully vectorised over alphas.  Fractions with
    gcd(a, q) > 1 are subsumed by their reduced form (smaller q), so the
    gcd condition need not be checked explicitly.
    """
    alphas = np.asarray(alphas, dtype=np.float64) % 1.0
    result = np.zeros(len(alphas), dtype=bool)
    q_max = int(min(math.floor(P), N))

    for q in range(1, q_max + 1):
        radius = P / (N * q)
        if radius >= 0.5:
            # This q's arc already covers the whole circle; so does every
            # larger q (radius grows with q).  Return immediately.
            result[:] = True
            return result

        diff = np.abs(alphas - np.round(alphas * q) / q)
        result |= diff <= radius

        if result.all():
            return result

    return result


def is_major_arc(alpha: float, P: float, N: int) -> bool:
    """Return True if *alpha* lies in a major arc."""
    return bool(_major_arc_flags(np.array([alpha]), P, N)[0])


def aliasing_count(alpha: float, W: int, P: float, N: int) -> int:
    """B(alpha) = #{t in {0,...,W-1} : alpha + t/W (mod 1) in a major arc}."""
    alphas = (alpha + np.arange(W, dtype=np.float64) / W) % 1.0
    return int(_major_arc_flags(alphas, P, N).sum())


# ---------------------------------------------------------------------------
# Minor-arc sampling
# ---------------------------------------------------------------------------

def sample_minor_arc_alpha(
    P: float, N: int, rng: np.random.Generator, max_attempts: int = 100_000
) -> float:
    """Draw a uniform random alpha in [0,1) from the minor arcs.

    Uses rejection sampling: draws uniformly from [0,1) and discards any
    point that lies in a major arc.

    Raises RuntimeError if no minor-arc point is found within *max_attempts*
    (typically indicates P >= N, i.e. all of [0,1) is covered by major arcs).
    """
    for _ in range(max_attempts):
        alpha = rng.uniform(0.0, 1.0)
        if not is_major_arc(alpha, P, N):
            return alpha
    raise RuntimeError(
        f"No minor-arc alpha found in {max_attempts} attempts "
        f"(P={P:.4g}, N={N}). "
        "This usually means P >= N — reduce A or increase N."
    )


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_phi_array(alpha: float, W: int, primes: np.ndarray) -> np.ndarray:
    """Return phi_vals[t] = Phi(alpha + t/W) for t = 0,...,W-1.

    Derivation
    ----------
    Phi(alpha + t/W) = sum_{p<=N} log(p) * exp(2*pi*i*(alpha + t/W)*p)
                     = sum_p  c_p * exp(2*pi*i * t * (p mod W) / W)

    where  c_p = log(p) * exp(2*pi*i*alpha*p).

    Collecting by residue class modulo W:

        g[k] = sum_{p <= N, p ≡ k (mod W)} c_p

    gives  Phi(alpha + t/W) = sum_k g[k] * exp(2*pi*i*k*t/W)
                            = W * IFFT(g)[t].

    Cost: O(pi(N) + W log W), independent of W*pi(N).

    Returns a complex128 array of length W.
    """
    log_p = np.log(primes.astype(np.float64))
    base_phases = 2.0 * math.pi * alpha * primes.astype(np.float64)
    c = log_p * np.exp(1j * base_phases)

    k_indices = primes % W
    g = np.zeros(W, dtype=np.complex128)
    np.add.at(g, k_indices, c)

    # Phi(alpha + t/W) = sum_k g[k] * exp(2*pi*i*k*t/W) = W * IFFT(g)[t]
    return np.fft.ifft(g) * W


# ---------------------------------------------------------------------------
# Experiment runner
# ---------------------------------------------------------------------------

def run_experiment(
    N: int = 1_000_000,
    W: int = 1024,
    A: float = 6.0,
    top_k: int = 10,
    seed: int = 42,
    save_csv: str | None = None,
    verbose: bool = True,
) -> dict:
    """Run the full Goldbach fiber-sum experiment and return a results dict.

    Steps
    -----
    1. Sieve all primes p <= N.
    2. Compute P = (ln N)^A.
    3. Sample a minor-arc alpha via rejection sampling.
    4. Compute Phi(alpha + t/W) for t = 0,...,W-1 via the g-vector/IFFT method.
    5. Form  F[t] = Phi(alpha + t/W)^2.
    6. Compute Fhat = FFT(F).
    7. Report |Fhat[k]| for the top-*top_k* indices and for k = N mod W.
    8. Compute B(alpha) = aliasing count.
    9. Optionally write all |Fhat[k]| to a CSV file.
    """
    rng = np.random.default_rng(seed)

    # 1. Primes
    if verbose:
        print(f"[1/7] Sieving primes up to N = {N:,} ...", flush=True)
    primes = sieve_primes(N)
    if verbose:
        print(f"      pi({N:,}) = {len(primes):,} primes")

    # 2. P
    P = major_arc_threshold(N, A)
    if verbose:
        print(f"[2/7] P = (ln {N})^{A} = {P:.6g}")
    if P >= N:
        print(
            f"      WARNING: P = {P:.4g} >= N = {N}. "
            "Every alpha lies in a major arc; minor-arc sampling will fall back "
            "to a random alpha. Use --A 2 (or smaller) for genuine minor-arc "
            "experiments with this N.",
            file=sys.stderr,
        )

    # 3. Sample minor-arc alpha
    if verbose:
        print("[3/7] Sampling minor-arc alpha ...", flush=True)
    fallback_used = False
    try:
        alpha = sample_minor_arc_alpha(P, N, rng)
    except RuntimeError as exc:
        print(f"      WARNING: {exc}", file=sys.stderr)
        alpha = rng.uniform(0.0, 1.0)
        fallback_used = True
    if verbose:
        label = " (fallback — not a genuine minor arc)" if fallback_used else ""
        print(f"      alpha = {alpha:.10f}{label}")

    # 4. Compute Phi(alpha + t/W) for all t
    if verbose:
        print("[4/7] Computing Phi(alpha + t/W) for t = 0 ... W-1 ...", flush=True)
    phi_vals = compute_phi_array(alpha, W, primes)

    # 5. F[t] = Phi(alpha + t/W)^2
    if verbose:
        print("[5/7] Forming F[t] = Phi(alpha + t/W)^2 ...", flush=True)
    F = phi_vals ** 2

    # 6. FFT
    if verbose:
        print("[6/7] Computing Fhat = FFT(F) ...", flush=True)
    Fhat = np.fft.fft(F)
    Fhat_abs = np.abs(Fhat)

    # 7. Aliasing count
    if verbose:
        print("[7/7] Computing aliasing count B(alpha) ...", flush=True)
    B = aliasing_count(alpha, W, P, N)

    k_N = int(N % W)
    results = {
        "N": N,
        "W": W,
        "A": A,
        "P": P,
        "alpha": alpha,
        "B": B,
        "fallback_used": fallback_used,
        "primes_count": len(primes),
        "phi_vals": phi_vals,
        "F": F,
        "Fhat": Fhat,
        "Fhat_abs": Fhat_abs,
        "k_N": k_N,
    }

    if verbose:
        _print_results(results, top_k)

    if save_csv:
        _write_csv(Fhat, Fhat_abs, W, save_csv, verbose)

    return results


def _print_results(r: dict, top_k: int) -> None:
    W = r["W"]
    Fhat = r["Fhat"]
    Fhat_abs = r["Fhat_abs"]
    k_N = r["k_N"]

    print("\n" + "=" * 60)
    print("GOLDBACH FIBER-SUM EXPERIMENT — RESULTS")
    print("=" * 60)
    print(f"  N         = {r['N']:,}")
    print(f"  W         = {r['W']}")
    print(f"  A         = {r['A']}")
    print(f"  P         = {r['P']:.6g}")
    print(f"  alpha     = {r['alpha']:.10f}")
    print(f"  B(alpha)  = {r['B']} / {W}  ({100 * r['B'] / W:.1f}% major-arc shifts)")
    print(f"  k_N = N mod W = {k_N}")
    print(f"  |Fhat[k_N]| = {Fhat_abs[k_N]:.6e}")
    print()

    top_indices = np.argsort(Fhat_abs)[::-1][:top_k]
    print(f"  Top-{top_k} |Fhat[k]| entries:")
    print(f"  {'k':>6}  {'|Fhat[k]|':>14}  {'phase (deg)':>12}")
    print(f"  {'─' * 6}  {'─' * 14}  {'─' * 12}")
    for idx in top_indices:
        phase_deg = math.degrees(math.atan2(Fhat[idx].imag, Fhat[idx].real))
        print(f"  {idx:>6}  {Fhat_abs[idx]:>14.6e}  {phase_deg:>12.2f}")

    print()
    print("  Summary statistics on |Fhat|:")
    print(f"    mean = {Fhat_abs.mean():.6e}")
    print(f"    std  = {Fhat_abs.std():.6e}")
    print(f"    max  = {Fhat_abs.max():.6e}  at k = {int(Fhat_abs.argmax())}")
    print(f"    min  = {Fhat_abs.min():.6e}  at k = {int(Fhat_abs.argmin())}")
    print("=" * 60)


def _write_csv(
    Fhat: np.ndarray, Fhat_abs: np.ndarray, W: int, path: str, verbose: bool
) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["k", "Fhat_real", "Fhat_imag", "Fhat_abs"])
        for k in range(W):
            writer.writerow([k, Fhat[k].real, Fhat[k].imag, Fhat_abs[k]])
    if verbose:
        print(f"\n  CSV written to: {path}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Goldbach fiber-sum experiment (Hardy-Littlewood circle method).\n\n"
            "Runs automatically with defaults N=1_000_000, W=1024, A=6.\n"
            "All parameters are overridable via CLI flags.\n\n"
            "NOTE: With default A=6 and N=1_000_000, P=(ln N)^A > N, so every\n"
            "alpha is in a major arc and a random fallback alpha is used.\n"
            "Use --A 2 for genuine minor-arc sampling at this N."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--N", type=int, default=1_000_000,
        metavar="N", help="Upper prime bound (must be even) [default: 1,000,000]",
    )
    p.add_argument(
        "--W", type=int, default=1024,
        metavar="W", help="Window / frequency-bin count (power of 2 recommended) [default: 1024]",
    )
    p.add_argument(
        "--A", type=float, default=6.0,
        metavar="A", help="Major-arc exponent: P = (ln N)^A [default: 6.0]",
    )
    p.add_argument(
        "--top-k", type=int, default=10,
        metavar="K", help="Number of top |Fhat| entries to display [default: 10]",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for alpha sampling [default: 42]",
    )
    p.add_argument(
        "--save-csv", type=str, default=None,
        metavar="PATH", help="Save full Fhat spectrum to PATH (CSV)",
    )
    p.add_argument(
        "--quiet", action="store_true",
        help="Suppress progress messages (results are still printed)",
    )
    return p


def main(argv=None):
    args = _build_parser().parse_args(argv)

    if args.N % 2 != 0:
        print("Error: N must be even.", file=sys.stderr)
        sys.exit(1)
    if args.W < 1:
        print("Error: W must be a positive integer.", file=sys.stderr)
        sys.exit(1)
    if args.A <= 0:
        print("Error: A must be positive.", file=sys.stderr)
        sys.exit(1)

    run_experiment(
        N=args.N,
        W=args.W,
        A=args.A,
        top_k=args.top_k,
        seed=args.seed,
        save_csv=args.save_csv,
        verbose=not args.quiet,
    )


if __name__ == "__main__":
    main()
