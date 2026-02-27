#!/usr/bin/env python3
"""Goldbach fiber-sum experiment tool.

Conventions
-----------
* N  : target even integer.
* Phi(alpha) = sum_{p <= N} log(p) * exp(2 pi i alpha p)
* Major arcs : P = (log N)^A.  alpha is *major* if there exists q <= P and
  a with gcd(a,q)=1 such that |alpha - a/q| <= P / (N * q).
* Fiber DFT : F_alpha(t) = Phi(alpha + t/W)^2  for  t mod W.
  DFT coefficient for k* = N mod W :
      F_hat(k*) = sum_{t mod W} F_alpha(t) * exp(-2 pi i k* t / W)
"""

import argparse
import csv
import math
import sys
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Sieve
# ---------------------------------------------------------------------------

def sieve_primes(n: int) -> np.ndarray:
    """Return array of all primes <= n using the Sieve of Eratosthenes."""
    if n < 2:
        return np.array([], dtype=np.int64)
    composite = np.zeros(n + 1, dtype=bool)
    composite[0] = composite[1] = True
    for i in range(2, int(n**0.5) + 1):
        if not composite[i]:
            composite[i * i :: i] = True
    return np.where(~composite)[0].astype(np.int64)


# ---------------------------------------------------------------------------
# Major-arc test
# ---------------------------------------------------------------------------

def is_major_arc(alpha: float, n: int, A: float) -> bool:
    """Return True if *alpha* lies on a major arc for the given N and A.

    A real alpha in [0,1) is *major* if there exist integers q, a with
    1 <= q <= P, gcd(a,q)=1, and |alpha - a/q| <= P / (N * q), where
    P = (log N)^A.
    """
    if n < 2:
        return False
    P = math.log(n) ** A
    max_q = int(P)
    for q in range(1, max_q + 1):
        for a in range(0, q):
            if math.gcd(a, q) != 1:
                continue
            frac = a / q
            dist = abs(alpha - frac)
            # wrap distance on the circle [0,1)
            dist = min(dist, 1.0 - dist)
            if dist <= P / (n * q):
                return True
    return False


def is_major_arc_shift(alpha_shifted: float, n: int, A: float) -> bool:
    """Wrap alpha_shifted into [0,1) then apply the major-arc test."""
    alpha_shifted = alpha_shifted % 1.0
    return is_major_arc(alpha_shifted, n, A)


# ---------------------------------------------------------------------------
# Phi computation (vectorised)
# ---------------------------------------------------------------------------

def compute_phi_vector(alphas: np.ndarray, primes: np.ndarray) -> np.ndarray:
    """Compute Phi(alpha) for each alpha in *alphas*.

    Phi(alpha) = sum_{p <= N} log(p) * exp(2 pi i alpha p)

    Parameters
    ----------
    alphas : shape (M,)   – evaluation points
    primes : shape (K,)   – primes up to N

    Returns
    -------
    phi : complex array of shape (M,)
    """
    log_p = np.log(primes.astype(np.float64))  # shape (K,)
    # phases[m, k] = 2 pi i * alphas[m] * primes[k]
    phases = 2.0 * math.pi * np.outer(alphas, primes)  # (M, K)
    weighted = log_p[np.newaxis, :] * np.exp(1j * phases)  # (M, K)
    return weighted.sum(axis=1)  # (M,)


def compute_fiber(alpha: float, W: int, primes: np.ndarray) -> np.ndarray:
    """Compute the fiber vector F_alpha(t) = Phi(alpha + t/W)^2 for t=0..W-1."""
    t_vals = np.arange(W, dtype=np.float64)
    shift_alphas = alpha + t_vals / W  # shape (W,)
    phi_vec = compute_phi_vector(shift_alphas, primes)  # shape (W,)
    return phi_vec ** 2  # F(t) = Phi(alpha+t/W)^2


# ---------------------------------------------------------------------------
# Aliasing count
# ---------------------------------------------------------------------------

def count_aliasing(alpha: float, W: int, n: int, A: float) -> int:
    """Count how many fiber shifts alpha+t/W (t=0..W-1) are on a major arc."""
    count = 0
    for t in range(W):
        if is_major_arc_shift(alpha + t / W, n, A):
            count += 1
    return count


# ---------------------------------------------------------------------------
# Self-check
# ---------------------------------------------------------------------------

def self_check(n: int, W: int, A: float) -> None:
    """Run validation checks for small N / W.

    1. Verify that the DFT via numpy FFT matches a direct DFT sum.
    2. Verify that F_hat(k*) matches an explicit sum over primes satisfying
       the congruence p1 + p2 ≡ N (mod W) (weighted by log(p1)*log(p2)).
    """
    print(f"\n--- Self-check: N={n}, W={W}, A={A} ---")
    primes = sieve_primes(n)
    if primes.size == 0:
        print("No primes found; skipping self-check.")
        return

    # Use a fixed alpha for the check (choose a clear minor-arc point)
    alpha = 0.0

    F = compute_fiber(alpha, W, primes)
    k_star = n % W

    # --- Check 1: FFT vs direct DFT ---
    F_hat_fft = np.fft.fft(F)
    t_vals = np.arange(W, dtype=np.float64)
    F_hat_direct = np.array([
        np.sum(F * np.exp(-2j * math.pi * k * t_vals / W))
        for k in range(W)
    ])
    max_err_1 = np.max(np.abs(F_hat_fft - F_hat_direct))
    print(f"  Check 1 (FFT vs direct DFT): max |error| = {max_err_1:.3e}",
          "✓" if max_err_1 < 1e-6 else "✗ FAIL")

    # --- Check 2: F_hat(k*) vs explicit prime-pair sum ---
    # F_hat(k*) = sum_t [ Phi(alpha+t/W)^2 ] * exp(-2pi i k* t / W)
    #           = sum_{p1,p2 <= N} log(p1) log(p2) * exp(2pi i alpha(p1+p2))
    #             * sum_t exp(2pi i (p1+p2-N) t/W) / W ... but simpler:
    # Because the inner DFT sum is W * [k == (p1+p2) mod W], we get:
    #   F_hat(k*) = W * sum_{p1+p2 ≡ N (mod W)} log(p1) log(p2)
    #               * exp(2pi i alpha (p1+p2))
    # For alpha=0 the exponential is 1, so:
    explicit = 0.0
    for p1 in primes:
        target_residue = int(n - p1) % W
        # find all p2 <= N with p2 ≡ target_residue (mod W)
        mask = (primes % W) == target_residue
        explicit += math.log(p1) * np.sum(np.log(primes[mask]))
    explicit *= W

    F_hat_k_fft = F_hat_fft[k_star]
    rel_err = abs(F_hat_k_fft - explicit) / (abs(explicit) + 1e-30)
    print(f"  Check 2 (F̂(k*) vs prime-pair sum): rel |error| = {rel_err:.3e}",
          "✓" if rel_err < 1e-4 else "✗ FAIL")

    print("--- End self-check ---\n")


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run_experiment(
    N: int,
    W: int,
    A: float,
    samples: int,
    seed: Optional[int],
    top_k: int,
    csv_path: Optional[str],
    self_check_mode: bool,
) -> None:
    rng = np.random.default_rng(seed)

    if self_check_mode:
        self_check(N, W, A)
        return

    print(f"Parameters: N={N}, W={W}, A={A}, samples={samples}, seed={seed}")

    primes = sieve_primes(N)
    print(f"Number of primes up to N: {primes.size}")
    if primes.size == 0:
        print("No primes found. Exiting.")
        return

    k_star = N % W
    print(f"k* = N mod W = {k_star}\n")

    rows = []
    sample_idx = 0
    attempts = 0
    # Allow up to this many random draws per requested sample before giving up.
    _MAX_ATTEMPTS_FACTOR = 10_000
    max_attempts = samples * _MAX_ATTEMPTS_FACTOR

    while sample_idx < samples and attempts < max_attempts:
        attempts += 1
        alpha = float(rng.uniform(0.0, 1.0))
        if is_major_arc(alpha, N, A):
            continue  # keep only minor-arc alphas

        # --- Compute Phi for all fiber shifts ---
        t_vals = np.arange(W, dtype=np.float64)
        shift_alphas = alpha + t_vals / W
        phi_vec = compute_phi_vector(shift_alphas, primes)

        # --- Fiber vector and FFT ---
        F = phi_vec ** 2
        F_hat = np.fft.fft(F)
        F_hat_abs = np.abs(F_hat)

        # --- Metrics ---
        aliasing = count_aliasing(alpha, W, N, A)
        fhat_kstar = float(F_hat_abs[k_star])

        # Top-K indices by |F̂(k)|
        top_indices = np.argsort(F_hat_abs)[::-1][:top_k]
        top_vals = F_hat_abs[top_indices]

        # Distribution stats of |Phi(alpha + t/W)|
        phi_abs = np.abs(phi_vec)

        phi_mean = float(phi_abs.mean())
        phi_std = float(phi_abs.std())
        phi_max = float(phi_abs.max())

        print(
            f"[{sample_idx + 1:3d}] alpha={alpha:.6f}  "
            f"B(alpha)={aliasing:4d}  "
            f"|F̂(k*)|={fhat_kstar:.4e}  "
            f"|Phi| mean/std/max={phi_mean:.3e}/{phi_std:.3e}/{phi_max:.3e}"
        )
        top_str = "  top-K: " + ", ".join(
            f"k={int(ki)}:{v:.3e}" for ki, v in zip(top_indices, top_vals)
        )
        print(top_str)

        rows.append({
            "sample": sample_idx + 1,
            "alpha": alpha,
            "aliasing_B": aliasing,
            "fhat_kstar": fhat_kstar,
            "phi_mean": phi_mean,
            "phi_std": phi_std,
            "phi_max": phi_max,
            "top_k_indices": ";".join(str(int(i)) for i in top_indices),
            "top_k_values": ";".join(f"{v:.4e}" for v in top_vals),
        })

        sample_idx += 1

    if sample_idx < samples:
        print(
            f"\nWarning: only collected {sample_idx}/{samples} minor-arc samples "
            f"after {attempts} attempts."
        )

    if csv_path:
        fieldnames = [
            "sample", "alpha", "aliasing_B", "fhat_kstar",
            "phi_mean", "phi_std", "phi_max",
            "top_k_indices", "top_k_values",
        ]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nResults written to {csv_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Goldbach fiber-sum experiment: major/minor-arc sampling and fiber DFT."
    )
    p.add_argument("--N", type=int, default=1_000_000,
                   help="Target even integer (default: 1_000_000)")
    p.add_argument("--W", type=int, default=1024,
                   help="Fiber width / DFT modulus (default: 1024)")
    p.add_argument("--A", type=float, default=6.0,
                   help="Major-arc exponent: P = (log N)^A (default: 6)")
    p.add_argument("--samples", type=int, default=50,
                   help="Number of minor-arc alpha samples (default: 50)")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducibility")
    p.add_argument("--top-k", type=int, default=5,
                   help="Number of top |F̂(k)| values to report (default: 5)")
    p.add_argument("--csv", type=str, default=None,
                   help="Write results to this CSV file")
    p.add_argument("--self-check", action="store_true",
                   help="Run self-check mode for small N/W then exit")
    return p


def main() -> None:
    args = build_parser().parse_args()

    if args.N <= 0 or args.N % 2 != 0:
        print("Error: N must be a positive even integer.", file=sys.stderr)
        sys.exit(1)
    if args.W <= 0:
        print("Error: W must be a positive integer.", file=sys.stderr)
        sys.exit(1)
    if args.A <= 0:
        print("Error: A must be positive.", file=sys.stderr)
        sys.exit(1)
    if args.samples <= 0:
        print("Error: samples must be positive.", file=sys.stderr)
        sys.exit(1)

    run_experiment(
        N=args.N,
        W=args.W,
        A=args.A,
        samples=args.samples,
        seed=args.seed,
        top_k=args.top_k,
        csv_path=args.csv,
        self_check_mode=args.self_check,
    )


if __name__ == "__main__":
    main()
