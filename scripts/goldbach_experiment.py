"""
goldbach_experiment.py
======================
Numerical experiment for the binary Goldbach conjecture using the circle method.

Convention (default)
--------------------
  N           – target even integer (e.g. 28)
  cutoff_max  = N          (sum runs n = 1 … N)
  W           = N          (DFT grid size)
  k           = N % W = 0  (DFT frequency bin for the target)

The von-Mangoldt weighted representation count is approximated as:

    r_Λ(N) ≈ DFT[ S · S ][k]

where S[j] = Σ_{n=1}^{cutoff_max} Λ(n) · e(2πi · n · j / W).

See docs/goldbach_convention.md for full derivation and the alternative
convention (use_half_convention=True).
"""

import math
import cmath


# ---------------------------------------------------------------------------
# von Mangoldt function
# ---------------------------------------------------------------------------

def von_mangoldt(n: int) -> float:
    """Return Λ(n): log(p) if n = p^r for some prime p and r ≥ 1, else 0."""
    if n < 2:
        return 0.0
    # trial-divide to find whether n is a prime power
    for p in range(2, int(math.isqrt(n)) + 1):
        if n % p == 0:
            # p is the smallest prime factor; check if n is a power of p
            while n % p == 0:
                n //= p
            return math.log(p) if n == 1 else 0.0
    # n is prime
    return math.log(n)


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_representation_count(
    target_even: int,
    use_half_convention: bool = False,
) -> float:
    """
    Approximate r_Λ(target_even) via a length-W DFT.

    Parameters
    ----------
    target_even : int
        The even integer N to represent as a sum of two primes.
    use_half_convention : bool
        If False (default): cutoff_max = N, W = N, k = N % W.
        If True:            cutoff_max = N//2, W = N, k = N % W.
        Both conventions pick DFT bin k = 0 for W = N.
        See docs/goldbach_convention.md §4 for details.

    Returns
    -------
    float
        Real part of DFT[S·S][k] (imaginary part should be ≈ 0).
    """
    N = target_even

    if N < 4 or N % 2 != 0:
        raise ValueError(f"target_even must be an even integer ≥ 4, got {N}")

    # --- convention parameters ---
    if use_half_convention:
        cutoff_max = N // 2          # alternative: sum only to N/2
        label = f"M = N//2 = {cutoff_max}  [half-convention]"
    else:
        cutoff_max = N               # default: sum to N
        label = f"N = {N}  [default convention]"

    W = N                            # DFT grid size
    k = N % W                        # frequency bin for target → 0

    # --- build S[j] for j = 0, …, W-1 ---
    S = [0.0 + 0.0j] * W
    for j in range(W):
        for n in range(1, cutoff_max + 1):
            lam = von_mangoldt(n)
            if lam == 0.0:
                continue
            S[j] += lam * cmath.exp(2j * math.pi * n * j / W)

    # --- DFT of S·S at bin k ---
    SS_k = 0.0 + 0.0j
    for j in range(W):
        SS_k += (S[j] ** 2) * cmath.exp(-2j * math.pi * k * j / W)
    SS_k /= W

    return SS_k.real, label, cutoff_max, k


# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------

def print_formula(target_even: int, use_half_convention: bool = False) -> None:
    N = target_even
    cutoff = N // 2 if use_half_convention else N
    W = N
    k = N % W
    print("=" * 60)
    print(f"  Target even integer  N = {N}")
    print(f"  Convention           : {'half (target=2M)' if use_half_convention else 'default (target=N)'}")
    print(f"  cutoff_max           = {cutoff}")
    print(f"  DFT grid W           = {W}")
    print(f"  DFT frequency bin k  = {k}  (k ≡ N mod W)")
    print()
    print(f"  S(α) = Σ_{{n=1}}^{{{cutoff}}} Λ(n)·e(αn)")
    print(f"  r_Λ({N}) ≈ DFT[S·S][k={k}]")
    print("=" * 60)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_experiment(target_even: int, use_half_convention: bool = False) -> None:
    print_formula(target_even, use_half_convention)
    value, label, cutoff_max, k = compute_representation_count(
        target_even, use_half_convention
    )
    print(f"\n  cutoff_max = {cutoff_max}  ({label})")
    print(f"  r_Λ({target_even}) ≈ {value:.6f}")
    print()
    # Cross-check: direct prime-pair count
    direct = sum(
        1
        for p in range(2, target_even)
        if _is_prime(p) and _is_prime(target_even - p)
    )
    print(f"  Direct prime-pair count (unweighted): {direct}")
    print()


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(math.isqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Goldbach circle-method experiment (N = target even integer)."
    )
    parser.add_argument(
        "N",
        type=int,
        nargs="?",
        default=28,
        help="Target even integer N (default: 28).",
    )
    parser.add_argument(
        "--half-convention",
        action="store_true",
        default=False,
        help=(
            "Use the alternative half-convention: cutoff_max = N//2 instead of N. "
            "See docs/goldbach_convention.md §4."
        ),
    )
    args = parser.parse_args()
    run_experiment(args.N, use_half_convention=args.half_convention)
