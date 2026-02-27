"""
goldbach_arcs.py
================
Helpers for the Hardy-Littlewood circle-method experiment on Goldbach's conjecture.

Major/minor arcs definition (case 1)
-------------------------------------
Given a target even integer N and exponent A, set

    P = (log N)^A

The *major arc* around the rational a/q is

    M(a, q) = { alpha in [0,1) : |alpha - a/q| <= P / (N * q) }

Major arcs M are the union of M(a, q) over all integers a, q with
    1 <= q <= P,  0 <= a < q,  gcd(a, q) = 1   (plus a=0, q=1).

The *minor arcs* are the complement m = [0,1) \\ M.

The exponential sum studied is
    Phi(alpha) = sum_{p <= N} log(p) * e(alpha * p)
where e(x) = exp(2*pi*i*x) and the sum runs over primes p.

W-trick / fiber DFT
--------------------
With smoothing modulus W the fiber DFT coefficient index is
    k = N mod W
because we want representations N = p + q.

Aliasing count
--------------
For fixed alpha and modulus W, the aliasing count is the number of
    t in {0, 1, ..., W-1}  such that  alpha + t/W  is in M.
"""

import math
import cmath
import argparse
from math import gcd, log, pi


# ---------------------------------------------------------------------------
# Primality helpers
# ---------------------------------------------------------------------------

def _sieve(n: int) -> list:
    """Return a boolean list `is_prime` of length n+1 (Sieve of Eratosthenes)."""
    if n < 2:
        return [False] * (n + 1)
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = False
    return is_prime


# ---------------------------------------------------------------------------
# Major arc classifier
# ---------------------------------------------------------------------------

def _arc_parameter(N: int, A: float) -> float:
    """Return P = (log N)^A."""
    return log(N) ** A


def is_major_arc(alpha: float, N: int, A: float) -> bool:
    """
    Return True if *alpha* (real, any value) belongs to the major arcs M
    defined by P = (log N)^A.

    The test is done modulo 1 so alpha is first reduced to [0, 1).

    A point alpha is in M if there exist integers a, q with
        1 <= q <= P,  gcd(a, q) = 1,  0 <= a < q,  |alpha - a/q| <= P / (N * q)
    """
    alpha = alpha % 1.0
    P = _arc_parameter(N, A)
    P_floor = int(P)

    for q in range(1, P_floor + 1):
        for a in range(0, q):
            if gcd(a, q) != 1:
                continue
            center = a / q
            # Distance on the circle [0,1): use minimum of |alpha-center|
            # and 1-|alpha-center| to handle wrap-around.
            dist = abs(alpha - center)
            dist = min(dist, 1.0 - dist)
            threshold = P / (N * q)
            if dist <= threshold:
                return True
    return False


def classify_alpha(alpha: float, N: int, A: float) -> str:
    """Return 'major' or 'minor' for the given alpha."""
    return "major" if is_major_arc(alpha, N, A) else "minor"


def closest_rational(alpha: float, N: int, A: float):
    """
    Return the tuple (a, q, dist, threshold) for the major-arc rational
    closest to alpha (mod 1), or None if alpha is in the minor arcs.
    """
    alpha = alpha % 1.0
    P = _arc_parameter(N, A)
    P_floor = int(P)

    best = None
    for q in range(1, P_floor + 1):
        for a in range(0, q):
            if gcd(a, q) != 1:
                continue
            center = a / q
            dist = abs(alpha - center)
            dist = min(dist, 1.0 - dist)
            threshold = P / (N * q)
            if dist <= threshold:
                if best is None or dist < best[2]:
                    best = (a, q, dist, threshold)
    return best


# ---------------------------------------------------------------------------
# Fiber DFT coefficient
# ---------------------------------------------------------------------------

def fiber_dft_k(N: int, W: int) -> int:
    """
    Return the fiber DFT coefficient index k = N mod W.

    When using the W-trick with smoothing modulus W, the relevant Fourier
    coefficient for counting representations N = p + q is the one at
    frequency k ≡ N (mod W).
    """
    return N % W


# ---------------------------------------------------------------------------
# Aliasing count
# ---------------------------------------------------------------------------

def count_aliasing(alpha: float, W: int, N: int, A: float) -> int:
    """
    Return the number of shifts t in {0, 1, ..., W-1} such that
    alpha + t/W lies in the major arcs M(N, A).
    """
    count = 0
    for t in range(W):
        if is_major_arc(alpha + t / W, N, A):
            count += 1
    return count


# ---------------------------------------------------------------------------
# Exponential sum Phi(alpha)
# ---------------------------------------------------------------------------

def phi_alpha(alpha: float, N: int) -> complex:
    """
    Compute  Phi(alpha) = sum_{p <= N} log(p) * e(alpha * p)
    where e(x) = exp(2*pi*i*x).
    """
    is_prime = _sieve(N)
    two_pi_i = 2j * pi
    total = 0.0 + 0j
    for p in range(2, N + 1):
        if is_prime[p]:
            total += log(p) * cmath.exp(two_pi_i * alpha * p)
    return total


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Classify a frequency alpha as major or minor arc "
            "for the Goldbach Hardy-Littlewood experiment."
        )
    )
    parser.add_argument(
        "--N",
        type=int,
        required=True,
        help="Target even integer N.",
    )
    parser.add_argument(
        "--A",
        type=float,
        default=2.0,
        help="Exponent A in P = (log N)^A  (default: 2.0).",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        required=True,
        help="Value of alpha to classify.",
    )
    parser.add_argument(
        "--W",
        type=int,
        default=1,
        help="Smoothing modulus W for the W-trick (default: 1).",
    )
    parser.add_argument(
        "--phi",
        action="store_true",
        help="Also compute Phi(alpha) = sum_{p<=N} log(p) e(alpha p). "
             "Warning: slow for large N.",
    )
    return parser


def main():
    parser = _build_parser()
    args = parser.parse_args()

    N = args.N
    A = args.A
    alpha = args.alpha
    W = args.W

    if N < 4 or N % 2 != 0:
        parser.error("N must be an even integer >= 4.")
    if A <= 0:
        parser.error("A must be positive.")
    if W < 1:
        parser.error("W must be >= 1.")

    P = _arc_parameter(N, A)
    label = classify_alpha(alpha, N, A)
    k = fiber_dft_k(N, W)
    alias = count_aliasing(alpha, W, N, A)

    print(f"N = {N},  A = {A},  P = (log N)^A = {P:.4f}")
    print(f"alpha = {alpha}")
    print(f"Classification: {label.upper()}")

    rat = closest_rational(alpha, N, A)
    if rat is not None:
        a, q, dist, threshold = rat
        print(
            f"  Closest major-arc rational: {a}/{q},  "
            f"|alpha - a/q| = {dist:.4e},  "
            f"threshold P/(Nq) = {threshold:.4e}"
        )

    print(f"Fiber DFT coefficient index: k = N mod W = {N} mod {W} = {k}")
    print(
        f"Aliasing count (t in [0,{W}) s.t. alpha+t/W in major arcs): {alias}"
    )

    if args.phi:
        val = phi_alpha(alpha, N)
        print(f"Phi(alpha) = {val.real:.6f} + {val.imag:.6f}i  "
              f"(|Phi| = {abs(val):.6f})")


if __name__ == "__main__":
    main()
