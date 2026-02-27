"""
Goldbach Fiber Experiments
==========================
Numerical exploration of the circle method's fiber/W-trick for the
Goldbach conjecture.  Computes Phi(alpha + t/W) on a discrete grid,
measures minor-arc aliasing counts, and prints basic estimates.

Usage
-----
    python tools/goldbach_fiber_experiments.py [--N 200] [--plot]

Requirements
------------
    Python >= 3.8, numpy, (optionally matplotlib for --plot)
"""

import argparse
import cmath
import math
from typing import List, Tuple


# ---------------------------------------------------------------------------
# Arithmetic helpers
# ---------------------------------------------------------------------------

def primes_up_to(n: int) -> List[int]:
    """Return all primes <= n via a simple sieve of Eratosthenes."""
    if n < 2:
        return []
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i in range(2, n + 1) if sieve[i]]


def primorial(w: int) -> int:
    """Return W = product of primes <= w."""
    result = 1
    for p in primes_up_to(w):
        result *= p
    return result


def von_mangoldt(n: int) -> float:
    """Return Lambda(n): log(p) if n is a prime power p^k, else 0."""
    if n < 2:
        return 0.0
    for p in primes_up_to(n):
        if p * p > n:
            break
        if n % p == 0:
            while n % p == 0:
                n //= p
            return math.log(p) if n == 1 else 0.0
    return math.log(n)


# ---------------------------------------------------------------------------
# Core exponential sum
# ---------------------------------------------------------------------------

def S_alpha(alpha: float, N: int) -> complex:
    """
    Compute S(alpha) = sum_{n=1}^{N} Lambda(n) * e(n * alpha)
    where e(x) = exp(2*pi*i*x).
    """
    two_pi_i = 2j * math.pi
    total = 0j
    for n in range(2, N + 1):
        lam = von_mangoldt(n)
        if lam:
            total += lam * cmath.exp(two_pi_i * n * alpha)
    return total


def S_alpha_vec(alpha: float, N: int) -> complex:
    """Vectorised version using numpy if available, else fallback."""
    try:
        import numpy as np  # optional fast path
        ns = np.arange(2, N + 1, dtype=float)
        lam = np.array([von_mangoldt(int(n)) for n in ns])
        return complex(np.sum(lam * np.exp(2j * math.pi * ns * alpha)))
    except ImportError:
        return S_alpha(alpha, N)


# ---------------------------------------------------------------------------
# W-trick residue classes
# ---------------------------------------------------------------------------

def coprime_residues(W: int) -> List[int]:
    """Return all b in [1, W) with gcd(b, W) == 1."""
    return [b for b in range(1, W) if math.gcd(b, W) == 1]


def S_fiber(b: int, W: int, alpha: float, N: int) -> complex:
    """
    Fiber sum S_b(alpha) = sum_{n<=N, n≡b (W)} Lambda(n) * e(n*alpha).
    """
    two_pi_i = 2j * math.pi
    total = 0j
    for n in range(b, N + 1, W):
        lam = von_mangoldt(n)
        if lam:
            total += lam * cmath.exp(two_pi_i * n * alpha)
    return total


# ---------------------------------------------------------------------------
# Grid experiment
# ---------------------------------------------------------------------------

def aliasing_count(t: int, W: int, N: int) -> int:
    """
    A(t) = #{(b,k) : b*k ≡ t (mod W), gcd(b,W)=1, 1<=b*k<=N}
    Counts constructive overlap of fiber frequencies at grid point t.
    """
    count = 0
    for b in coprime_residues(W):
        k = 1
        while b * k <= N:
            if (b * k) % W == t % W:
                count += 1
            k += 1
    return count


def run_grid_experiment(N: int, w_param: int = 6) -> Tuple[List[Tuple], int]:
    """
    Evaluate |S(t/W)| for t=0,...,W-1 and collect aliasing info.

    Parameters
    ----------
    N       : even integer to study
    w_param : smooth bound for primorial (default: 6 → W=2*3*5=30)

    Returns
    -------
    List of (t, alpha, |S(alpha)|, A(t)) tuples, sorted by |S| descending.
    """
    W = primorial(w_param)
    print(f"N = {N},  W = primorial({w_param}) = {W}")
    print(f"  phi(W) = {len(coprime_residues(W))} coprime residues")

    results = []
    for t in range(W):
        alpha = t / W
        s_val = abs(S_alpha_vec(alpha, N))
        a_count = aliasing_count(t, W, N)
        results.append((t, alpha, s_val, a_count))

    results.sort(key=lambda x: x[2], reverse=True)
    return results, W


def goldbach_count(N: int) -> int:
    """
    Count r(N) = #{(p,q) : p+q=N, p,q prime}, direct brute-force.
    """
    ps = set(primes_up_to(N))
    return sum(1 for p in ps if (N - p) in ps and p <= N - p)


# ---------------------------------------------------------------------------
# Singular series estimate
# ---------------------------------------------------------------------------

def singular_series_estimate(N: int, prime_limit: int = 200) -> float:
    """
    Approximate the Hardy–Littlewood singular series S(N).
    S(N) ≈ 2*C_2 * prod_{p|N, p>2} (p-1)/(p-2)
    C_2 = prod_{p>2} [1 - 1/(p-1)^2]  (twin-prime constant)
    """
    ps = primes_up_to(prime_limit)
    # twin-prime-like constant
    C2 = 1.0
    for p in ps:
        if p > 2:
            C2 *= 1.0 - 1.0 / (p - 1) ** 2
    # local factors for primes dividing N
    local = 1.0
    for p in ps:
        if p > 2 and N % p == 0:
            local *= (p - 1) / (p - 2)
    return 2.0 * C2 * local


def hardy_littlewood_estimate(N: int) -> float:
    """
    Hardy–Littlewood asymptotic estimate for R(N):
      R_HL(N) ≈ S(N) * N / (log N)^2
    """
    sg = singular_series_estimate(N)
    return sg * N / (math.log(N) ** 2)


# ---------------------------------------------------------------------------
# Minor arc numerical estimate
# ---------------------------------------------------------------------------

def minor_arc_l2_estimate(N: int, W: int, num_points: int = 100) -> float:
    """
    Rough numerical estimate of integral_m |S(alpha)|^2 d alpha
    via a random sample on [0,1] excluding major arcs.

    Only grid points with t/W not close to a/q for small q are used.
    This is illustrative, not rigorous.
    """
    import random
    random.seed(42)
    q_bound = int(math.log(N) ** 2) + 1
    major_set = set()
    for q in range(1, q_bound + 1):
        for a in range(q + 1):
            major_set.add((a, q))

    def is_major(alpha: float) -> bool:
        for q in range(1, q_bound + 1):
            for a in range(q + 1):
                if abs(alpha - a / q) < (math.log(N) ** 2) / N:
                    return True
        return False

    sample_sum = 0.0
    count = 0
    for _ in range(num_points):
        alpha = random.random()
        if not is_major(alpha):
            sample_sum += abs(S_alpha_vec(alpha, N)) ** 2
            count += 1

    return sample_sum / count if count else 0.0


# ---------------------------------------------------------------------------
# Pretty-print and plot
# ---------------------------------------------------------------------------

def print_results(results: List[Tuple], W: int, N: int, top_k: int = 10) -> None:
    print(f"\n{'='*60}")
    print(f"Top {top_k} grid points by |S(t/W)|  (N={N}, W={W})")
    print(f"{'t':>5}  {'alpha':>8}  {'|S|':>12}  {'A(t)':>6}")
    print("-" * 40)
    for t, alpha, s_val, a_count in results[:top_k]:
        print(f"{t:>5}  {alpha:8.4f}  {s_val:12.4f}  {a_count:>6}")

    hl = hardy_littlewood_estimate(N)
    r_direct = goldbach_count(N)
    print(f"\nHardy–Littlewood estimate R_HL({N}) ≈ {hl:.2f}")
    print(f"Direct count R({N})            = {r_direct}")


def plot_results(results: List[Tuple], W: int, N: int) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping plot.")
        return

    ts = [r[0] for r in results]
    s_vals = [r[2] for r in results]
    a_counts = [r[3] for r in results]

    # Sort back by t for plotting
    order = sorted(range(len(ts)), key=lambda i: ts[i])
    ts_sorted = [ts[i] for i in order]
    s_sorted = [s_vals[i] for i in order]
    a_sorted = [a_counts[i] for i in order]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax1.bar(ts_sorted, s_sorted, color="steelblue")
    ax1.set_ylabel("|S(t/W)|")
    ax1.set_title(f"Fiber grid experiment  N={N}, W={W}")

    ax2.bar(ts_sorted, a_sorted, color="coral")
    ax2.set_ylabel("Aliasing count A(t)")
    ax2.set_xlabel("t  (alpha = t/W)")

    plt.tight_layout()
    out = f"/tmp/goldbach_N{N}_W{W}.png"
    plt.savefig(out)
    print(f"Plot saved to {out}")
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Goldbach circle method fiber/W-trick grid experiment"
    )
    parser.add_argument("--N", type=int, default=200,
                        help="Even integer to study (default: 200)")
    parser.add_argument("--w", type=int, default=6,
                        help="Primorial parameter: W = prod of primes <= w (default: 6 → W=30)")
    parser.add_argument("--top", type=int, default=10,
                        help="Number of top grid points to display (default: 10)")
    parser.add_argument("--plot", action="store_true",
                        help="Show/save a bar chart (requires matplotlib)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    N = args.N
    if N < 4 or N % 2 != 0:
        N = max(4, N + (N % 2))
        print(f"Warning: N adjusted to {N} (must be even and >= 4).")

    results, W = run_grid_experiment(N, w_param=args.w)
    print_results(results, W, N, top_k=args.top)

    if args.plot:
        plot_results(results, W, N)


if __name__ == "__main__":
    main()
