#!/usr/bin/env python3
"""
Goldbach Experiment Tooling
============================
Explores Goldbach's conjecture for a target even integer N by computing
fiber DFT coefficients on the reduced residue system modulo W (wheel).

Key relationship:
    k ≡ N (mod W)  — the residue class that governs which fiber carries
                      the Goldbach signal in the DFT spectrum.

Usage examples
--------------
# Default run (N=100, W=30, A=1.0)
python goldbach_experiment.py

# Custom target
python goldbach_experiment.py --N 200 --W 30 --A 1.0

# Larger wheel (W = 2·3·5·7 = 210)
python goldbach_experiment.py --N 1000 --W 210 --A 2.5

# Batch scan over many even integers, saving CSV output
python goldbach_experiment.py --N-start 4 --N-end 200 --W 30 --output results.csv
"""

import argparse
import csv
import math
import sys
from typing import List, Tuple


# ---------------------------------------------------------------------------
# Number-theory helpers
# ---------------------------------------------------------------------------

def is_prime(n: int) -> bool:
    """Return True if *n* is a prime number (trial-division, sufficient for
    the demo-scale values used in Goldbach experiments)."""
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


def reduced_residues(W: int) -> List[int]:
    """Return the reduced residue system modulo *W*: all r in [1, W) with
    gcd(r, W) = 1.  These form the "fiber" over which the DFT is computed."""
    return [r for r in range(1, W) if math.gcd(r, W) == 1]


def goldbach_pairs(N: int) -> List[Tuple[int, int]]:
    """Return all Goldbach decompositions of *N* as ordered pairs (p, q) with
    p ≤ q and p + q = N."""
    if N <= 2 or N % 2 != 0:
        return []
    pairs = []
    for p in range(2, N // 2 + 1):
        if is_prime(p) and is_prime(N - p):
            pairs.append((p, N - p))
    return pairs


# ---------------------------------------------------------------------------
# Fiber DFT
# ---------------------------------------------------------------------------

def fiber_dft_coefficient(N: int, W: int, A: float) -> complex:
    """Compute the fiber DFT Goldbach coefficient for residue class k ≡ N (mod W).

    The coefficient is defined as::

        C(N, W, A) = A · Σ_{r ∈ R(W)}  exp(2πi · r · k / φ(W))

    where

    * R(W) is the reduced residue system mod W (the "fiber"),
    * k = N mod W  is the target residue,
    * φ(W) = |R(W)| is Euler's totient of W,
    * A   is a user-supplied amplitude (scaling factor).

    The real part of the sum counts the constructive interference of the
    additive characters on the fiber, providing a spectral fingerprint of
    which residue class N belongs to.

    Parameters
    ----------
    N : int
        Target even integer (Goldbach candidate).
    W : int
        Wheel modulus (product of small primes, e.g. 30 = 2·3·5).
    A : float
        Amplitude / scaling coefficient.

    Returns
    -------
    complex
        The DFT coefficient C(N, W, A).
    """
    k = N % W  # k ≡ N (mod W)  — the residue that drives the coefficient
    fiber = reduced_residues(W)
    phi_W = len(fiber)
    if phi_W == 0:
        return complex(0)
    total = sum(
        math.e ** (2j * math.pi * r * k / phi_W)
        for r in fiber
    )
    return A * total


def fiber_dft_spectrum(N: int, W: int, A: float) -> List[Tuple[int, complex]]:
    """Return the full fiber DFT spectrum for all residues k in [0, W).

    Each entry is (k, C_k) where C_k is the DFT coefficient for residue k.
    The entry at k = N mod W is the Goldbach-relevant coefficient.
    """
    fiber = reduced_residues(W)
    phi_W = len(fiber)
    spectrum = []
    for k in range(W):
        coeff = A * sum(
            math.e ** (2j * math.pi * r * k / phi_W)
            for r in fiber
        )
        spectrum.append((k, coeff))
    return spectrum


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="goldbach_experiment",
        description=(
            "Goldbach experiment tooling: explore prime-pair decompositions "
            "and fiber DFT coefficients for a target even integer N.\n\n"
            "The key parameter relationship is:  k ≡ N (mod W)\n"
            "where W is the wheel modulus and k is the residue class that "
            "governs the Goldbach DFT coefficient."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # --- target integer(s) ---
    target = parser.add_mutually_exclusive_group()
    target.add_argument(
        "--N",
        type=int,
        default=100,
        metavar="N",
        help=(
            "Target even integer for a single Goldbach experiment "
            "(default: %(default)s). "
            "Must be an even integer ≥ 4."
        ),
    )
    target.add_argument(
        "--N-start",
        type=int,
        dest="N_start",
        metavar="N_START",
        help="Start of an inclusive even-integer batch range (use with --N-end).",
    )

    parser.add_argument(
        "--N-end",
        type=int,
        dest="N_end",
        metavar="N_END",
        help="End of an inclusive even-integer batch range (use with --N-start).",
    )

    # --- wheel ---
    parser.add_argument(
        "--W",
        type=int,
        default=30,
        metavar="W",
        help=(
            "Wheel modulus — product of the first few primes "
            "(default: %(default)s = 2·3·5). "
            "Common values: 6, 30, 210, 2310."
        ),
    )

    # --- amplitude ---
    parser.add_argument(
        "--A",
        type=float,
        default=1.0,
        metavar="A",
        help=(
            "Amplitude / scaling coefficient applied to every DFT term "
            "(default: %(default)s)."
        ),
    )

    # --- output ---
    parser.add_argument(
        "--output",
        metavar="FILE",
        default=None,
        help="Write results to a CSV file instead of stdout.",
    )
    parser.add_argument(
        "--spectrum",
        action="store_true",
        default=False,
        help="Print the full fiber DFT spectrum (all k in [0, W)).",
    )
    parser.add_argument(
        "--pairs",
        action="store_true",
        default=True,
        help="Print Goldbach prime pairs for N (default: True).",
    )
    parser.add_argument(
        "--no-pairs",
        action="store_false",
        dest="pairs",
        help="Suppress Goldbach prime-pair output.",
    )

    return parser


def validate_args(args: argparse.Namespace) -> None:
    """Raise SystemExit with a human-readable message on bad input."""
    if args.N_start is not None:
        if args.N_end is None:
            sys.exit("error: --N-end is required when --N-start is used.")
        if args.N_start > args.N_end:
            sys.exit("error: --N-start must be ≤ --N-end.")
        if args.N_start < 4:
            sys.exit("error: --N-start must be ≥ 4.")
        for label, val in [("--N-start", args.N_start), ("--N-end", args.N_end)]:
            if val % 2 != 0:
                sys.exit(f"error: {label} must be even (got {val}).")
    else:
        if args.N < 4 or args.N % 2 != 0:
            sys.exit(f"error: --N must be an even integer ≥ 4 (got {args.N}).")

    if args.W < 2:
        sys.exit(f"error: --W must be ≥ 2 (got {args.W}).")
    if args.A == 0.0:
        sys.exit("error: --A must be non-zero.")


def run_single(N: int, W: int, A: float, show_pairs: bool, show_spectrum: bool) -> dict:
    """Run a single Goldbach experiment and return a result dict."""
    k = N % W
    coeff = fiber_dft_coefficient(N, W, A)
    pairs = goldbach_pairs(N)  # always computed so num_pairs is accurate in CSV

    result = {
        "N": N,
        "W": W,
        "A": A,
        "k (= N mod W)": k,
        "phi(W)": len(reduced_residues(W)),
        "DFT_real": coeff.real,
        "DFT_imag": coeff.imag,
        "DFT_abs": abs(coeff),
        "num_pairs": len(pairs),
    }

    print(f"\n{'='*60}")
    print(f"  Goldbach Experiment")
    print(f"{'='*60}")
    print(f"  N        = {N}  (target even integer)")
    print(f"  W        = {W}  (wheel modulus)")
    print(f"  A        = {A}  (amplitude)")
    print(f"  k=N%W    = {k}  (k ≡ N (mod W) — governs DFT coefficient)")
    print(f"  φ(W)     = {result['phi(W)']}  (size of reduced residue fiber)")
    print(f"  DFT coeff= {coeff.real:.6f} + {coeff.imag:.6f}j  (|C| = {abs(coeff):.6f})")
    print(f"  # pairs  = {result['num_pairs']}")

    if show_pairs and pairs:
        print(f"\n  Goldbach pairs for N={N}:")
        for p, q in pairs:
            print(f"    {p} + {q} = {N}")

    if show_spectrum:
        print(f"\n  Fiber DFT spectrum (N={N}, W={W}, A={A}):")
        print(f"  {'k':>4}  {'Re(C)':>12}  {'Im(C)':>12}  {'|C|':>12}  (← marks k = N mod W)")
        for sk, sc in fiber_dft_spectrum(N, W, A):
            marker = "  ← k = N mod W" if sk == k else ""
            print(f"  {sk:>4}  {sc.real:>12.6f}  {sc.imag:>12.6f}  {abs(sc):>12.6f}{marker}")

    return result


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args)

    # Determine the list of N values to process
    if args.N_start is not None:
        n_values = list(range(args.N_start, args.N_end + 1, 2))
    else:
        n_values = [args.N]

    results = []
    for n in n_values:
        r = run_single(
            N=n,
            W=args.W,
            A=args.A,
            show_pairs=args.pairs and len(n_values) == 1,
            show_spectrum=args.spectrum and len(n_values) == 1,
        )
        results.append(r)

    # Batch summary
    if len(n_values) > 1:
        print(f"\n{'='*60}")
        print(f"  Batch summary: {len(n_values)} values, W={args.W}, A={args.A}")
        print(f"  {'N':>6}  {'k=N%W':>7}  {'|DFT C|':>12}  {'# pairs':>8}")
        print(f"  {'-'*6}  {'-'*7}  {'-'*12}  {'-'*8}")
        for r in results:
            print(
                f"  {r['N']:>6}  {r['k (= N mod W)']:>7}  "
                f"{r['DFT_abs']:>12.6f}  {r['num_pairs']:>8}"
            )

    # CSV output
    if args.output:
        with open(args.output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults written to: {args.output}")


if __name__ == "__main__":
    main()
