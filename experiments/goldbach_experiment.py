"""
Goldbach experiment: count representations of 2N as p+q (p,q prime)
using a discrete approximation of the Hardy-Littlewood circle method.

Cutoff convention
-----------------
By default ``cutoff_max = 2 * N``, which ensures that the exponential sum
S(α) = Σ_{p ≤ cutoff_max} e(p α) covers the full convolution support for
the representation identity

    r(2N) = ∫₀¹ S(α)² e(-2N α) dα.

Using ``cutoff_max = N`` restricts both primes to [2, N], which misses pairs
(p, q) where one prime exceeds N.  If you use this symmetric shortcut you
must multiply the returned count by 2 (see --help for details).

See docs/goldbach_circle_method.md for a full explanation.
"""

import argparse
import math


def sieve_of_eratosthenes(limit):
    """Return a sorted list of all primes up to *limit* (inclusive)."""
    if limit < 2:
        return []
    is_prime = bytearray([1] * (limit + 1))
    is_prime[0] = is_prime[1] = 0
    for i in range(2, math.isqrt(limit) + 1):
        if is_prime[i]:
            is_prime[i * i :: i] = bytearray(len(is_prime[i * i :: i]))
    return [i for i, v in enumerate(is_prime) if v]


def count_goldbach_representations(N, cutoff_max=None):
    """Count the number of ordered pairs (p, q) of primes with p + q = 2N.

    Parameters
    ----------
    N : int
        Target half-sum; representations of ``2*N`` are counted.
    cutoff_max : int or None
        Upper bound for the prime sieve.  Defaults to ``2 * N``, which is the
        correct support for the circle-method convolution.  Passing ``N``
        applies the symmetric shortcut and will under-count r(2N) by roughly
        half unless you multiply the result by 2.

    Returns
    -------
    count : int
        Number of ordered pairs (p, q) with p + q = 2*N, p ≤ cutoff_max,
        q ≤ cutoff_max.
    primes_used : list[int]
        Sorted list of primes up to *cutoff_max*.
    """
    if cutoff_max is None:
        cutoff_max = 2 * N  # default: full convolution support

    target = 2 * N
    primes = sieve_of_eratosthenes(cutoff_max)
    prime_set = set(primes)

    count = 0
    pairs = []
    for p in primes:
        q = target - p
        if q >= 2 and q in prime_set:
            count += 1
            pairs.append((p, q))

    return count, primes, pairs


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Count Goldbach representations of 2N as a sum of two primes. "
            "See docs/goldbach_circle_method.md for the cutoff convention."
        )
    )
    parser.add_argument(
        "N",
        type=int,
        help="Half-sum target; representations of 2*N are counted.",
    )
    parser.add_argument(
        "--cutoff-max",
        type=int,
        default=None,
        dest="cutoff_max",
        help=(
            "Upper bound for the prime sieve (default: 2*N). "
            "Setting this to N applies the symmetric shortcut: only pairs "
            "with both primes ≤ N are counted, roughly halving the result. "
            "If you use --cutoff-max=N you should multiply the output count "
            "by 2 to recover r(2N).  See docs/goldbach_circle_method.md."
        ),
    )
    parser.add_argument(
        "--show-pairs",
        action="store_true",
        help="Print every (p, q) pair found.",
    )
    args = parser.parse_args()

    N = args.N
    cutoff_max = args.cutoff_max if args.cutoff_max is not None else 2 * N

    if cutoff_max == N:
        print(
            f"[WARNING] cutoff_max={cutoff_max} equals N={N}. "
            "This is the symmetric shortcut: pairs where one prime > N are "
            "excluded.  Multiply the count by 2 to estimate r(2N).  "
            "See docs/goldbach_circle_method.md for details."
        )

    count, primes, pairs = count_goldbach_representations(N, cutoff_max=cutoff_max)

    print(f"Target : 2N = {2 * N}")
    print(f"Cutoff : cutoff_max = {cutoff_max}  (default is 2*N = {2 * N})")
    print(f"Primes sieved up to {cutoff_max}: {len(primes)} primes found")
    print(f"r({2 * N}) = {count}  (ordered pairs (p, q) with p + q = {2 * N})")

    if args.show_pairs:
        print("Pairs:")
        for p, q in pairs:
            print(f"  {p} + {q} = {2 * N}")


if __name__ == "__main__":
    main()
