"""
Unit tests for goldbach_experiment.py
"""

import sys
import os

import numpy as np
import pytest

# Allow running from the repo root without installing the package
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from goldbach_experiment import (
    sieve_primes,
    compute_phi,
    run_experiment,
    self_check,
    random_sample_sweep,
)


# ---------------------------------------------------------------------------
# sieve_primes
# ---------------------------------------------------------------------------

class TestSievePrimes:
    def test_empty_below_2(self):
        assert sieve_primes(0).tolist() == []
        assert sieve_primes(1).tolist() == []

    def test_small_values(self):
        assert sieve_primes(2).tolist() == [2]
        assert sieve_primes(10).tolist() == [2, 3, 5, 7]

    def test_known_list_up_to_30(self):
        expected = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        assert sieve_primes(30).tolist() == expected

    def test_all_prime(self):
        primes = sieve_primes(100)
        for p in primes:
            assert p >= 2
            for d in range(2, int(p**0.5) + 1):
                assert p % d != 0, f"{p} is not prime"

    def test_no_composites_included(self):
        primes = set(sieve_primes(50).tolist())
        for n in range(2, 51):
            is_prime = all(n % d != 0 for d in range(2, n))
            assert (n in primes) == is_prime

    def test_count_up_to_1000(self):
        # π(1000) = 168
        assert len(sieve_primes(1000)) == 168


# ---------------------------------------------------------------------------
# compute_phi
# ---------------------------------------------------------------------------

class TestComputePhi:
    def test_phi_at_zero_weighted(self):
        primes = sieve_primes(30)
        phi = compute_phi(primes, np.array([0.0]))
        expected = np.log(primes).sum()
        assert abs(phi[0].real - expected) < 1e-9
        assert abs(phi[0].imag) < 1e-9

    def test_phi_at_zero_unweighted(self):
        primes = sieve_primes(30)
        phi = compute_phi(primes, np.array([0.0]), weighted=False)
        assert abs(phi[0].real - len(primes)) < 1e-9
        assert abs(phi[0].imag) < 1e-9

    def test_phi_shape(self):
        primes = sieve_primes(20)
        thetas = np.linspace(0, 1, 16, endpoint=False)
        phi = compute_phi(primes, thetas)
        assert phi.shape == (16,)

    def test_phi_conjugate_symmetry(self):
        # Φ(-θ) = conj(Φ(θ)) because primes are integers
        primes = sieve_primes(30)
        theta = 0.123456
        phi_pos = compute_phi(primes, np.array([theta]))[0]
        phi_neg = compute_phi(primes, np.array([-theta]))[0]
        assert abs(phi_pos - phi_neg.conjugate()) < 1e-9


# ---------------------------------------------------------------------------
# run_experiment
# ---------------------------------------------------------------------------

class TestRunExperiment:
    def test_returns_dict(self):
        r = run_experiment(30, 8, 1, 6, 0.0, verbose=False)
        assert isinstance(r, dict)

    def test_goldbach_k_formula(self):
        for N, W in [(50, 8), (100, 16), (1000, 32)]:
            r = run_experiment(N, W, 1, 6, 0.0, verbose=False)
            assert r["goldbach_k"] == (2 * N) % W

    def test_top_k_length(self):
        r = run_experiment(50, 8, 1, 6, 0.0, topk=3, verbose=False)
        assert len(r["top_k_indices"]) == 3
        assert len(r["top_k_magnitudes"]) == 3

    def test_top_k_sorted_descending(self):
        r = run_experiment(50, 16, 1, 6, 0.0, topk=5, verbose=False)
        mags = r["top_k_magnitudes"]
        assert list(mags) == sorted(mags, reverse=True)

    def test_phi_stats_keys(self):
        r = run_experiment(50, 8, 1, 6, 0.0, verbose=False)
        for key in ("min", "max", "mean", "median", "std"):
            assert key in r["phi_stats"]

    def test_parseval_identity(self):
        r = run_experiment(50, 8, 1, 6, 0.0, verbose=False)
        F, F_hat, W = r["F"], r["F_hat"], r["W"]
        lhs = np.sum(np.abs(F_hat) ** 2)
        rhs = W * np.sum(np.abs(F) ** 2)
        assert np.isclose(lhs, rhs, rtol=1e-9)

    def test_goldbach_magnitude_non_negative(self):
        r = run_experiment(100, 16, 1, 6, 0.0, verbose=False)
        assert r["goldbach_magnitude"] >= 0.0

    def test_aliasing_count_non_negative(self):
        r = run_experiment(100, 16, 1, 6, 0.0, verbose=False)
        assert r["aliasing_count"] >= 0

    def test_f_equals_phi_squared(self):
        r = run_experiment(30, 8, 1, 6, 0.0, verbose=False)
        assert np.allclose(r["F"], r["phi"] ** 2)


# ---------------------------------------------------------------------------
# random_sample_sweep
# ---------------------------------------------------------------------------

class TestRandomSampleSweep:
    def test_returns_correct_count(self):
        results = random_sample_sweep(30, 8, 1, 6, n_samples=5, verbose=False)
        assert len(results) == 5

    def test_deterministic_with_seed(self):
        r1 = random_sample_sweep(30, 8, 1, 6, n_samples=3, seed=7, verbose=False)
        r2 = random_sample_sweep(30, 8, 1, 6, n_samples=3, seed=7, verbose=False)
        for a, b in zip(r1, r2):
            assert a["beta"] == b["beta"]


# ---------------------------------------------------------------------------
# self_check
# ---------------------------------------------------------------------------

class TestSelfCheck:
    def test_self_check_passes(self):
        assert self_check() is True
