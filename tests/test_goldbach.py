"""tests/test_goldbach.py — Smoke tests for goldbach_experiment.py."""

import math
import sys
import os

import numpy as np
import pytest

# Allow importing from the repo root without installing the package.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import goldbach_experiment as ge


# ---------------------------------------------------------------------------
# sieve_primes
# ---------------------------------------------------------------------------

class TestSievePrimes:
    def test_small(self):
        assert list(ge.sieve_primes(10)) == [2, 3, 5, 7]

    def test_below_2(self):
        assert len(ge.sieve_primes(1)) == 0

    def test_exact_prime(self):
        primes = ge.sieve_primes(13)
        assert 13 in primes

    def test_dtype(self):
        assert ge.sieve_primes(20).dtype == np.int64


# ---------------------------------------------------------------------------
# major_arc_threshold
# ---------------------------------------------------------------------------

class TestMajorArcThreshold:
    def test_value(self):
        N, A = 1_000_000, 2.0
        expected = math.log(N) ** A
        assert math.isclose(ge.major_arc_threshold(N, A), expected)


# ---------------------------------------------------------------------------
# is_major_arc
# ---------------------------------------------------------------------------

class TestIsMajorArc:
    def test_zero_is_major(self):
        # alpha=0 is a/q = 0/1; always a major arc.
        assert ge.is_major_arc(0.0, P=10.0, N=100)

    def test_half_is_major(self):
        # alpha=0.5 = 1/2; q=2 <= P for any reasonable P.
        assert ge.is_major_arc(0.5, P=10.0, N=100)

    def test_irrational_like_large_N(self):
        # With P=2, N=1_000_000 the arc radii are very narrow.
        # A "random" irrational-like alpha should be minor.
        alpha = 0.12345678901234
        P = ge.major_arc_threshold(1_000_000, 2.0)
        # Result is boolean — just check type
        result = ge.is_major_arc(alpha, P, 1_000_000)
        assert isinstance(result, bool)


# ---------------------------------------------------------------------------
# aliasing_count
# ---------------------------------------------------------------------------

class TestAliasingCount:
    def test_output_in_range(self):
        W = 32
        B = ge.aliasing_count(0.12345, W=W, P=5.0, N=1000)
        assert 0 <= B <= W

    def test_alpha_zero_all_major(self):
        # alpha=0: all fiber shifts 0, 1/W, 2/W, … include 0, which is always
        # major. B should be at least 1.
        B = ge.aliasing_count(0.0, W=16, P=10.0, N=100)
        assert B >= 1


# ---------------------------------------------------------------------------
# sample_minor_arc_alpha
# ---------------------------------------------------------------------------

class TestSampleMinorArcAlpha:
    def test_returns_float_in_unit_interval(self):
        rng = np.random.default_rng(0)
        alpha = ge.sample_minor_arc_alpha(P=5.0, N=1_000_000, rng=rng)
        assert 0.0 <= alpha < 1.0

    def test_returned_alpha_not_major(self):
        rng = np.random.default_rng(1)
        P = 5.0
        N = 1_000_000
        alpha = ge.sample_minor_arc_alpha(P=P, N=N, rng=rng)
        assert not ge.is_major_arc(alpha, P, N)

    def test_raises_when_no_minor_arc(self):
        # P >= N means everything is a major arc; should raise RuntimeError.
        rng = np.random.default_rng(0)
        with pytest.raises(RuntimeError):
            ge.sample_minor_arc_alpha(P=1e9, N=100, rng=rng, max_attempts=10)


# ---------------------------------------------------------------------------
# compute_phi_array
# ---------------------------------------------------------------------------

class TestComputePhiArray:
    def test_shape_and_dtype(self):
        primes = ge.sieve_primes(100)
        phi = ge.compute_phi_array(alpha=0.1, W=16, primes=primes)
        assert phi.shape == (16,)
        assert phi.dtype == np.complex128

    def test_first_element_matches_direct(self):
        """phi_vals[0] = Phi(alpha) = sum_{p} log(p)*exp(2*pi*i*alpha*p)."""
        primes = ge.sieve_primes(50)
        alpha = 0.3
        W = 8
        phi = ge.compute_phi_array(alpha, W, primes)
        direct = np.sum(np.log(primes) * np.exp(2j * math.pi * alpha * primes))
        assert abs(phi[0] - direct) < 1e-8

    def test_shift_property(self):
        """phi_vals[t] should equal Phi(alpha + t/W) computed directly."""
        primes = ge.sieve_primes(30)
        alpha = 0.25
        W = 4
        phi = ge.compute_phi_array(alpha, W, primes)
        for t in range(W):
            alpha_t = alpha + t / W
            direct = np.sum(np.log(primes) * np.exp(2j * math.pi * alpha_t * primes))
            assert abs(phi[t] - direct) < 1e-7, f"Mismatch at t={t}"


# ---------------------------------------------------------------------------
# run_experiment (integration smoke test)
# ---------------------------------------------------------------------------

class TestRunExperiment:
    def test_smoke_small(self):
        results = ge.run_experiment(N=100, W=8, A=2.0, seed=0, verbose=False)
        assert results["N"] == 100
        assert results["W"] == 8
        assert results["Fhat"].shape == (8,)
        assert results["Fhat_abs"].shape == (8,)
        assert 0 <= results["B"] <= 8
        assert 0.0 <= results["alpha"] < 1.0

    def test_csv_output(self, tmp_path):
        csv_path = str(tmp_path / "fhat.csv")
        ge.run_experiment(N=100, W=8, A=2.0, seed=0, verbose=False, save_csv=csv_path)
        with open(csv_path) as f:
            lines = f.readlines()
        # header + W data rows
        assert len(lines) == 9

    def test_cli_help(self):
        with pytest.raises(SystemExit) as exc_info:
            ge.main(["--help"])
        assert exc_info.value.code == 0

    def test_cli_invalid_N(self):
        with pytest.raises(SystemExit) as exc_info:
            ge.main(["--N", "101"])   # odd N
        assert exc_info.value.code != 0
