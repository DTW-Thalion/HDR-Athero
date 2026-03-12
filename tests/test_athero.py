"""Tests for the HDR-Athero parameter set and simulation."""
from __future__ import annotations

import numpy as np
import pytest

from hdr_athero.model.parameters import (
    A, B, J, Q, TAU, N_STATE, N_CONTROL, N_BASINS,
    Q_COST, R_COST, TARGET_RADIUS, P_TRANS,
)
from hdr_athero.model.slds import AtheroSLDS
from hdr_athero.simulation.controller import CONTROLLERS
from hdr_athero.analysis.spectral import (
    spectral_analysis, controllability_report,
)


class TestParameterShapes:
    """Verify all parameter matrices have correct dimensions."""

    def test_tau_shape(self):
        assert TAU.shape == (N_BASINS, N_STATE)

    def test_A_shapes(self):
        for k in range(N_BASINS):
            assert A[k].shape == (N_STATE, N_STATE)

    def test_B_shapes(self):
        for k in range(N_BASINS):
            assert B[k].shape == (N_STATE, N_CONTROL)

    def test_J_shapes(self):
        for k in range(N_BASINS):
            assert J[k].shape == (N_STATE, N_STATE)

    def test_Q_shapes(self):
        for k in range(N_BASINS):
            assert Q[k].shape == (N_STATE, N_STATE)

    def test_cost_matrices(self):
        assert Q_COST.shape == (N_STATE, N_STATE)
        assert R_COST.shape == (N_CONTROL, N_CONTROL)


class TestParameterProperties:
    """Verify biological constraints on parameters."""

    def test_tau_positive(self):
        assert np.all(TAU > 0), "All decay constants must be positive"

    def test_J_zero_diagonal(self):
        for k in range(N_BASINS):
            assert np.allclose(np.diag(J[k]), 0.0), \
                f"Basin {k}: J diagonal must be zero"

    def test_J_nonnegative_offdiag(self):
        """Atherosclerosis coupling is predominantly positive."""
        for k in range(N_BASINS):
            offdiag = J[k] - np.diag(np.diag(J[k]))
            assert np.all(offdiag >= 0), \
                f"Basin {k}: off-diagonal J entries should be non-negative"

    def test_B_therapeutic_sign(self):
        """Most B_k entries should be ≤ 0 (therapeutic direction)."""
        for k in range(N_BASINS):
            assert np.all(B[k] <= 0.0 + 1e-10), \
                f"Basin {k}: B_k entries should be non-positive"

    def test_Q_spd(self):
        """Process noise covariance must be symmetric positive definite."""
        for k in range(N_BASINS):
            eigvals = np.linalg.eigvalsh(Q[k])
            assert np.all(eigvals > 0), f"Basin {k}: Q must be SPD"

    def test_P_trans_rows_valid(self):
        """Transition probabilities: off-diagonal sum ≤ 1."""
        for k in range(N_BASINS):
            off_diag_sum = P_TRANS[k].sum()
            assert off_diag_sum <= 1.0 + 1e-10, \
                f"Basin {k}: transition probs sum to {off_diag_sum}"

    def test_basin1_tau_larger(self):
        """Basin 1 (vulnerable) should have larger τ than Basin 0."""
        for i in range(N_STATE):
            assert TAU[1, i] >= TAU[0, i], \
                f"Axis {i}: Basin 1 τ should be ≥ Basin 0 τ"


class TestStability:
    """Verify ρ(A_k) < 1 for all basins."""

    def test_all_basins_stable(self):
        spec = spectral_analysis()
        for name, r in spec.items():
            assert r["stable"], \
                f"Basin {name}: ρ(A_k) = {r['spectral_radius']:.4f} ≥ 1"

    def test_spectral_ordering(self):
        """ρ(A_2) < ρ(A_0) < ρ(A_1)."""
        spec = spectral_analysis()
        rho = {name: r["spectral_radius"] for name, r in spec.items()}
        assert rho["post_ACS"] < rho["subclinical"], \
            f"Expected ρ(post_ACS) < ρ(subclinical), got {rho}"
        assert rho["subclinical"] < rho["vulnerable"], \
            f"Expected ρ(subclinical) < ρ(vulnerable), got {rho}"

    def test_vulnerable_near_unit_root(self):
        """Basin 1 should have ρ ∈ [0.88, 0.998]."""
        spec = spectral_analysis()
        rho = spec["vulnerable"]["spectral_radius"]
        assert 0.88 <= rho <= 0.998, \
            f"Basin 1 ρ = {rho:.4f}, expected [0.88, 0.998]"


class TestControllability:
    """Verify controllability of the (A_k, B_k) pairs."""

    def test_all_basins_controllable(self):
        report = controllability_report()
        for name, r in report.items():
            assert r["full_rank"], \
                f"Basin {name}: controllability rank = {r['rank']}/{N_STATE}"

    def test_gramian_positive_definite(self):
        report = controllability_report()
        for name, r in report.items():
            assert r["min_eigval"] > 0, \
                f"Basin {name}: Gramian min eigenvalue = {r['min_eigval']}"


class TestSimulation:
    """End-to-end simulation smoke tests."""

    def test_open_loop_runs(self):
        model = AtheroSLDS(seed=42)
        x0 = np.array([1.0, 1.0, 1.0, 1.0])
        episode = model.simulate(x0, 0, CONTROLLERS["open_loop"], T=30)
        assert episode.states.shape == (31, N_STATE)
        assert np.all(np.isfinite(episode.states))
        assert np.all(np.isfinite(episode.costs))

    def test_lqr_basin_runs(self):
        model = AtheroSLDS(seed=42)
        x0 = np.array([1.0, 1.0, 1.0, 1.0])
        episode = model.simulate(x0, 1, CONTROLLERS["lqr_basin"], T=30)
        assert episode.states.shape == (31, N_STATE)
        assert np.all(np.isfinite(episode.states))

    def test_lqr_reduces_cost_vs_open_loop(self):
        """LQR should outperform open-loop in Basin 0.

        Note: Basin 1 comparison is unreliable because the one-sided
        control constraint (u ≥ 0) breaks standard LQR optimality.
        The near-unit-root dynamics amplify the mismatch.  This is a
        known limitation — a proper nonnegative MPC solver is needed
        for Basin 1.  We test on Basin 0 where the LQR is effective.
        """
        costs_ol = []
        costs_lqr = []
        for ep in range(10):
            model_ol = AtheroSLDS(seed=ep)
            model_lqr = AtheroSLDS(seed=ep)
            x0 = np.abs(model_ol.rng.normal(scale=0.8, size=N_STATE))

            ep_ol = model_ol.simulate(x0, 0, CONTROLLERS["open_loop"], T=90)
            ep_lqr = model_lqr.simulate(x0, 0, CONTROLLERS["lqr_basin"], T=90)

            costs_ol.append(ep_ol.mean_cost)
            costs_lqr.append(ep_lqr.mean_cost)

        assert np.mean(costs_lqr) < np.mean(costs_ol), \
            f"LQR ({np.mean(costs_lqr):.3f}) should beat open-loop ({np.mean(costs_ol):.3f})"

    def test_basin_transitions_occur(self):
        """Over a long simulation, basin transitions should occur."""
        model = AtheroSLDS(seed=123)
        x0 = np.zeros(N_STATE)
        episode = model.simulate(x0, 0, CONTROLLERS["open_loop"], T=5000)
        unique_basins = len(set(episode.basins))
        assert unique_basins > 1, "Expected basin transitions over 5000 steps"

    def test_controls_within_bounds(self):
        """All controllers should produce u ∈ [0, u_max]."""
        model = AtheroSLDS(seed=42)
        x0 = np.array([2.0, 2.0, 2.0, 2.0])
        for ctrl_name, ctrl_fn in CONTROLLERS.items():
            episode = model.simulate(x0, 1, ctrl_fn, T=30)
            assert np.all(episode.controls >= -1e-10), \
                f"{ctrl_name}: negative control detected"
            assert np.all(episode.controls <= 1.0 + 1e-10), \
                f"{ctrl_name}: control exceeds bounds"

    def test_recovery_surrogate_amplification(self):
        """Basin 1 recovery surrogate should be much larger than Basin 0."""
        model = AtheroSLDS(seed=0)
        x = np.array([1.0, 1.0, 1.0, 1.0])
        tau_0 = model.recovery_surrogate(x, 0)
        tau_1 = model.recovery_surrogate(x, 1)
        assert tau_1 > 3 * tau_0, \
            f"Basin 1 τ̃ ({tau_1:.2f}) should be >> Basin 0 τ̃ ({tau_0:.2f})"
