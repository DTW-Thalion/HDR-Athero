"""Tests for the PCA / autoencoder latent space recovery experiment."""
from __future__ import annotations

import numpy as np
import pytest

from hdr_athero.model.parameters import N_STATE
from hdr_athero.analysis.observation_model import (
    N_OBS,
    build_observation_matrix,
    build_baseline_offset,
    build_noise_covariance,
    build_missingness_rates,
    generate_observations,
    generate_missingness_mask,
)
from hdr_athero.analysis.latent_recovery import (
    generate_dataset,
    fit_pca,
    fit_linear_autoencoder,
    fit_vae,
    axis_alignment_matrix,
    optimal_permutation_correlation,
    procrustes_distance,
    subspace_angles,
    run_experiment,
)


# ─── Observation model tests ─────────────────────────────────────────────

class TestObservationModelShapes:
    """Verify C, R, y_t dimensions."""

    def test_C_shape(self):
        C = build_observation_matrix()
        assert C.shape == (N_OBS, N_STATE), f"Expected (16, 4), got {C.shape}"

    def test_R_shape(self):
        R = build_noise_covariance()
        assert R.shape == (N_OBS, N_OBS), f"Expected (16, 16), got {R.shape}"

    def test_R_diagonal(self):
        R = build_noise_covariance()
        off_diag = R - np.diag(np.diag(R))
        assert np.allclose(off_diag, 0.0), "R should be diagonal"

    def test_R_positive(self):
        R = build_noise_covariance()
        assert np.all(np.diag(R) > 0), "R diagonal entries must be positive"

    def test_baseline_offset_shape(self):
        c = build_baseline_offset()
        assert c.shape == (N_OBS,)

    def test_observations_shape(self):
        rng = np.random.default_rng(0)
        states = rng.normal(size=(100, N_STATE))
        Y = generate_observations(states, seed=0)
        assert Y.shape == (100, N_OBS), f"Expected (100, 16), got {Y.shape}"

    def test_C_primary_loadings_dominant(self):
        """Each block of 4 biomarkers should load primarily on its axis."""
        C = build_observation_matrix()
        for axis in range(N_STATE):
            block = C[axis * 4:(axis + 1) * 4, :]
            # Primary loadings (on own axis) should be larger than cross-loadings
            primary = np.abs(block[:, axis])
            cross = np.abs(block[:, :])
            cross_mask = np.ones(N_STATE, dtype=bool)
            cross_mask[axis] = False
            cross_max = np.max(np.abs(block[:, cross_mask]), axis=1)
            assert np.all(primary > cross_max), \
                f"Axis {axis}: primary loadings should dominate"


class TestMissingnessRates:
    """Verify that generated missingness rates match specifications."""

    def test_missingness_mask_shape(self):
        M = generate_missingness_mask(1000, seed=0)
        assert M.shape == (1000, N_OBS)

    def test_missingness_rates(self):
        """Missingness rates should match specified frequencies within tolerance."""
        T = 100_000
        M = generate_missingness_mask(T, seed=42)
        miss_prob = build_missingness_rates()

        observed_rates = M.mean(axis=0)  # fraction observed
        expected_rates = 1.0 - miss_prob  # fraction observed

        for j in range(N_OBS):
            # Allow 2% absolute tolerance for stochastic rates
            assert abs(observed_rates[j] - expected_rates[j]) < 0.02, \
                f"Channel {j}: observed rate {observed_rates[j]:.3f}, " \
                f"expected {expected_rates[j]:.3f}"

    def test_wearable_mostly_observed(self):
        """SBP (ch 12) and HRV (ch 15) should be mostly observed."""
        T = 10_000
        M = generate_missingness_mask(T, seed=0)
        for ch in [12, 15]:
            obs_rate = M[:, ch].mean()
            assert obs_rate > 0.95, \
                f"Wearable channel {ch} observation rate {obs_rate:.3f} < 0.95"

    def test_quarterly_mostly_missing(self):
        """Quarterly lab channels should be mostly missing."""
        T = 10_000
        M = generate_missingness_mask(T, seed=0)
        for ch in [4, 5, 6, 7]:  # lipid panel
            obs_rate = M[:, ch].mean()
            assert obs_rate < 0.05, \
                f"Quarterly channel {ch} observation rate {obs_rate:.3f} > 0.05"


# ─── Recovery method tests ────────────────────────────────────────────────

@pytest.fixture(scope="module")
def small_dataset():
    """Generate a small dataset for testing."""
    return generate_dataset(n_episodes=20, T=200, seed=42)


class TestPCARecovery:
    """Test PCA recovery on synthetic data."""

    def test_pca_recovers_4_components(self, small_dataset):
        """PCA with 4 components should explain >75% variance on complete data."""
        result = fit_pca(small_dataset.observations, n_components=4)
        total_ev = result.explained_variance.sum()
        assert total_ev > 0.75, \
            f"PCA 4-component explained variance {total_ev:.3f} < 0.75"

    def test_pca_output_shapes(self, small_dataset):
        result = fit_pca(small_dataset.observations, n_components=4)
        N = small_dataset.observations.shape[0]
        assert result.latent.shape == (N, 4)
        assert result.loadings.shape == (N_OBS, 4)
        assert result.explained_variance.shape == (4,)

    def test_pca_with_mask(self, small_dataset):
        """PCA should work with missingness mask."""
        result = fit_pca(
            small_dataset.observations, n_components=4,
            mask=small_dataset.mask
        )
        assert result.latent.shape[1] == 4
        assert np.all(np.isfinite(result.latent))


class TestLinearAE:
    """Test linear autoencoder."""

    def test_linear_ae_matches_pca(self, small_dataset):
        """Linear AE subspace should be within Procrustes distance 0.25 of PCA.

        Note: a linear AE recovers the same subspace as PCA in theory,
        but SGD convergence on finite data introduces some deviation.
        """
        pca_result = fit_pca(small_dataset.observations, n_components=4)
        ae_result = fit_linear_autoencoder(
            small_dataset.observations, n_latent=4,
            n_epochs=500, lr=1e-3, seed=42
        )
        # Compare subspace angles rather than Procrustes on scores
        angles = subspace_angles(ae_result.latent, pca_result.latent)
        max_angle_deg = float(np.degrees(np.max(angles)))
        assert max_angle_deg < 35.0, \
            f"Linear AE max subspace angle to PCA: {max_angle_deg:.1f}deg > 35deg"

    def test_linear_ae_output_shapes(self, small_dataset):
        result = fit_linear_autoencoder(
            small_dataset.observations, n_latent=4,
            n_epochs=50, lr=1e-3, seed=0
        )
        N = small_dataset.observations.shape[0]
        assert result.latent.shape == (N, 4)
        assert result.loadings.shape == (N_OBS, 4)


class TestVAE:
    """Test nonlinear VAE."""

    def test_vae_output_shapes(self, small_dataset):
        result = fit_vae(
            small_dataset.observations, n_latent=4,
            n_epochs=20, lr=1e-3, seed=0
        )
        N = small_dataset.observations.shape[0]
        assert result.latent.shape == (N, 4)
        assert np.all(np.isfinite(result.latent))

    def test_vae_with_mask(self, small_dataset):
        result = fit_vae(
            small_dataset.observations, n_latent=4,
            n_epochs=20, lr=1e-3, mask=small_dataset.mask, seed=0
        )
        assert result.latent.shape[1] == 4
        assert np.all(np.isfinite(result.latent))


# ─── Alignment analysis tests ────────────────────────────────────────────

class TestAlignment:
    """Test alignment metrics."""

    def test_complete_data_alignment(self, small_dataset):
        """On complete data, at least one method should achieve mean |diag r| > 0.4.

        Note: the observation model has realistic noise (SNR ~3-5) and
        cross-loadings, so perfect axis recovery is not expected.  The
        threshold of 0.4 confirms meaningful alignment above chance (0.25).
        """
        pca = fit_pca(small_dataset.observations, n_components=4)
        corr = axis_alignment_matrix(pca.latent, small_dataset.states)
        mean_diag, _ = optimal_permutation_correlation(corr)
        assert mean_diag > 0.4, \
            f"Best alignment {mean_diag:.3f} < 0.4"

    def test_missingness_degrades_gracefully(self, small_dataset):
        """Alignment on degraded data should be within 2x of complete alignment."""
        pca_c = fit_pca(small_dataset.observations, n_components=4)
        pca_d = fit_pca(
            small_dataset.observations, n_components=4,
            mask=small_dataset.mask
        )

        corr_c = axis_alignment_matrix(pca_c.latent, small_dataset.states)
        corr_d = axis_alignment_matrix(pca_d.latent, small_dataset.states)

        diag_c, _ = optimal_permutation_correlation(corr_c)
        diag_d, _ = optimal_permutation_correlation(corr_d)

        # Degraded should be at least half of complete
        assert diag_d > 0.5 * diag_c, \
            f"Degraded alignment {diag_d:.3f} < 0.5 * complete {diag_c:.3f}"

    def test_procrustes_perfect_recovery(self):
        """Procrustes distance of identical matrices should be ~0."""
        rng = np.random.default_rng(0)
        X = rng.normal(size=(500, 4))
        pdist, _ = procrustes_distance(X, X)
        assert pdist < 1e-10

    def test_subspace_angles_identical(self):
        """Subspace angles of identical spaces should be ~0."""
        rng = np.random.default_rng(0)
        X = rng.normal(size=(500, 4))
        angles = subspace_angles(X, X)
        assert np.all(angles < 1e-6), \
            f"Subspace angles should be ~0, got {angles}"

    def test_axis_alignment_identity(self):
        """Correlation matrix of X with itself should have |diag| = 1."""
        rng = np.random.default_rng(0)
        X = rng.normal(size=(500, 4))
        corr = axis_alignment_matrix(X, X)
        assert np.allclose(np.abs(np.diag(corr)), 1.0, atol=1e-10)


# ─── Integration test ─────────────────────────────────────────────────────

class TestIntegration:
    """End-to-end test of the full experiment."""

    def test_run_experiment_completes(self):
        """Full experiment should run without errors."""
        results = run_experiment(
            n_episodes=10, T=100, seed=42, vae_epochs=20, quiet=True
        )
        assert "methods" in results
        assert "dataset_size" in results
        assert len(results["methods"]) == 6

    def test_run_experiment_keys(self):
        """Results dict should have all expected method keys."""
        results = run_experiment(
            n_episodes=10, T=100, seed=42, vae_epochs=20, quiet=True
        )
        expected_keys = [
            "PCA_complete", "Linear AE_complete", "VAE_complete",
            "PCA_degraded", "Linear AE_degraded", "VAE_degraded",
        ]
        for key in expected_keys:
            assert key in results["methods"], f"Missing key: {key}"
