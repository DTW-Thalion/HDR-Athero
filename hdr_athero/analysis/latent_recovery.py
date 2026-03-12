"""
PCA / Autoencoder latent space recovery experiment.

Tests whether the hand-defined physiological axes (VI, LD, ED, HS) align
with the dominant directions of variation discovered by data-driven methods
(PCA, linear autoencoder, nonlinear VAE) on synthetic observations.

Run with:
    python -c "from hdr_athero.analysis.latent_recovery import run_experiment; run_experiment()"
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.linalg import orthogonal_procrustes
from scipy.optimize import linear_sum_assignment

from hdr_athero.model.parameters import (
    A, N_STATE, N_BASINS, AXIS_NAMES, BASIN_NAMES,
)
from hdr_athero.model.slds import AtheroSLDS
from hdr_athero.simulation.controller import CONTROLLERS
from hdr_athero.analysis.observation_model import (
    N_OBS,
    build_observation_matrix,
    build_baseline_offset,
    build_noise_covariance,
    build_missingness_rates,
    generate_observations,
    generate_missingness_mask,
)


# ═══════════════════════════════════════════════════════════════════════════
# Stage 1: Synthetic dataset generation
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class SyntheticDataset:
    """Container for a synthetic observation dataset."""
    states: np.ndarray          # (N, 4) ground-truth latent states
    basins: np.ndarray          # (N,) basin labels
    observations: np.ndarray    # (N, 16) complete observations
    mask: np.ndarray            # (N, 16) missingness mask (1=observed)
    C: np.ndarray               # (16, 4) observation matrix
    R: np.ndarray               # (16, 16) noise covariance


def generate_dataset(
    n_episodes: int = 100,
    T: int = 500,
    seed: int = 42,
) -> SyntheticDataset:
    """Generate a synthetic dataset by simulating the SLDS across basins.

    Runs multiple episodes with varied initial conditions, starting basins,
    and controllers (open-loop and LQR) to produce diverse state trajectories.

    Parameters
    ----------
    n_episodes : number of episodes to simulate.
    T : length of each episode in days.
    seed : base random seed.

    Returns
    -------
    SyntheticDataset with ~n_episodes * T observation vectors.
    """
    rng = np.random.default_rng(seed)
    C = build_observation_matrix()
    c = build_baseline_offset()
    R = build_noise_covariance()

    all_states = []
    all_basins = []

    ctrl_keys = ["open_loop", "lqr_basin"]

    for ep in range(n_episodes):
        z0 = ep % N_BASINS
        ctrl_name = ctrl_keys[ep % len(ctrl_keys)]
        ctrl_fn = CONTROLLERS[ctrl_name]

        model = AtheroSLDS(seed=seed + ep * 137)
        x0 = np.abs(rng.normal(scale=0.8, size=N_STATE))

        episode = model.simulate(x0, z0, ctrl_fn, T)

        # Use states[:-1] to match basins[:-1] (T states with known basins)
        all_states.append(episode.states[:-1])
        all_basins.append(episode.basins[:-1])

    states = np.concatenate(all_states, axis=0)  # (N, 4)
    basins = np.concatenate(all_basins, axis=0)   # (N,)

    # Generate observations
    observations = generate_observations(
        states, C=C, c=c, R=R, seed=seed + 999
    )

    # Generate missingness mask
    mask = generate_missingness_mask(
        states.shape[0], seed=seed + 888
    )

    return SyntheticDataset(
        states=states,
        basins=basins,
        observations=observations,
        mask=mask,
        C=C,
        R=R,
    )


# ═══════════════════════════════════════════════════════════════════════════
# Stage 2: Recovery methods
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class RecoveryResult:
    """Result of a latent space recovery method."""
    name: str
    latent: np.ndarray           # (N, 4) recovered latent coordinates
    loadings: np.ndarray         # (16, 4) or equivalent mapping
    explained_variance: np.ndarray | None  # (4,) variance ratios (PCA only)
    extra: dict = field(default_factory=dict)


def fit_pca(
    Y: np.ndarray,
    n_components: int = 4,
    mask: np.ndarray | None = None,
) -> RecoveryResult:
    """Fit PCA on observations using SVD.

    If mask is provided, impute missing entries with column means before fitting.

    Parameters
    ----------
    Y : observations (N, 16).
    n_components : number of components to extract.
    mask : optional missingness mask (N, 16), 1=observed.

    Returns
    -------
    RecoveryResult with PCA loadings and scores.
    """
    Y_work = Y.copy()

    if mask is not None:
        # Impute missing values with column means of observed entries
        col_means = np.zeros(Y.shape[1])
        for j in range(Y.shape[1]):
            observed = mask[:, j].astype(bool)
            if observed.any():
                col_means[j] = Y[observed, j].mean()
        Y_work[~mask.astype(bool)] = np.tile(
            col_means, (Y.shape[0], 1)
        )[~mask.astype(bool)]

    # Center
    mean = Y_work.mean(axis=0)
    Y_centered = Y_work - mean

    # SVD
    U, S, Vt = np.linalg.svd(Y_centered, full_matrices=False)

    # Loadings (columns of V = rows of Vt transposed)
    loadings = Vt[:n_components].T  # (16, n_components)

    # Scores
    scores = Y_centered @ loadings  # (N, n_components)

    # Explained variance ratios
    total_var = np.sum(S ** 2)
    explained = S[:n_components] ** 2 / total_var

    return RecoveryResult(
        name="PCA",
        latent=scores,
        loadings=loadings,
        explained_variance=explained,
    )


def fit_linear_autoencoder(
    Y: np.ndarray,
    n_latent: int = 4,
    n_epochs: int = 200,
    lr: float = 1e-3,
    mask: np.ndarray | None = None,
    seed: int = 0,
) -> RecoveryResult:
    """Fit a linear autoencoder using gradient descent.

    Encoder: W_enc (n_latent x p), Decoder: W_dec (p x n_latent).
    Loss: MSE reconstruction, optionally masked.

    Parameters
    ----------
    Y : observations (N, p).
    n_latent : bottleneck dimension.
    n_epochs : training epochs.
    lr : learning rate.
    mask : optional missingness mask (N, p), 1=observed.
    seed : random seed.

    Returns
    -------
    RecoveryResult
    """
    rng = np.random.default_rng(seed)
    N, p = Y.shape

    # Center
    if mask is not None:
        col_means = np.zeros(p)
        for j in range(p):
            obs = mask[:, j].astype(bool)
            if obs.any():
                col_means[j] = Y[obs, j].mean()
        Y_imp = Y.copy()
        Y_imp[~mask.astype(bool)] = np.tile(
            col_means, (N, 1)
        )[~mask.astype(bool)]
    else:
        Y_imp = Y.copy()

    mean = Y_imp.mean(axis=0)
    Y_c = Y_imp - mean

    # Initialize weights
    W_enc = rng.normal(0, 0.1, size=(n_latent, p)).astype(np.float64)
    W_dec = rng.normal(0, 0.1, size=(p, n_latent)).astype(np.float64)

    # Mini-batch SGD
    batch_size = min(256, N)

    for epoch in range(n_epochs):
        # Shuffle indices
        idx = rng.permutation(N)
        for start in range(0, N, batch_size):
            batch_idx = idx[start:start + batch_size]
            Y_batch = Y_c[batch_idx]  # (B, p)

            # Forward
            Z_batch = Y_batch @ W_enc.T  # (B, n_latent)
            Y_hat = Z_batch @ W_dec.T    # (B, p)
            residual = Y_hat - Y_batch   # (B, p)

            if mask is not None:
                M_batch = mask[batch_idx].astype(np.float64)
                residual = residual * M_batch

            B = len(batch_idx)

            # Gradients
            # dL/dW_dec: (p, n_latent)
            grad_W_dec = (residual.T @ Z_batch) / B
            # dL/dW_enc: chain through decoder
            # dL/dZ = residual @ W_dec  → (B, n_latent)
            dL_dZ = residual @ W_dec
            grad_W_enc = (dL_dZ.T @ Y_batch) / B

            W_dec -= lr * grad_W_dec
            W_enc -= lr * grad_W_enc

    # Final encoding
    latent = Y_c @ W_enc.T

    return RecoveryResult(
        name="Linear AE",
        latent=latent,
        loadings=W_enc.T,  # (p, n_latent) to match PCA convention
        explained_variance=None,
    )


def _relu(x: np.ndarray) -> np.ndarray:
    """Element-wise ReLU."""
    return np.maximum(0, x)


def _relu_grad(x: np.ndarray) -> np.ndarray:
    """Gradient of ReLU."""
    return (x > 0).astype(np.float64)


def fit_vae(
    Y: np.ndarray,
    n_latent: int = 4,
    n_epochs: int = 200,
    lr: float = 1e-3,
    beta: float = 0.1,
    batch_size: int = 256,
    mask: np.ndarray | None = None,
    seed: int = 0,
) -> RecoveryResult:
    """Fit a nonlinear VAE using numpy-based gradient descent.

    Architecture:
        Encoder: p → 32 (ReLU) → 16 (ReLU) → n_latent (mean) + n_latent (log-var)
        Decoder: n_latent → 16 (ReLU) → 32 (ReLU) → p

    Loss: ELBO = MSE_reconstruction + beta * KL_divergence

    Parameters
    ----------
    Y : observations (N, p).
    n_latent : latent dimension.
    n_epochs : training epochs.
    lr : learning rate.
    beta : KL weight (beta-VAE).
    batch_size : mini-batch size.
    mask : optional missingness mask (N, p), 1=observed.
    seed : random seed.

    Returns
    -------
    RecoveryResult
    """
    rng = np.random.default_rng(seed)
    N, p = Y.shape

    # Pre-process: center + impute if masked
    if mask is not None:
        col_means = np.zeros(p)
        for j in range(p):
            obs = mask[:, j].astype(bool)
            if obs.any():
                col_means[j] = Y[obs, j].mean()
        Y_imp = Y.copy()
        Y_imp[~mask.astype(bool)] = np.tile(
            col_means, (N, 1)
        )[~mask.astype(bool)]
    else:
        Y_imp = Y.copy()

    mean = Y_imp.mean(axis=0)
    std = Y_imp.std(axis=0)
    std[std < 1e-8] = 1.0
    Y_c = (Y_imp - mean) / std

    # Xavier initialization helper
    def _init(fan_in: int, fan_out: int) -> np.ndarray:
        scale = np.sqrt(2.0 / (fan_in + fan_out))
        return rng.normal(0, scale, size=(fan_out, fan_in))

    # Encoder weights
    W_e1 = _init(p, 32)       # (32, p)
    b_e1 = np.zeros(32)
    W_e2 = _init(32, 16)      # (16, 32)
    b_e2 = np.zeros(16)
    W_mu = _init(16, n_latent)  # (n_latent, 16)
    b_mu = np.zeros(n_latent)
    W_lv = _init(16, n_latent)  # (n_latent, 16)
    b_lv = np.zeros(n_latent)

    # Decoder weights
    W_d1 = _init(n_latent, 16)  # (16, n_latent)
    b_d1 = np.zeros(16)
    W_d2 = _init(16, 32)        # (32, 16)
    b_d2 = np.zeros(32)
    W_d3 = _init(32, p)         # (p, 32)
    b_d3 = np.zeros(p)

    # Adam state
    params = [W_e1, b_e1, W_e2, b_e2, W_mu, b_mu, W_lv, b_lv,
              W_d1, b_d1, W_d2, b_d2, W_d3, b_d3]
    m_adam = [np.zeros_like(p_) for p_ in params]
    v_adam = [np.zeros_like(p_) for p_ in params]
    beta1, beta2, eps_adam = 0.9, 0.999, 1e-8

    batch_size = min(batch_size, N)
    t_step = 0

    for epoch in range(n_epochs):
        idx = rng.permutation(N)
        for start in range(0, N, batch_size):
            t_step += 1
            bi = idx[start:start + batch_size]
            B = len(bi)
            x_batch = Y_c[bi]  # (B, p)

            # ── Forward pass ──
            # Encoder
            h1_pre = x_batch @ W_e1.T + b_e1  # (B, 32)
            h1 = _relu(h1_pre)
            h2_pre = h1 @ W_e2.T + b_e2       # (B, 16)
            h2 = _relu(h2_pre)

            mu = h2 @ W_mu.T + b_mu            # (B, n_latent)
            log_var = h2 @ W_lv.T + b_lv       # (B, n_latent)
            # Clamp log_var
            log_var = np.clip(log_var, -10, 10)
            sd = np.exp(0.5 * log_var)

            # Reparameterize
            eps = rng.normal(0, 1, size=(B, n_latent))
            z = mu + sd * eps                   # (B, n_latent)

            # Decoder
            d1_pre = z @ W_d1.T + b_d1         # (B, 16)
            d1 = _relu(d1_pre)
            d2_pre = d1 @ W_d2.T + b_d2        # (B, 32)
            d2 = _relu(d2_pre)
            y_hat = d2 @ W_d3.T + b_d3         # (B, p)

            # ── Loss ──
            recon_err = y_hat - x_batch         # (B, p)
            if mask is not None:
                M_batch = mask[bi].astype(np.float64)
                recon_err = recon_err * M_batch
                n_obs = M_batch.sum() + 1e-8
            else:
                n_obs = float(B * p)

            mse = np.sum(recon_err ** 2) / n_obs
            kl = -0.5 * np.sum(1 + log_var - mu ** 2 - np.exp(log_var)) / B

            # ── Backward pass ──
            # d(MSE)/d(y_hat)
            d_yhat = 2.0 * recon_err / n_obs    # (B, p)

            # Decoder backward
            d_d2 = d_yhat @ W_d3                 # (B, 32)
            gW_d3 = d_yhat.T @ d2                # (p, 32)
            gb_d3 = d_yhat.sum(axis=0)

            d_d2 = d_d2 * _relu_grad(d2_pre)
            d_d1 = d_d2 @ W_d2                   # (B, 16)
            gW_d2 = d_d2.T @ d1                   # (32, 16)
            gb_d2 = d_d2.sum(axis=0)

            d_d1 = d_d1 * _relu_grad(d1_pre)
            d_z = d_d1 @ W_d1                     # (B, n_latent)
            gW_d1 = d_d1.T @ z                    # (16, n_latent)
            gb_d1 = d_d1.sum(axis=0)

            # KL gradients w.r.t. mu and log_var
            d_mu_kl = mu / B                      # (B, n_latent)
            d_lv_kl = 0.5 * (np.exp(log_var) - 1) / B

            # d_z → d_mu and d_log_var via reparameterization
            d_mu_total = d_z + beta * d_mu_kl
            d_lv_total = d_z * eps * 0.5 * sd + beta * d_lv_kl

            # Encoder backward from mu/log_var
            d_h2 = d_mu_total @ W_mu + d_lv_total @ W_lv  # (B, 16)
            gW_mu = d_mu_total.T @ h2             # (n_latent, 16)
            gb_mu = d_mu_total.sum(axis=0)
            gW_lv = d_lv_total.T @ h2
            gb_lv = d_lv_total.sum(axis=0)

            d_h2 = d_h2 * _relu_grad(h2_pre)
            d_h1 = d_h2 @ W_e2                    # (B, 32)
            gW_e2 = d_h2.T @ h1                   # (16, 32)
            gb_e2 = d_h2.sum(axis=0)

            d_h1 = d_h1 * _relu_grad(h1_pre)
            gW_e1 = d_h1.T @ x_batch              # (32, p)
            gb_e1 = d_h1.sum(axis=0)

            # ── Adam update ──
            grads = [gW_e1, gb_e1, gW_e2, gb_e2, gW_mu, gb_mu, gW_lv, gb_lv,
                     gW_d1, gb_d1, gW_d2, gb_d2, gW_d3, gb_d3]

            for i, (p_, g) in enumerate(zip(params, grads)):
                m_adam[i] = beta1 * m_adam[i] + (1 - beta1) * g
                v_adam[i] = beta2 * v_adam[i] + (1 - beta2) * g ** 2
                m_hat = m_adam[i] / (1 - beta1 ** t_step)
                v_hat = v_adam[i] / (1 - beta2 ** t_step)
                p_ -= lr * m_hat / (np.sqrt(v_hat) + eps_adam)

    # ── Final encoding (use mean, no sampling) ──
    h1 = _relu(Y_c @ W_e1.T + b_e1)
    h2 = _relu(h1 @ W_e2.T + b_e2)
    latent = h2 @ W_mu.T + b_mu

    # Approximate linear loadings via Jacobian at zero
    h1_0 = _relu(np.zeros(p) @ W_e1.T + b_e1)
    h2_0 = _relu(h1_0 @ W_e2.T + b_e2)
    # For nonlinear models, loadings are the effective first-order mapping
    # Use W_mu @ diag(relu_grad) @ W_e2 @ diag(relu_grad) @ W_e1 as approx
    # But simpler: just use correlation-based alignment
    approx_loadings = W_mu @ W_e2 @ W_e1  # (n_latent, p)

    return RecoveryResult(
        name="VAE",
        latent=latent,
        loadings=approx_loadings.T,  # (p, n_latent)
        explained_variance=None,
        extra={"beta": beta},
    )


# ═══════════════════════════════════════════════════════════════════════════
# Stage 3: Alignment analysis
# ═══════════════════════════════════════════════════════════════════════════

def axis_alignment_matrix(
    recovered: np.ndarray,
    ground_truth: np.ndarray,
) -> np.ndarray:
    """Compute the 4x4 correlation matrix between recovered and ground-truth axes.

    Parameters
    ----------
    recovered : (N, 4) recovered latent coordinates.
    ground_truth : (N, 4) ground-truth latent states.

    Returns
    -------
    corr : (4, 4) correlation matrix where corr[i, j] = correlation
           between recovered dim i and ground-truth dim j.
    """
    n_dims = ground_truth.shape[1]
    corr = np.zeros((n_dims, n_dims))
    for i in range(n_dims):
        for j in range(n_dims):
            r = np.corrcoef(recovered[:, i], ground_truth[:, j])[0, 1]
            corr[i, j] = r if np.isfinite(r) else 0.0
    return corr


def optimal_permutation_correlation(corr: np.ndarray) -> tuple[float, np.ndarray]:
    """Find the optimal permutation that maximises diagonal |correlation|.

    Uses the Hungarian algorithm on |corr|.

    Returns
    -------
    mean_diag : mean |correlation| on the optimally-permuted diagonal.
    perm : permutation array.
    """
    cost = -np.abs(corr)
    row_ind, col_ind = linear_sum_assignment(cost)
    perm = col_ind
    diag_corrs = np.array([np.abs(corr[i, perm[i]]) for i in range(len(perm))])
    return float(diag_corrs.mean()), perm


def procrustes_distance(
    recovered: np.ndarray,
    ground_truth: np.ndarray,
) -> tuple[float, np.ndarray]:
    """Compute Procrustes distance between recovered and ground-truth latent spaces.

    Standardises both matrices, then finds the optimal orthogonal rotation.

    Returns
    -------
    distance : residual Frobenius norm after alignment (normalised).
    R_opt : optimal rotation matrix (4x4).
    """
    # Standardise columns
    def _std(X: np.ndarray) -> np.ndarray:
        m = X.mean(axis=0)
        s = X.std(axis=0)
        s[s < 1e-10] = 1.0
        return (X - m) / s

    Z_rec = _std(recovered)
    Z_gt = _std(ground_truth)

    # Orthogonal Procrustes: find R minimising ||Z_rec @ R - Z_gt||_F
    R_opt, scale = orthogonal_procrustes(Z_rec, Z_gt)
    aligned = Z_rec @ R_opt

    residual = np.linalg.norm(aligned - Z_gt, "fro")
    normaliser = np.linalg.norm(Z_gt, "fro")
    distance = residual / normaliser if normaliser > 1e-10 else residual

    return float(distance), R_opt


def subspace_angles(
    recovered: np.ndarray,
    ground_truth: np.ndarray,
) -> np.ndarray:
    """Compute principal angles between the 4D subspaces.

    Uses SVD of Q_1^T Q_2 where Q_1, Q_2 are orthonormal bases.

    Returns
    -------
    angles : (4,) principal angles in radians.
    """
    def _orth_basis(X: np.ndarray) -> np.ndarray:
        X_c = X - X.mean(axis=0)
        Q, _ = np.linalg.qr(X_c)
        return Q[:, :X.shape[1]]

    Q1 = _orth_basis(recovered)
    Q2 = _orth_basis(ground_truth)

    _, S, _ = np.linalg.svd(Q1.T @ Q2)
    # Clamp to [0, 1] for numerical safety
    S = np.clip(S, 0.0, 1.0)
    angles = np.arccos(S)
    return angles


def fit_var1(X: np.ndarray) -> np.ndarray:
    """Fit a VAR(1) model: X[t+1] = A_hat @ X[t] + noise.

    Parameters
    ----------
    X : (T, n) time series.

    Returns
    -------
    A_hat : (n, n) estimated transition matrix.
    """
    X_lag = X[:-1]   # (T-1, n)
    X_lead = X[1:]   # (T-1, n)
    # OLS: A_hat = (X_lead^T X_lag) (X_lag^T X_lag)^{-1}
    A_hat, _, _, _ = np.linalg.lstsq(X_lag, X_lead, rcond=None)
    return A_hat.T  # lstsq gives (n, n) where A_hat[j, i], transpose to standard form


def dynamics_recovery_error(
    recovered_latent: np.ndarray,
    ground_truth_states: np.ndarray,
    ground_truth_basins: np.ndarray,
    R_procrustes: np.ndarray,
) -> dict[str, float]:
    """Compute dynamics recovery error per basin.

    Aligns recovered latent space using Procrustes rotation, fits VAR(1)
    in the aligned coordinates, and compares with true A_k.

    Returns
    -------
    errors : dict mapping basin name to relative Frobenius error.
    """
    # Align recovered latent to ground-truth coordinates
    rec_std = recovered_latent - recovered_latent.mean(axis=0)
    s = recovered_latent.std(axis=0)
    s[s < 1e-10] = 1.0
    rec_std = rec_std / s
    aligned = rec_std @ R_procrustes

    # Scale to match ground-truth variance
    gt_std = ground_truth_states.std(axis=0)
    al_std = aligned.std(axis=0)
    scale = gt_std / (al_std + 1e-10)
    aligned = aligned * scale

    errors = {}
    for k in range(N_BASINS):
        basin_mask = ground_truth_basins == k
        if basin_mask.sum() < 20:
            errors[BASIN_NAMES[k]] = float("nan")
            continue

        # Find contiguous segments in this basin
        indices = np.where(basin_mask)[0]
        # Use only consecutive pairs
        consecutive = np.where(np.diff(indices) == 1)[0]
        if len(consecutive) < 10:
            errors[BASIN_NAMES[k]] = float("nan")
            continue

        X_seg = aligned[indices[consecutive]]
        X_seg_next = aligned[indices[consecutive] + 1]

        # Fit VAR(1) on this basin segment
        A_hat, _, _, _ = np.linalg.lstsq(X_seg, X_seg_next, rcond=None)
        A_hat = A_hat.T

        A_true = A[k]
        rel_error = np.linalg.norm(A_hat - A_true, "fro") / np.linalg.norm(A_true, "fro")
        errors[BASIN_NAMES[k]] = float(rel_error)

    return errors


def cross_loading_detection(corr: np.ndarray, perm: np.ndarray) -> dict[str, Any]:
    """Check whether PCA merges the VI-LD axes.

    Looks at the correlation matrix after optimal permutation to detect
    whether VI and LD load onto the same recovered component.

    Returns
    -------
    info : dict with flags about axis merging.
    """
    # After permutation: corr[i, perm[i]] should be dominant
    # Check if VI (axis 0) and LD (axis 1) map to the same recovered dim
    vi_component = perm[0] if 0 in perm else -1
    ld_component = perm[1] if 1 in perm else -1

    # Actually: perm[i] = which ground-truth axis maps to recovered dim i
    # We want: does any recovered dim have strong loading on both VI and LD?
    merged = False
    for i in range(corr.shape[0]):
        vi_load = abs(corr[i, 0])
        ld_load = abs(corr[i, 1])
        if vi_load > 0.4 and ld_load > 0.4:
            merged = True
            break

    # Separation quality: how well are VI and LD on different components?
    # Best case: each on its own component with high correlation
    vi_best = int(np.argmax(np.abs(corr[:, 0])))
    ld_best = int(np.argmax(np.abs(corr[:, 1])))
    separated = vi_best != ld_best

    return {
        "vi_ld_merged": merged,
        "vi_ld_separated": separated,
        "vi_best_component": vi_best,
        "ld_best_component": ld_best,
        "vi_max_corr": float(np.max(np.abs(corr[:, 0]))),
        "ld_max_corr": float(np.max(np.abs(corr[:, 1]))),
    }


# ═══════════════════════════════════════════════════════════════════════════
# Full experiment
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class MethodResult:
    """Full results for one method × one data condition."""
    method_name: str
    data_condition: str  # "complete" or "degraded"
    alignment_corr: np.ndarray       # (4, 4) correlation matrix
    mean_diag_corr: float
    perm: np.ndarray
    procrustes_dist: float
    R_procrustes: np.ndarray
    subspace_angles: np.ndarray
    per_basin_alignment: dict[str, float]
    dynamics_errors: dict[str, float]
    cross_loading: dict[str, Any]
    explained_variance: np.ndarray | None


def _analyse_recovery(
    result: RecoveryResult,
    dataset: SyntheticDataset,
    data_condition: str,
) -> MethodResult:
    """Run all alignment analyses on a recovery result."""
    gt = dataset.states

    # Axis alignment
    corr = axis_alignment_matrix(result.latent, gt)
    mean_diag, perm = optimal_permutation_correlation(corr)

    # Procrustes
    pdist, R_proc = procrustes_distance(result.latent, gt)

    # Subspace angles
    angles = subspace_angles(result.latent, gt)

    # Per-basin alignment
    per_basin = {}
    for k in range(N_BASINS):
        mask_k = dataset.basins == k
        if mask_k.sum() < 10:
            per_basin[BASIN_NAMES[k]] = float("nan")
            continue
        corr_k = axis_alignment_matrix(result.latent[mask_k], gt[mask_k])
        diag_k, _ = optimal_permutation_correlation(corr_k)
        per_basin[BASIN_NAMES[k]] = diag_k

    # Dynamics recovery
    dyn_errors = dynamics_recovery_error(
        result.latent, gt, dataset.basins, R_proc
    )

    # Cross-loading detection
    xl = cross_loading_detection(corr, perm)

    return MethodResult(
        method_name=result.name,
        data_condition=data_condition,
        alignment_corr=corr,
        mean_diag_corr=mean_diag,
        perm=perm,
        procrustes_dist=pdist,
        R_procrustes=R_proc,
        subspace_angles=angles,
        per_basin_alignment=per_basin,
        dynamics_errors=dyn_errors,
        cross_loading=xl,
        explained_variance=result.explained_variance,
    )


def run_experiment(
    n_episodes: int = 100,
    T: int = 500,
    seed: int = 42,
    vae_epochs: int = 200,
    quiet: bool = False,
) -> dict[str, Any]:
    """Run the full PCA / AE / VAE latent space recovery experiment.

    Parameters
    ----------
    n_episodes : episodes for dataset generation.
    T : steps per episode.
    seed : random seed.
    vae_epochs : training epochs for VAE.
    quiet : suppress printed report.

    Returns
    -------
    results : dict with all metrics and method results.
    """
    # ── Generate data ──
    if not quiet:
        print("\n" + "=" * 72)
        print("  LATENT SPACE RECOVERY EXPERIMENT")
        print("  Do PCA / autoencoder axes align with the domain-defined")
        print("  VI, LD, ED, HS axes?")
        print("=" * 72)
        print(f"\n  Generating synthetic dataset ({n_episodes} episodes x "
              f"{T} days)...")

    dataset = generate_dataset(n_episodes=n_episodes, T=T, seed=seed)
    N_total = dataset.states.shape[0]

    if not quiet:
        print(f"  Dataset: {N_total} observation vectors, "
              f"p={N_OBS} biomarkers, n={N_STATE} latent dims")
        obs_rate = dataset.mask.mean(axis=0)
        print(f"  Overall observation rate: {obs_rate.mean():.1%}")

    # ── Fit methods on complete data ──
    if not quiet:
        print("\n  Fitting methods on COMPLETE data...")

    pca_complete = fit_pca(dataset.observations, n_components=4)
    ae_complete = fit_linear_autoencoder(
        dataset.observations, n_latent=4, n_epochs=200, lr=1e-3, seed=seed + 1
    )
    vae_complete = fit_vae(
        dataset.observations, n_latent=4, n_epochs=vae_epochs,
        lr=1e-3, beta=0.1, seed=seed + 2
    )

    # ── Fit methods on degraded (missingness) data ──
    if not quiet:
        print("  Fitting methods on DEGRADED data (realistic missingness)...")

    pca_degraded = fit_pca(
        dataset.observations, n_components=4, mask=dataset.mask
    )
    ae_degraded = fit_linear_autoencoder(
        dataset.observations, n_latent=4, n_epochs=200, lr=1e-3,
        mask=dataset.mask, seed=seed + 3
    )
    vae_degraded = fit_vae(
        dataset.observations, n_latent=4, n_epochs=vae_epochs,
        lr=1e-3, beta=0.1, mask=dataset.mask, seed=seed + 4
    )

    # ── Analyse all 6 fits ──
    if not quiet:
        print("  Analysing alignment...")

    results_list = []
    for result, cond in [
        (pca_complete, "complete"), (ae_complete, "complete"),
        (vae_complete, "complete"),
        (pca_degraded, "degraded"), (ae_degraded, "degraded"),
        (vae_degraded, "degraded"),
    ]:
        mr = _analyse_recovery(result, dataset, cond)
        results_list.append(mr)

    # ── Build results dict ──
    results = {
        "dataset_size": N_total,
        "methods": {
            f"{mr.method_name}_{mr.data_condition}": mr
            for mr in results_list
        },
    }

    # ── Print report ──
    if not quiet:
        _print_report(results)

    return results


def _print_report(results: dict[str, Any]) -> None:
    """Print a formatted experiment report."""
    methods = results["methods"]

    # ── Main comparison table ──
    print("\n" + "=" * 72)
    print("  RESULTS: Latent Space Alignment")
    print("=" * 72)

    print(f"\n  {'Method':<16s} {'Data':<10s} {'MeanDiag|r|':>12s} "
          f"{'Procrustes':>10s} {'MaxAngle':>10s} {'ExplVar':>10s}")
    print("  " + "-" * 68)

    for key in ["PCA_complete", "Linear AE_complete", "VAE_complete",
                "PCA_degraded", "Linear AE_degraded", "VAE_degraded"]:
        mr = methods[key]
        ev_str = (f"{mr.explained_variance.sum():.3f}"
                  if mr.explained_variance is not None else "N/A")
        max_angle = float(np.max(mr.subspace_angles))
        print(f"  {mr.method_name:<16s} {mr.data_condition:<10s} "
              f"{mr.mean_diag_corr:12.3f} "
              f"{mr.procrustes_dist:10.4f} "
              f"{np.degrees(max_angle):9.1f}deg "
              f"{ev_str:>10s}")

    # ── Per-basin breakdown ──
    print("\n  Per-basin alignment (mean |diagonal correlation|):")
    print(f"  {'Method':<16s} {'Data':<10s}", end="")
    for bn in BASIN_NAMES:
        print(f"  {bn:>14s}", end="")
    print()
    print("  " + "-" * 68)

    for key in ["PCA_complete", "Linear AE_complete", "VAE_complete",
                "PCA_degraded", "Linear AE_degraded", "VAE_degraded"]:
        mr = methods[key]
        print(f"  {mr.method_name:<16s} {mr.data_condition:<10s}", end="")
        for bn in BASIN_NAMES:
            val = mr.per_basin_alignment[bn]
            if np.isfinite(val):
                print(f"  {val:14.3f}", end="")
            else:
                print(f"  {'N/A':>14s}", end="")
        print()

    # ── Dynamics recovery ──
    print("\n  Dynamics recovery (relative Frobenius error ||A_hat - A_k||/||A_k||):")
    print(f"  {'Method':<16s} {'Data':<10s}", end="")
    for bn in BASIN_NAMES:
        print(f"  {bn:>14s}", end="")
    print()
    print("  " + "-" * 68)

    for key in ["PCA_complete", "VAE_complete",
                "PCA_degraded", "VAE_degraded"]:
        mr = methods[key]
        print(f"  {mr.method_name:<16s} {mr.data_condition:<10s}", end="")
        for bn in BASIN_NAMES:
            val = mr.dynamics_errors.get(bn, float("nan"))
            if np.isfinite(val):
                print(f"  {val:14.3f}", end="")
            else:
                print(f"  {'N/A':>14s}", end="")
        print()

    # ── Cross-loading detection ──
    print("\n  Cross-loading detection (VI-LD separation):")
    for key in ["PCA_complete", "VAE_complete"]:
        mr = methods[key]
        xl = mr.cross_loading
        status = "MERGED" if xl["vi_ld_merged"] else "SEPARATED"
        print(f"  {mr.method_name:<16s} [{mr.data_condition}]: "
              f"VI-LD {status}  "
              f"(VI max|r|={xl['vi_max_corr']:.3f}, "
              f"LD max|r|={xl['ld_max_corr']:.3f})")

    # ── PCA explained variance ──
    mr_pca = methods["PCA_complete"]
    if mr_pca.explained_variance is not None:
        ev = mr_pca.explained_variance
        print(f"\n  PCA explained variance ratios: "
              f"{', '.join(f'{v:.3f}' for v in ev)}")
        print(f"  Total (4 components): {ev.sum():.3f}")

    # ── Interpretation ──
    print("\n" + "=" * 72)
    print("  INTERPRETATION FOR HDR")
    print("=" * 72)

    pca_c = methods["PCA_complete"]
    pca_d = methods["PCA_degraded"]
    vae_c = methods["VAE_complete"]

    # Flag (a): does PCA merge VI-LD?
    merged = pca_c.cross_loading["vi_ld_merged"]
    if merged:
        print("\n  (a) PCA MERGES the VI-LD axes. The inflammation-lipid")
        print("      coupling (J matrix) is strong enough that PCA cannot")
        print("      separate them. This supports the domain-defined axis")
        print("      decomposition, which imposes separation based on")
        print("      mechanistic knowledge.")
    else:
        print("\n  (a) PCA SEPARATES the VI-LD axes. The cross-coupling in")
        print("      J is not strong enough to merge them in the first")
        print("      principal components. Domain-defined and data-driven")
        print("      axes are compatible.")

    # Flag (b): does VAE achieve better per-basin separation?
    vae_basin_mean = np.nanmean(list(vae_c.per_basin_alignment.values()))
    pca_basin_mean = np.nanmean(list(pca_c.per_basin_alignment.values()))
    if vae_basin_mean > pca_basin_mean + 0.02:
        print(f"\n  (b) VAE achieves BETTER per-basin separation than PCA")
        print(f"      (VAE: {vae_basin_mean:.3f} vs PCA: {pca_basin_mean:.3f}).")
        print(f"      Nonlinear disentanglement captures basin-specific structure.")
    else:
        print(f"\n  (b) VAE does NOT significantly outperform PCA on per-basin")
        print(f"      separation (VAE: {vae_basin_mean:.3f} vs PCA: {pca_basin_mean:.3f}).")
        print(f"      The linear subspace is sufficient for this model.")

    # Flag (c): missingness degradation
    pca_drop = pca_c.mean_diag_corr - pca_d.mean_diag_corr
    if pca_d.mean_diag_corr < 0.5 * pca_c.mean_diag_corr:
        print(f"\n  (c) Missingness CATASTROPHICALLY degrades recovery")
        print(f"      (complete: {pca_c.mean_diag_corr:.3f} -> "
              f"degraded: {pca_d.mean_diag_corr:.3f}).")
        print(f"      Data-driven latent space methods require high-frequency")
        print(f"      observation infrastructure to be viable.")
    elif pca_drop > 0.1:
        print(f"\n  (c) Missingness MODERATELY degrades recovery")
        print(f"      (complete: {pca_c.mean_diag_corr:.3f} -> "
              f"degraded: {pca_d.mean_diag_corr:.3f}, "
              f"drop: {pca_drop:.3f}).")
    else:
        print(f"\n  (c) Missingness has MINIMAL impact on recovery")
        print(f"      (complete: {pca_c.mean_diag_corr:.3f} -> "
              f"degraded: {pca_d.mean_diag_corr:.3f}).")

    print("\n  Disclaimer: All results are from synthetic data generated by")
    print("  the SLDS model itself. External validation would require real")
    print("  biomarker panel data.\n")
