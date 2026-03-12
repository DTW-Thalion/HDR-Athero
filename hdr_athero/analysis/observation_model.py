"""
Observation model for the 4-axis atherosclerosis SLDS.

Maps the 4D latent state (VI, LD, ED, HS) to a 16D observed biomarker
panel, with realistic cross-loadings, noise levels, and missingness
patterns reflecting clinical measurement schedules.

Observation model:
    y_t = C @ x_t + c + v_t,   v_t ~ N(0, R)
"""
from __future__ import annotations

import numpy as np

from hdr_athero.model.parameters import N_STATE

# ─── Dimensions ──────────────────────────────────────────────────────────
N_OBS = 16  # 4 biomarkers per axis

# ─── Biomarker labels ────────────────────────────────────────────────────
BIOMARKER_NAMES = [
    # Axis 0 (VI): inflammatory
    "hsCRP", "IL-6", "VCAM-1", "GlycA",
    # Axis 1 (LD): lipid
    "LDL-C", "ApoB", "ox-LDL", "triglycerides",
    # Axis 2 (ED): endothelial
    "FMD", "ADMA", "endothelin-1", "NO_metabolites",
    # Axis 3 (HS): haemodynamic
    "SBP", "pulse_pressure", "PWV", "HRV",
]


def build_observation_matrix() -> np.ndarray:
    """Build the 16x4 observation (factor loadings) matrix C.

    Each row corresponds to a biomarker; each column to a latent axis.
    Primary loadings are strong (~1.0-1.5); cross-loadings are weaker
    (~0.1-0.4) and reflect known biological couplings.

    Returns
    -------
    C : ndarray of shape (16, 4)
    """
    C = np.zeros((N_OBS, N_STATE), dtype=np.float64)

    # ── Axis 0 (VI) → inflammatory biomarkers ──
    # y0: hsCRP — primary VI, small cross-loading LD
    C[0, 0] = 1.2
    C[0, 1] = 0.15
    # y1: IL-6 — primary VI
    C[1, 0] = 1.0
    # y2: VCAM-1 — primary VI, small cross-loading ED
    C[2, 0] = 0.9
    C[2, 2] = 0.20
    # y3: GlycA — primary VI, small metabolic cross-loading (mapped to LD)
    C[3, 0] = 0.8
    C[3, 1] = 0.10

    # ── Axis 1 (LD) → lipid biomarkers ──
    # y4: LDL-C — strong primary LD
    C[4, 1] = 1.4
    # y5: ApoB — strong primary LD
    C[5, 1] = 1.3
    # y6: ox-LDL — primary LD, moderate cross-loading VI
    C[6, 1] = 1.0
    C[6, 0] = 0.35
    # y7: triglycerides — primary LD, small cross-loading HS
    C[7, 1] = 0.9
    C[7, 3] = 0.15

    # ── Axis 2 (ED) → endothelial biomarkers ──
    # y8: FMD (sign-inverted) — primary ED
    C[8, 2] = 1.1
    # y9: ADMA — primary ED, small cross-loading HS
    C[9, 2] = 1.0
    C[9, 3] = 0.20
    # y10: endothelin-1 — primary ED, cross-loading HS
    C[10, 2] = 0.9
    C[10, 3] = 0.30
    # y11: NO metabolites (sign-inverted) — primary ED
    C[11, 2] = 1.0

    # ── Axis 3 (HS) → haemodynamic biomarkers ──
    # y12: SBP — strong primary HS
    C[12, 3] = 1.5
    # y13: pulse pressure — primary HS
    C[13, 3] = 1.2
    # y14: PWV — primary HS, cross-loading ED and VI
    C[14, 3] = 1.0
    C[14, 2] = 0.25
    C[14, 0] = 0.20
    # y15: HRV (sign-inverted) — primary HS, cross-loading VI
    C[15, 3] = 0.9
    C[15, 0] = 0.15

    return C


def build_baseline_offset() -> np.ndarray:
    """Build the baseline offset vector c (R^16).

    These represent mean biomarker values when the latent state is at
    the normative reference (x=0).  Values are in standardised units
    so set to zero for simplicity.

    Returns
    -------
    c : ndarray of shape (16,)
    """
    return np.zeros(N_OBS, dtype=np.float64)


def build_noise_covariance() -> np.ndarray:
    """Build the diagonal observation noise covariance R (16x16).

    R_ii reflects known test-retest variability (CV) for each biomarker.
    Higher CV → larger R_ii.

    Returns
    -------
    R : ndarray of shape (16, 16), diagonal
    """
    # Noise standard deviations (in standardised units)
    # Calibrated so SNR = signal_loading / noise_sd ≈ 3-5 for primaries
    noise_sd = np.array([
        # Inflammatory (hsCRP ~20% CV, IL-6 ~25%, VCAM-1 ~15%, GlycA ~10%)
        0.35, 0.30, 0.25, 0.20,
        # Lipid (LDL-C ~8%, ApoB ~8%, ox-LDL ~15%, TG ~15%)
        0.30, 0.28, 0.30, 0.25,
        # Endothelial (FMD ~18%, ADMA ~12%, ET-1 ~15%, NO ~20%)
        0.30, 0.25, 0.25, 0.30,
        # Haemodynamic (SBP ~5%, PP ~8%, PWV ~10%, HRV ~15%)
        0.30, 0.25, 0.25, 0.25,
    ], dtype=np.float64)

    return np.diag(noise_sd ** 2)


def build_missingness_rates() -> np.ndarray:
    """Build the per-channel daily missingness probability vector.

    Reflects realistic clinical measurement schedules at daily resolution:
    - Wearable channels (SBP, HRV): ~2% missing
    - Monthly labs (hsCRP, IL-6): ~97% missing (≈ blood draw every ~30 days)
    - Quarterly labs (LDL-C, ApoB, ox-LDL, TG): ~99% missing
    - Semi-annual specialist (FMD, ADMA, ET-1, NO): ~99.5% missing

    Returns
    -------
    miss_prob : ndarray of shape (16,)
        Probability that each channel is missing on any given day.
    """
    miss_prob = np.array([
        # Inflammatory — monthly labs
        0.967, 0.967, 0.967, 0.967,
        # Lipid — quarterly labs
        0.989, 0.989, 0.989, 0.989,
        # Endothelial — semi-annual specialist
        0.995, 0.995, 0.995, 0.995,
        # Haemodynamic — wearable (SBP, HRV) and quarterly (PP, PWV)
        0.02, 0.989, 0.989, 0.02,
    ], dtype=np.float64)

    return miss_prob


def generate_observations(
    states: np.ndarray,
    *,
    C: np.ndarray | None = None,
    c: np.ndarray | None = None,
    R: np.ndarray | None = None,
    seed: int = 0,
) -> np.ndarray:
    """Generate complete (no missingness) observations from latent states.

    Parameters
    ----------
    states : ndarray of shape (T, 4)
        Latent state trajectory.
    C : observation matrix (16x4), default built from build_observation_matrix.
    c : baseline offset (16,), default zeros.
    R : noise covariance (16x16), default from build_noise_covariance.
    seed : random seed.

    Returns
    -------
    Y : ndarray of shape (T, 16)
        Complete observation panel.
    """
    rng = np.random.default_rng(seed)

    if C is None:
        C = build_observation_matrix()
    if c is None:
        c = build_baseline_offset()
    if R is None:
        R = build_noise_covariance()

    T = states.shape[0]
    noise_sd = np.sqrt(np.diag(R))
    noise = rng.normal(0.0, 1.0, size=(T, N_OBS)) * noise_sd[None, :]

    Y = states @ C.T + c[None, :] + noise
    return Y


def generate_missingness_mask(
    T: int,
    *,
    miss_prob: np.ndarray | None = None,
    seed: int = 0,
) -> np.ndarray:
    """Generate a missingness mask M ∈ {0,1}^{T x 16}.

    M[t, j] = 1 means channel j is OBSERVED at time t.
    M[t, j] = 0 means channel j is MISSING at time t.

    Parameters
    ----------
    T : number of time steps.
    miss_prob : per-channel missingness probability (16,).
    seed : random seed.

    Returns
    -------
    M : ndarray of shape (T, 16), dtype bool
    """
    rng = np.random.default_rng(seed)

    if miss_prob is None:
        miss_prob = build_missingness_rates()

    # M[t, j] = 1 if observed (i.e., random draw > missingness prob)
    U = rng.random((T, N_OBS))
    M = U >= miss_prob[None, :]

    return M
