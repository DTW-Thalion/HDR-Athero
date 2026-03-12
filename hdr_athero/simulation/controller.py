"""
Controllers for the atherosclerosis SLDS.

Implements:
  - Open-loop (no intervention) baseline
  - Static guideline-directed therapy (constant u)
  - Per-basin LQR (Mode A approximation)
  - Pooled LQR (single gain, no basin awareness)
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import solve_discrete_are

from hdr_athero.model.parameters import (
    A, B, N_STATE, N_CONTROL, N_BASINS,
    Q_COST, R_COST, LAMBDA_U, U_MAX_BASIN,
    BURDEN_BUDGET,
)


def open_loop(x: np.ndarray, z: int, t: int) -> np.ndarray:
    """No intervention — natural disease progression."""
    return np.zeros(N_CONTROL)


def static_guideline(x: np.ndarray, z: int, t: int) -> np.ndarray:
    """Guideline-directed medical therapy: moderate fixed doses.

    Represents standard-of-care:
    - Moderate statin (u₁ = 0.5)
    - No anti-inflammatory (u₂ = 0.0, not yet standard of care)
    - Moderate antihypertensive (u₃ = 0.4)
    - Moderate exercise (u₄ = 0.4)
    """
    u = np.array([0.5, 0.0, 0.4, 0.4])
    return np.clip(u, 0.0, U_MAX_BASIN[z])


def static_aggressive(x: np.ndarray, z: int, t: int) -> np.ndarray:
    """Aggressive fixed therapy (post-ACS level).

    - High-intensity statin (u₁ = 0.9)
    - Low-dose colchicine (u₂ = 0.5)
    - Moderate-high antihypertensive (u₃ = 0.6)
    - Moderate exercise (u₄ = 0.4)
    """
    u = np.array([0.9, 0.5, 0.6, 0.4])
    return np.clip(u, 0.0, U_MAX_BASIN[z])


def _solve_lqr(A_k: np.ndarray, B_k: np.ndarray) -> np.ndarray:
    """Solve the infinite-horizon discrete LQR for basin k.

    Returns gain K such that u = -K x minimises
    sum x'Qx + λ u'Ru.
    """
    Q_lqr = Q_COST.copy()
    R_lqr = R_COST * LAMBDA_U
    try:
        P = solve_discrete_are(A_k, B_k, Q_lqr, R_lqr)
        K = np.linalg.solve(
            R_lqr + B_k.T @ P @ B_k,
            B_k.T @ P @ A_k,
        )
        return K
    except (np.linalg.LinAlgError, ValueError):
        # Fallback: zero gain
        return np.zeros((N_CONTROL, N_STATE))


# Pre-compute per-basin LQR gains
_K_BASIN = [_solve_lqr(A[k], B[k]) for k in range(N_BASINS)]

# Pooled LQR: average dynamics
_A_pool = np.mean(A, axis=0)
_B_pool = np.mean(B, axis=0)
_K_pool = _solve_lqr(_A_pool, _B_pool)


def lqr_basin_aware(x: np.ndarray, z: int, t: int) -> np.ndarray:
    """Per-basin LQR — the HDR Mode A approximation.

    Uses the correct basin-specific gain K_z.
    """
    u = -_K_BASIN[z] @ x
    u = np.clip(u, 0.0, U_MAX_BASIN[z])
    return u


def lqr_pooled(x: np.ndarray, z: int, t: int) -> np.ndarray:
    """Pooled LQR — single gain, basin-unaware.

    Baseline comparator: what you get without basin identification.
    """
    u = -_K_pool @ x
    u = np.clip(u, 0.0, U_MAX_BASIN[z])
    return u


# Registry for easy access
CONTROLLERS = {
    "open_loop": open_loop,
    "guideline": static_guideline,
    "aggressive": static_aggressive,
    "lqr_basin": lqr_basin_aware,
    "lqr_pooled": lqr_pooled,
}
