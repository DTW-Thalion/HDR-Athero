"""
Spectral and controllability analysis for the atherosclerosis SLDS.
"""
from __future__ import annotations

import numpy as np

from hdr_athero.model.parameters import A, B, N_STATE, N_BASINS, AXIS_NAMES, BASIN_NAMES


def spectral_analysis() -> dict:
    """Compute spectral properties of all basin dynamics matrices.

    Returns dict with eigenvalues, spectral radii, and recovery amplification.
    """
    results = {}
    for k in range(N_BASINS):
        eigvals = np.linalg.eigvals(A[k])
        rho = float(np.max(np.abs(eigvals)))
        # Recovery amplification factor: 1/(1-ρ²)
        if rho < 1.0:
            amplification = 1.0 / (1.0 - rho**2)
        else:
            amplification = float("inf")
        results[BASIN_NAMES[k]] = {
            "eigenvalues": eigvals,
            "spectral_radius": rho,
            "recovery_amplification": amplification,
            "stable": rho < 1.0,
        }
    return results


def controllability_gramian(basin: int, T: int = 50) -> np.ndarray:
    """Compute the T-step controllability Gramian for basin k.

    W_c(T) = Σ_{t=0}^{T-1} A_k^t B_k B_k' (A_k')^t
    """
    Ak = A[basin]
    Bk = B[basin]
    W = np.zeros((N_STATE, N_STATE))
    Ak_power = np.eye(N_STATE)
    for _ in range(T):
        W += Ak_power @ Bk @ Bk.T @ Ak_power.T
        Ak_power = Ak_power @ Ak
    return W


def controllability_report() -> dict:
    """Generate controllability analysis for all basins.

    Returns dict with Gramian eigenvalues, rank, and per-axis reachability.
    """
    results = {}
    for k in range(N_BASINS):
        W = controllability_gramian(k, T=50)
        eigvals = np.linalg.eigvalsh(W)
        rank = int(np.sum(eigvals > 1e-10))
        # Per-axis: diagonal of Gramian indicates reachability energy
        diag = np.diag(W)

        results[BASIN_NAMES[k]] = {
            "gramian_eigenvalues": np.sort(eigvals)[::-1],
            "rank": rank,
            "full_rank": rank == N_STATE,
            "per_axis_energy": {
                AXIS_NAMES[i]: float(diag[i]) for i in range(N_STATE)
            },
            "min_eigval": float(eigvals.min()),
            "condition_number": float(eigvals.max() / max(eigvals.min(), 1e-15)),
        }
    return results


def print_spectral_report():
    """Print a formatted spectral analysis report."""
    results = spectral_analysis()
    print("\n" + "=" * 65)
    print("  SPECTRAL ANALYSIS — Atherosclerosis SLDS (4-axis)")
    print("=" * 65)
    for name, r in results.items():
        print(f"\n  Basin: {name}")
        print(f"    ρ(A_k) = {r['spectral_radius']:.4f}  "
              f"{'(STABLE)' if r['stable'] else '*** UNSTABLE ***'}")
        print(f"    Recovery amplification 1/(1-ρ²) = "
              f"{r['recovery_amplification']:.2f}×")
        eigvals = r["eigenvalues"]
        print(f"    Eigenvalues: {', '.join(f'{e:.4f}' for e in sorted(np.abs(eigvals), reverse=True))}")
    print()


def print_controllability_report():
    """Print a formatted controllability report."""
    results = controllability_report()
    print("\n" + "=" * 65)
    print("  CONTROLLABILITY ANALYSIS — 50-step Gramian")
    print("=" * 65)
    for name, r in results.items():
        print(f"\n  Basin: {name}")
        print(f"    Rank: {r['rank']}/{N_STATE}  "
              f"{'(CONTROLLABLE)' if r['full_rank'] else '*** NOT FULL RANK ***'}")
        print(f"    Gramian eigenvalues: "
              f"{', '.join(f'{e:.4f}' for e in r['gramian_eigenvalues'])}")
        print(f"    Condition number: {r['condition_number']:.1f}")
        print(f"    Per-axis reachability energy:")
        for axis, energy in r["per_axis_energy"].items():
            bar = "█" * min(int(energy * 10), 40)
            print(f"      {axis:4s}: {energy:8.4f}  {bar}")
    print()
