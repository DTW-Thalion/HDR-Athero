"""
Main simulation entry point.

Run with:  python -m hdr_athero.simulation.run

Executes:
  1. Spectral analysis of all A_k matrices
  2. Controllability analysis (Gramian)
  3. Comparative simulation: 5 controllers × 3 starting basins
  4. Recovery surrogate demonstration
  5. Summary report
"""
from __future__ import annotations

import sys
import time

import numpy as np

from hdr_athero.model.parameters import (
    N_STATE, N_BASINS, AXIS_NAMES, BASIN_NAMES, CONTROL_NAMES,
)
from hdr_athero.model.slds import AtheroSLDS
from hdr_athero.simulation.controller import CONTROLLERS
from hdr_athero.analysis.spectral import (
    print_spectral_report,
    print_controllability_report,
    spectral_analysis,
    controllability_report,
)


def _run_comparative_simulation(
    n_episodes: int = 20,
    T: int = 180,
    seed: int = 42,
) -> dict:
    """Run all controllers across all starting basins.

    Parameters
    ----------
    n_episodes : number of episodes per (controller, basin) pair
    T : episode length in days (default 180 = ~6 months)
    seed : base random seed

    Returns
    -------
    results : nested dict[controller][basin] → {mean_cost, std_cost, tib, ...}
    """
    results = {}

    for ctrl_name, ctrl_fn in CONTROLLERS.items():
        results[ctrl_name] = {}
        for z0 in range(N_BASINS):
            costs = []
            tibs = []
            for ep in range(n_episodes):
                model = AtheroSLDS(seed=seed + ep * 1000 + z0 * 100)
                # Initial state: moderately deviated
                x0 = model.rng.normal(scale=0.8, size=N_STATE)
                x0 = np.abs(x0)  # all positive (above target)

                episode = model.simulate(x0, z0, ctrl_fn, T)
                costs.append(episode.mean_cost)
                tibs.append(episode.time_in_target)

            results[ctrl_name][BASIN_NAMES[z0]] = {
                "mean_cost": float(np.mean(costs)),
                "std_cost": float(np.std(costs)),
                "mean_tib": float(np.mean(tibs)),
                "std_tib": float(np.std(tibs)),
            }

    return results


def _print_comparison_table(results: dict):
    """Print formatted comparison table."""
    print("\n" + "=" * 90)
    print("  COMPARATIVE SIMULATION — 20 episodes × 180 days per condition")
    print("=" * 90)

    for basin_name in BASIN_NAMES:
        print(f"\n  Starting Basin: {basin_name}")
        print(f"  {'Controller':<18s}  {'Mean Cost':>10s}  {'± Std':>8s}  "
              f"{'Time-in-Target':>14s}  {'± Std':>8s}")
        print("  " + "-" * 62)

        for ctrl_name in CONTROLLERS:
            r = results[ctrl_name][basin_name]
            print(f"  {ctrl_name:<18s}  {r['mean_cost']:10.3f}  "
                  f"{r['std_cost']:8.3f}  "
                  f"{r['mean_tib']:13.1%}  "
                  f"{r['std_tib']:8.3f}")


def _recovery_surrogate_demo():
    """Demonstrate the recovery surrogate across basins."""
    print("\n" + "=" * 65)
    print("  RECOVERY SURROGATE DEMONSTRATION")
    print("  Same state deviation, different basins")
    print("=" * 65)

    model = AtheroSLDS(seed=0)
    # Fixed moderate deviation
    x_test = np.array([1.0, 1.0, 1.0, 1.0])

    print(f"\n  State: x = [{', '.join(f'{v:.1f}' for v in x_test)}]")
    print(f"  (= 1 SD above target on all four axes)\n")

    for k in range(N_BASINS):
        tau = model.recovery_surrogate(x_test, k)
        rho = model.rho[k]
        amp = 1.0 / (1.0 - rho**2) if rho < 1.0 else float("inf")
        print(f"  Basin {k} ({BASIN_NAMES[k]}):")
        print(f"    ρ(A_k) = {rho:.4f}")
        print(f"    Amplification 1/(1-ρ²) = {amp:.2f}×")
        print(f"    τ̃(x, z={k}) = {tau:.3f}")
        print()


def _lqr_gain_analysis():
    """Show the per-basin LQR gain structure."""
    from hdr_athero.simulation.controller import _K_BASIN, _K_pool

    print("\n" + "=" * 65)
    print("  LQR GAIN ANALYSIS")
    print("=" * 65)

    for k in range(N_BASINS):
        K = _K_BASIN[k]
        print(f"\n  Basin {k} ({BASIN_NAMES[k]}) — K_{k}:")
        print(f"  {'':4s}", end="")
        for j, ax in enumerate(AXIS_NAMES):
            print(f"  {ax:>8s}", end="")
        print()
        for i, ctrl in enumerate(CONTROL_NAMES):
            print(f"  {ctrl:>14s}", end="")
            for j in range(N_STATE):
                print(f"  {K[i, j]:8.4f}", end="")
            print()

    print(f"\n  Pooled (basin-unaware) — K_pool:")
    print(f"  {'':4s}", end="")
    for j, ax in enumerate(AXIS_NAMES):
        print(f"  {ax:>8s}", end="")
    print()
    for i, ctrl in enumerate(CONTROL_NAMES):
        print(f"  {ctrl:>14s}", end="")
        for j in range(N_STATE):
            print(f"  {_K_pool[i, j]:8.4f}", end="")
        print()


def _ici_diagnostic():
    """Demonstrate the ICI insight: Basin 1 is hardest to identify."""
    print("\n" + "=" * 65)
    print("  ICI DIAGNOSTIC — Effective Sample Count")
    print("=" * 65)

    spec = spectral_analysis()
    print(f"\n  The autocorrelation time within each basin is ~1/(1-ρ):")
    for name, r in spec.items():
        rho = r["spectral_radius"]
        autocorr_time = 1.0 / (1.0 - rho) if rho < 1.0 else float("inf")
        # Assume 180-day dwell
        dwell = 180
        T_eff = dwell / autocorr_time
        print(f"    {name}: ρ={rho:.4f}  →  autocorr time = "
              f"{autocorr_time:.1f} days  →  T_eff(180d dwell) = "
              f"{T_eff:.1f} effective observations")

    print(f"\n  For n=4 state dimensions, minimum T_eff for reliable")
    print(f"  identification is ~40-80 (rule of thumb: 10-20× state dim).")
    print(f"  Basin 'vulnerable' is near or below this threshold.")
    print(f"  This is the ICI's fundamental insight: the hardest basin")
    print(f"  to control is also the hardest to identify.\n")


def main():
    """Run the full simulation suite."""
    t0 = time.time()

    print("\n" + "╔" + "═" * 63 + "╗")
    print("║  HDR-Athero: Literature-Calibrated Atherosclerosis SLDS" + " " * 7 + "║")
    print("║  4-axis model · 3 basins · 4 control channels" + " " * 16 + "║")
    print("╚" + "═" * 63 + "╝")

    # 1. Spectral analysis
    print_spectral_report()

    # 2. Controllability
    print_controllability_report()

    # 3. Recovery surrogate demo
    _recovery_surrogate_demo()

    # 4. LQR gains
    _lqr_gain_analysis()

    # 5. ICI diagnostic
    _ici_diagnostic()

    # 6. Comparative simulation
    print("\n  Running comparative simulation (20 episodes × 180 days × "
          "5 controllers × 3 basins)...")
    results = _run_comparative_simulation(n_episodes=20, T=180, seed=42)
    _print_comparison_table(results)

    # 7. HDR advantage analysis
    print("\n" + "=" * 65)
    print("  HDR ADVANTAGE ANALYSIS (LQR Basin-Aware vs Pooled)")
    print("=" * 65)
    for basin_name in BASIN_NAMES:
        hdr = results["lqr_basin"][basin_name]["mean_cost"]
        pooled = results["lqr_pooled"][basin_name]["mean_cost"]
        if pooled > 0:
            gain = (pooled - hdr) / pooled * 100
        else:
            gain = 0.0
        print(f"  {basin_name:<14s}: HDR={hdr:.3f}  Pooled={pooled:.3f}  "
              f"Gain={gain:+.1f}%")

    elapsed = time.time() - t0
    print(f"\n  Completed in {elapsed:.1f}s")
    print(f"\n  Disclaimer: All results are computational only.")
    print(f"  No clinical validity is asserted.\n")


if __name__ == "__main__":
    main()
