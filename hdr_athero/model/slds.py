"""
Switched Linear Dynamical System model for atherosclerosis.

Provides the core simulation loop:
    x_{t+1} = A_{z_t} x_t + B_{z_t} u_t + b_{z_t} + w_t
    z_{t+1} ~ HSMM(z_t, h_k(ℓ), P_trans)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from hdr_athero.model.parameters import (
    A, B, Q, P_TRANS, N_STATE, N_CONTROL, N_BASINS, DT, TAU,
    TARGET_RADIUS, Q_COST, W1, W2, LAMBDA_U, R_COST,
)


@dataclass
class StepResult:
    """Result of a single simulation step."""
    x: np.ndarray         # state after transition
    z: int                # basin index after transition
    u: np.ndarray         # control applied
    cost: float           # per-step cost
    recovery_surr: float  # τ̃ value
    rho: float            # spectral radius of current basin


@dataclass
class EpisodeResult:
    """Result of a full episode simulation."""
    states: np.ndarray       # (T+1, n)
    basins: np.ndarray       # (T+1,) int
    controls: np.ndarray     # (T, m)
    costs: np.ndarray        # (T,)
    recovery: np.ndarray     # (T,)
    rhos: np.ndarray         # (T,)
    total_cost: float
    mean_cost: float
    time_in_target: float    # fraction of steps with ||x|| ≤ target


class AtheroSLDS:
    """4-axis atherosclerosis SLDS model.

    Parameters
    ----------
    seed : int
        Random seed for reproducibility.
    """

    def __init__(self, seed: int = 42):
        self.rng = np.random.default_rng(seed)
        self.A = [a.copy() for a in A]
        self.B = [b.copy() for b in B]
        self.Q = [q.copy() for q in Q]

        # Pre-compute spectral radii
        self.rho = np.array([
            float(np.max(np.abs(np.linalg.eigvals(a)))) for a in self.A
        ])

    def spectral_radii(self) -> np.ndarray:
        """Return spectral radii for all basins."""
        return self.rho.copy()

    def recovery_surrogate(self, x: np.ndarray, z: int) -> float:
        """Compute τ̃(x, z) = ||x - Π(x, S*)||²_Q / (1 - ρ²).

        For simplicity, Π projects to the target box [-r, r]^n.
        """
        x_proj = np.clip(x, -TARGET_RADIUS, TARGET_RADIUS)
        deviation = x - x_proj
        quad = float(deviation @ Q_COST @ deviation)
        rho_k = self.rho[z]
        if rho_k >= 1.0:
            return float("inf")
        return quad / (1.0 - rho_k**2)

    def step_cost(self, x: np.ndarray, u: np.ndarray, z: int) -> float:
        """Compute per-step composite cost (eq. 7.5 of HDR v5.4)."""
        x_proj = np.clip(x, -TARGET_RADIUS, TARGET_RADIUS)
        deviation = x - x_proj

        normative = float(deviation @ Q_COST @ deviation)
        tau_surr = self.recovery_surrogate(x, z)
        effort = float(u @ R_COST @ u)

        return W1 * normative + W2 * tau_surr + LAMBDA_U * effort

    def _transition_basin(self, z: int) -> int:
        """Sample basin transition per HSMM transition matrix."""
        probs = P_TRANS[z].copy()
        p_stay = 1.0 - probs.sum()
        full_probs = np.zeros(N_BASINS)
        full_probs[:] = probs
        full_probs[z] = p_stay
        return int(self.rng.choice(N_BASINS, p=full_probs))

    def step(self, x: np.ndarray, z: int, u: np.ndarray) -> StepResult:
        """Advance one time step.

        x_{t+1} = A_z x_t + B_z u_t + w_t
        z_{t+1} ~ transition(z_t)
        """
        # State dynamics
        w = self.rng.multivariate_normal(np.zeros(N_STATE), self.Q[z])
        x_next = self.A[z] @ x + self.B[z] @ u + w

        # Basin transition
        z_next = self._transition_basin(z)

        # Cost
        cost = self.step_cost(x, u, z)
        tau = self.recovery_surrogate(x, z)

        return StepResult(
            x=x_next,
            z=z_next,
            u=u.copy(),
            cost=cost,
            recovery_surr=tau,
            rho=self.rho[z],
        )

    def simulate(
        self,
        x0: np.ndarray,
        z0: int,
        controller,
        T: int,
    ) -> EpisodeResult:
        """Run a full episode of T steps.

        Parameters
        ----------
        x0 : initial state (n,)
        z0 : initial basin index
        controller : callable(x, z, t) → u
        T : number of time steps
        """
        states = np.zeros((T + 1, N_STATE))
        basins = np.zeros(T + 1, dtype=int)
        controls = np.zeros((T, N_CONTROL))
        costs = np.zeros(T)
        recovery = np.zeros(T)
        rhos = np.zeros(T)

        states[0] = x0
        basins[0] = z0
        x, z = x0.copy(), z0

        for t in range(T):
            u = controller(x, z, t)
            result = self.step(x, z, u)

            controls[t] = result.u
            costs[t] = result.cost
            recovery[t] = result.recovery_surr
            rhos[t] = result.rho

            x = result.x
            z = result.z
            states[t + 1] = x
            basins[t + 1] = z

        # Time in target
        in_target = np.all(
            np.abs(states[:-1]) <= TARGET_RADIUS[None, :], axis=1
        )

        return EpisodeResult(
            states=states,
            basins=basins,
            controls=controls,
            costs=costs,
            recovery=recovery,
            rhos=rhos,
            total_cost=float(costs.sum()),
            mean_cost=float(costs.mean()),
            time_in_target=float(in_target.mean()),
        )
