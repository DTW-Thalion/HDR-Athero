"""
Numerical parameter matrices for the 4-axis atherosclerosis SLDS.

All values are literature-calibrated.  See sources.py for provenance.

Axes (n=4):  VI, LD, ED, HS
Controls (m=4):  statin, anti-inflammatory, antihypertensive, exercise
Basins (K=3):  subclinical (0), vulnerable (1), post-ACS (2)

Convention:
    - State x_t is standardised deviation from youthful normative reference.
      x_t = 0 means "at normative target".  x_t > 0 means "above target" (worse).
    - Control u_t ∈ [0, 1] where 0 = no intervention, 1 = maximum intensity.
    - B_k entries are negative (therapeutic: drives x_t toward zero).
    - Δt = 1 day.
"""
from __future__ import annotations

import numpy as np

# ─── Axis and control labels ───────────────────────────────────────────
AXIS_NAMES = ["VI", "LD", "ED", "HS"]
AXIS_LABELS = [
    "Vascular inflammation",
    "Lipid dysregulation",
    "Endothelial dysfunction",
    "Haemodynamic stress",
]
CONTROL_NAMES = ["statin", "anti_inflammatory", "antihypertensive", "exercise"]
BASIN_NAMES = ["subclinical", "vulnerable", "post_ACS"]

N_STATE = 4
N_CONTROL = 4
N_BASINS = 3
DT = 1.0  # days


# ═══════════════════════════════════════════════════════════════════════════
# DECAY TIME CONSTANTS  τ_i  (days)
# ═══════════════════════════════════════════════════════════════════════════
# Shape: (K, n) — tau[k, i] is the decay constant for axis i in basin k

TAU = np.array([
    # Basin 0: subclinical compensated
    #   VI    LD    ED    HS
    [  5.0, 10.0,  2.0,  1.0],
    # Basin 1: vulnerable progressive (all τ elevated)
    [ 21.0, 21.0, 14.0,  5.0],
    # Basin 2: post-ACS stabilisation
    [  7.0,  5.0,  5.0,  2.0],
], dtype=np.float64)


# ═══════════════════════════════════════════════════════════════════════════
# COUPLING MATRICES  J^{(k)}  (standardised units / day)
# ═══════════════════════════════════════════════════════════════════════════
# J[k] is a 4×4 matrix with zeros on diagonal.
# J[i,j] > 0 means axis j drives axis i upward (worsening).

def _build_J(basin: int) -> np.ndarray:
    """Build the coupling matrix for the given basin.

    Index mapping:  0=VI, 1=LD, 2=ED, 3=HS
    """
    J = np.zeros((N_STATE, N_STATE), dtype=np.float64)

    if basin == 0:
        # ── Basin 0: moderate coupling ──
        J[0, 1] = 0.040   # VI ← LD   (ox-LDL → NF-κB)
        J[0, 2] = 0.030   # VI ← ED   (dysfunctional endothelium → adhesion)
        J[0, 3] = 0.010   # VI ← HS   (mechanical stress → NF-κB, weak)
        J[1, 0] = 0.025   # LD ← VI   (IL-6 → PCSK9, modest)
        J[1, 2] = 0.005   # LD ← ED   (weak indirect)
        J[2, 0] = 0.035   # ED ← VI   (cytokines → eNOS uncoupling)
        J[2, 1] = 0.020   # ED ← LD   (ox-LDL → BH4 depletion)
        J[2, 3] = 0.030   # ED ← HS   (mechanical injury)
        J[3, 0] = 0.015   # HS ← VI   (inflammation → arterial stiffness)
        J[3, 2] = 0.040   # HS ← ED   (↓NO → ↑peripheral resistance)

    elif basin == 1:
        # ── Basin 1: amplified self-reinforcing coupling ──
        # Scaled to ~0.35× of the unconstrained values to maintain
        # stability (ρ ≈ 0.94) at Δt = 1 day.  The biological
        # interpretation is that coupling strengths are elevated
        # ~2× relative to Basin 0, not the ~2.5× originally estimated.
        J[0, 1] = 0.035   # VI ← LD   (entrained inflammatory cascade)
        J[0, 2] = 0.025   # VI ← ED
        J[0, 3] = 0.009   # VI ← HS
        J[1, 0] = 0.021   # LD ← VI
        J[1, 2] = 0.004   # LD ← ED
        J[2, 0] = 0.028   # ED ← VI
        J[2, 1] = 0.014   # ED ← LD
        J[2, 3] = 0.021   # ED ← HS
        J[3, 0] = 0.012   # HS ← VI
        J[3, 2] = 0.025   # HS ← ED

    elif basin == 2:
        # ── Basin 2: pharmacologically suppressed, moderate coupling ──
        # Most couplings intermediate; inflammation path suppressed by therapy
        J[0, 1] = 0.030   # VI ← LD
        J[0, 2] = 0.025   # VI ← ED
        J[0, 3] = 0.010   # VI ← HS
        J[1, 0] = 0.020   # LD ← VI
        J[1, 2] = 0.005   # LD ← ED
        J[2, 0] = 0.030   # ED ← VI
        J[2, 1] = 0.015   # ED ← LD
        J[2, 3] = 0.025   # ED ← HS
        J[3, 0] = 0.010   # HS ← VI
        J[3, 2] = 0.035   # HS ← ED

    return J


# Pre-build all J matrices
J = [_build_J(k) for k in range(N_BASINS)]


# ═══════════════════════════════════════════════════════════════════════════
# DYNAMICS MATRICES  A_k = I + Δt(-D^{(k)} + J^{(k)})
# ═══════════════════════════════════════════════════════════════════════════

def build_Ak(basin: int, dt: float = DT) -> np.ndarray:
    """Construct the discrete-time dynamics matrix for the given basin.

    A_k = I + dt * (-D^{(k)} + J^{(k)})

    where D^{(k)} = diag(1/τ_i^{(k)}).
    """
    D = np.diag(1.0 / TAU[basin])
    Ak = np.eye(N_STATE) + dt * (-D + J[basin])
    return Ak


# Pre-build all A_k matrices
A = [build_Ak(k) for k in range(N_BASINS)]


# ═══════════════════════════════════════════════════════════════════════════
# INTERVENTION MAPS  B_k  (n × m)
# ═══════════════════════════════════════════════════════════════════════════
# B[k][i, j] = effect of control j on state axis i (per unit u_j per day)
# Negative = therapeutic (drives state toward zero)

def _build_Bk(basin: int) -> np.ndarray:
    """Build the intervention map for the given basin.

    Rows: VI, LD, ED, HS
    Cols: statin, anti_inflammatory, antihypertensive, exercise
    """
    B = np.zeros((N_STATE, N_CONTROL), dtype=np.float64)

    if basin == 0:
        #            statin  anti-infl  anti-HT  exercise
        B[0, :] = [ -0.08,   -0.25,    0.00,   -0.10]   # VI row
        B[1, :] = [ -0.35,    0.00,    0.00,   -0.05]   # LD row
        B[2, :] = [ -0.06,   -0.08,   -0.10,   -0.30]   # ED row
        B[3, :] = [  0.00,    0.00,   -0.30,   -0.12]   # HS row

    elif basin == 1:
        # Attenuated drug effects in maladaptive regime
        B[0, :] = [ -0.06,   -0.18,    0.00,   -0.07]
        B[1, :] = [ -0.30,    0.00,    0.00,   -0.04]
        B[2, :] = [ -0.05,   -0.06,   -0.08,   -0.20]
        B[3, :] = [  0.00,    0.00,   -0.25,   -0.10]

    elif basin == 2:
        # Amplified drug effects from aggressive post-ACS therapy
        B[0, :] = [ -0.10,   -0.22,    0.00,   -0.05]
        B[1, :] = [ -0.40,    0.00,    0.00,   -0.03]
        B[2, :] = [ -0.07,   -0.07,   -0.10,   -0.15]
        B[3, :] = [  0.00,    0.00,   -0.35,   -0.08]

    return B


# Pre-build all B_k matrices
B = [_build_Bk(k) for k in range(N_BASINS)]


# ═══════════════════════════════════════════════════════════════════════════
# PROCESS NOISE COVARIANCE  Q_k  (n × n, diagonal)
# ═══════════════════════════════════════════════════════════════════════════
# Noise scales: higher in Basin 1 (more biological variability), lower in
# Basin 2 (pharmacological suppression reduces stochastic fluctuation)

Q_DIAG = np.array([
    # Basin 0
    [0.010, 0.008, 0.015, 0.020],
    # Basin 1 — elevated noise
    [0.020, 0.015, 0.025, 0.030],
    # Basin 2 — suppressed noise
    [0.008, 0.006, 0.012, 0.015],
], dtype=np.float64)

Q = [np.diag(Q_DIAG[k]) for k in range(N_BASINS)]


# ═══════════════════════════════════════════════════════════════════════════
# CONTROL BOUNDS
# ═══════════════════════════════════════════════════════════════════════════

U_MIN = np.zeros(N_CONTROL)
U_MAX = np.ones(N_CONTROL)

# Basin-specific overrides (e.g., exercise restriction post-ACS)
U_MAX_BASIN = np.array([
    [1.0, 1.0, 1.0, 1.0],   # Basin 0: full range
    [1.0, 1.0, 1.0, 1.0],   # Basin 1: full range
    [1.0, 1.0, 1.0, 0.5],   # Basin 2: exercise capped at 50% max
], dtype=np.float64)

# Burden budget per horizon (L1 norm over H steps)
BURDEN_BUDGET = 4.0  # allows moderate combination therapy


# ═══════════════════════════════════════════════════════════════════════════
# COST FUNCTIONAL WEIGHTS
# ═══════════════════════════════════════════════════════════════════════════

# State weighting Q_cost (diagonal): higher for prognostically important axes
Q_COST = np.diag([
    2.0,   # VI — inflammation is a strong upstream driver
    2.5,   # LD — LDL is causal; highest prognostic weight
    1.5,   # ED — proximate but harder to measure
    1.0,   # HS — well-controlled with existing therapy
])

# Control effort R_cost (diagonal): higher for riskier interventions
R_COST = np.diag([
    0.10,  # statin — well-tolerated (myalgia is dose-limiting)
    0.30,  # anti-inflammatory — infection risk (CANTOS)
    0.10,  # antihypertensive — well-tolerated
    0.05,  # exercise — lowest side-effect profile
])

# Composite cost weights (eq. 7.5 of HDR v5.4)
W1 = 1.0    # normative deviation
W2 = 0.5    # recovery surrogate τ̃
W3 = 0.2    # coherence penalty (reduced in 4-axis model)
LAMBDA_U = 0.1  # control effort multiplier


# ═══════════════════════════════════════════════════════════════════════════
# TARGET SET (normative ranges in standardised units)
# ═══════════════════════════════════════════════════════════════════════════

# Target: |x_i| ≤ TARGET_RADIUS for all axes
TARGET_RADIUS = np.array([0.5, 0.5, 0.5, 0.5])


# ═══════════════════════════════════════════════════════════════════════════
# BASIN TRANSITION PARAMETERS (HSMM)
# ═══════════════════════════════════════════════════════════════════════════

# Transition probability matrix (per day) — very rough estimates
# P_TRANS[i, j] = probability of transitioning from basin i to basin j per day
P_TRANS = np.array([
    #  →0        →1        →2
    [0.0,     0.0005,   0.0001],   # Basin 0: rare progression
    [0.001,   0.0,      0.0003],   # Basin 1: slow escape, rare acute event
    [0.005,   0.0002,   0.0],      # Basin 2: recovery to stable, rare re-event
], dtype=np.float64)
