# HDR-Athero Model Specification

## 1. System Definition

The HDR-Athero model is a **switched linear dynamical system (SLDS)** for atherosclerotic cardiovascular disease with:

- **n = 4** state dimensions (axes)
- **m = 4** control channels (interventions)
- **K = 3** basins (disease regimes)
- **dt = 1 day** (discrete time step)

### 1.1 State Axes

| Index | Symbol | Full Name | Primary Biomarker | Population SD |
|-------|--------|-----------|-------------------|---------------|
| 0 | x_VI | Vascular inflammation | hsCRP, IL-6 | ~2 mg/L |
| 1 | x_LD | Lipid dysregulation | LDL-C, ApoB | ~1 mmol/L |
| 2 | x_ED | Endothelial dysfunction | FMD (%) | ~3-4% |
| 3 | x_HS | Haemodynamic stress | Systolic BP | ~15 mmHg |

**Convention**: `x_t = 0` means "at youthful normative target". `x_t > 0` means "above target" (worse). All states are standardised deviations.

### 1.2 Control Channels

| Index | Symbol | Intervention | Examples |
|-------|--------|-------------|----------|
| 0 | u_1 | Lipid-lowering | Atorvastatin 80mg, rosuvastatin 40mg, PCSK9i |
| 1 | u_2 | Anti-inflammatory | Colchicine 0.5mg, canakinumab 150mg |
| 2 | u_3 | Antihypertensive | ACEi/ARB + CCB combination |
| 3 | u_4 | Exercise | Structured aerobic/resistance programme |

**Convention**: `u_t in [0, 1]` where 0 = no intervention, 1 = maximum tolerated intensity.

### 1.3 Basins

| Index | Name | Spectral Radius | Clinical Regime |
|-------|------|-----------------|-----------------|
| 0 | subclinical | rho = 0.9101 | Compensated atherosclerosis, stable |
| 1 | vulnerable | rho = 0.9909 | Progressive disease, self-reinforcing |
| 2 | post_ACS | rho = 0.8780 | Post-event stabilisation under therapy |

---

## 2. Dynamics Equations

### 2.1 State Transition

```
x_{t+1} = A_z * x_t + B_z * u_t + w_t
```

where:
- `z` is the current basin index
- `w_t ~ N(0, Q_z)` is basin-dependent process noise
- `A_z` is the dynamics matrix for basin z
- `B_z` is the intervention map for basin z

### 2.2 Dynamics Matrix Construction

```
A_k = I + dt * (-D^(k) + J^(k))
```

where:
- `D^(k) = diag(1/tau_i^(k))` is the natural decay matrix
- `J^(k)` is the inter-axis coupling matrix
- `I` is the 4x4 identity matrix

### 2.3 Decay Time Constants (tau, days)

| Axis | Basin 0 (subclinical) | Basin 1 (vulnerable) | Basin 2 (post_ACS) |
|------|----------------------|---------------------|-------------------|
| VI | 5.0 | 21.0 | 7.0 |
| LD | 10.0 | 21.0 | 5.0 |
| ED | 2.0 | 14.0 | 5.0 |
| HS | 1.0 | 5.0 | 2.0 |

Basin 1 time constants are uniformly elevated, reflecting slowed recovery in the maladaptive regime.

### 2.4 Coupling Matrix Structure

The coupling matrix `J^(k)` has:
- **Zero diagonal**: no self-coupling
- **Non-negative off-diagonal**: disease axes reinforce each other

Key coupling pathways:
- **VI <-- LD**: ox-LDL activates NF-kB and NLRP3 inflammasome
- **LD <-- VI**: IL-6/STAT3 upregulates hepatic PCSK9
- **ED <-- VI**: Inflammatory cytokines reduce NO bioavailability
- **ED <-- HS**: Elevated BP causes mechanical endothelial injury
- **HS <-- ED**: Endothelial dysfunction increases peripheral resistance

### 2.5 Intervention Map (B_k)

The B_k matrix encodes the per-day therapeutic effect of each control on each state axis. All entries are non-positive (therapeutic sign convention).

**Basin 0 (subclinical)**:

|  | statin | anti-infl | anti-HT | exercise |
|--|--------|-----------|---------|----------|
| VI | -0.08 | -0.25 | 0.00 | -0.10 |
| LD | -0.35 | 0.00 | 0.00 | -0.05 |
| ED | -0.06 | -0.08 | -0.10 | -0.30 |
| HS | 0.00 | 0.00 | -0.30 | -0.12 |

**Basin 1 (vulnerable)** -- attenuated efficacy:

|  | statin | anti-infl | anti-HT | exercise |
|--|--------|-----------|---------|----------|
| VI | -0.06 | -0.18 | 0.00 | -0.07 |
| LD | -0.30 | 0.00 | 0.00 | -0.04 |
| ED | -0.05 | -0.06 | -0.08 | -0.20 |
| HS | 0.00 | 0.00 | -0.25 | -0.10 |

**Basin 2 (post_ACS)** -- amplified drug response:

|  | statin | anti-infl | anti-HT | exercise |
|--|--------|-----------|---------|----------|
| VI | -0.10 | -0.22 | 0.00 | -0.05 |
| LD | -0.40 | 0.00 | 0.00 | -0.03 |
| ED | -0.07 | -0.07 | -0.10 | -0.15 |
| HS | 0.00 | 0.00 | -0.35 | -0.08 |

---

## 3. Basin Switching (HSMM)

Basin transitions follow a Hidden Semi-Markov Model with per-day transition probabilities:

```
P_TRANS[i, j] = probability of transitioning from basin i to basin j per day
```

|  | to Basin 0 | to Basin 1 | to Basin 2 |
|--|-----------|-----------|-----------|
| from Basin 0 | -- | 0.0005 | 0.0001 |
| from Basin 1 | 0.001 | -- | 0.0003 |
| from Basin 2 | 0.005 | 0.0002 | -- |

Expected dwell times:
- Basin 0: ~1667 days (4.6 years) before transition
- Basin 1: ~769 days (2.1 years) before transition
- Basin 2: ~192 days (6.3 months) before transition

---

## 4. Cost Functional

The per-step composite cost (HDR v5.4, eq. 7.5):

```
J(x, u, z) = W1 * ||x - Pi(x, S*)||^2_Q + W2 * tau_tilde(x, z) + lambda_u * u' R u
```

where:
- `W1 = 1.0` (normative deviation weight)
- `W2 = 0.5` (recovery surrogate weight)
- `lambda_u = 0.1` (control effort multiplier)
- `Q_COST = diag([2.0, 2.5, 1.5, 1.0])` (state weights)
- `R_COST = diag([0.10, 0.30, 0.10, 0.05])` (control effort weights)
- `TARGET_RADIUS = [0.5, 0.5, 0.5, 0.5]` (normative target bounds)

### 4.1 Recovery Surrogate

```
tau_tilde(x, z) = ||x - Pi(x, S*)||^2_Q / (1 - rho(A_z)^2)
```

The recovery surrogate amplifies state deviations by the basin's spectral proximity to the unit circle. In Basin 1 (rho = 0.9909), this amplification factor is 55.44x compared to 5.82x in Basin 0.

---

## 5. Process Noise

Diagonal covariance matrices reflecting biological variability:

| Axis | Basin 0 | Basin 1 | Basin 2 |
|------|---------|---------|---------|
| VI | 0.010 | 0.020 | 0.008 |
| LD | 0.008 | 0.015 | 0.006 |
| ED | 0.015 | 0.025 | 0.012 |
| HS | 0.020 | 0.030 | 0.015 |

Basin 1 has elevated noise (more biological variability in the vulnerable regime). Basin 2 has suppressed noise (pharmacological stabilisation).

---

## 6. Control Bounds

- `u_min = [0, 0, 0, 0]` (no negative interventions)
- `u_max = [1, 1, 1, 1]` for Basins 0 and 1
- `u_max = [1, 1, 1, 0.5]` for Basin 2 (exercise restricted post-ACS)
- `BURDEN_BUDGET = 4.0` (L1 norm constraint over the planning horizon)
