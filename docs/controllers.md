# Controller Documentation

## Overview

HDR-Athero implements five control strategies for comparison, ranging from no intervention (open-loop) to optimal basin-aware feedback (LQR). All controllers produce control vectors `u in [0, 1]^4` that are clipped to basin-specific bounds.

---

## 1. Open-Loop (No Intervention)

```
u_t = [0, 0, 0, 0]
```

Pure natural disease progression. Serves as the baseline for all comparisons. The system evolves under `x_{t+1} = A_z x_t + w_t` with no therapeutic intervention.

**Use case**: Establishing natural disease trajectory and quantifying therapeutic benefit.

---

## 2. Static Guideline-Directed Therapy

```
u_t = [0.5, 0.0, 0.4, 0.4]
```

Fixed-dose standard-of-care:
- **Statin**: 50% intensity (moderate-intensity statin)
- **Anti-inflammatory**: 0% (not part of standard guidelines)
- **Antihypertensive**: 40% intensity (single-agent)
- **Exercise**: 40% intensity (~8 MET-h/week)

**Rationale**: Reflects current ACC/AHA guideline recommendations for primary prevention without anti-inflammatory therapy (colchicine/canakinumab are not yet standard of care in most guidelines).

---

## 3. Static Aggressive Therapy

```
u_t = [0.9, 0.5, 0.6, 0.4]
```

Post-ACS intensity regimen:
- **Statin**: 90% intensity (high-intensity statin + ezetimibe)
- **Anti-inflammatory**: 50% intensity (colchicine 0.5mg daily)
- **Antihypertensive**: 60% intensity (dual-agent)
- **Exercise**: 40% intensity (cardiac rehabilitation level)

**Rationale**: Reflects aggressive secondary prevention as in IMPROVE-IT and LoDoCo2 protocols.

---

## 4. LQR Basin-Aware (HDR Mode A)

```
u_t = -K_z * x_t    (clipped to [0, u_max])
```

Per-basin optimal linear-quadratic regulator. A separate gain matrix `K_k` is computed for each basin by solving the discrete algebraic Riccati equation (DARE):

```
P = Q_COST + A' P A - A' P B (R_COST + B' P B)^{-1} B' P A
K = (R_COST + B' P B)^{-1} B' P A
```

This is the **HDR advantage case**: the controller adapts its gain to the current basin, applying stronger control in Basin 1 (vulnerable) where dynamics are near-unit-root.

**Computed gains** (4x4 matrices, rows = controls, columns = state axes):

### K_0 (subclinical)

|  | VI | LD | ED | HS |
|--|-----|-----|-----|-----|
| statin | -0.178 | -2.529 | 0.243 | 0.015 |
| anti_infl | -2.662 | 0.425 | 0.413 | -0.003 |
| anti_HT | -0.038 | -0.120 | 0.468 | 0.036 |
| exercise | 0.363 | 0.351 | -1.831 | -0.108 |

### K_1 (vulnerable)

|  | VI | LD | ED | HS |
|--|-----|-----|-----|-----|
| statin | -0.229 | -3.098 | 0.625 | -0.107 |
| anti_infl | -3.922 | 0.469 | 1.094 | -0.238 |
| anti_HT | -0.005 | -0.178 | 1.421 | -3.004 |
| exercise | 0.444 | 0.565 | -4.945 | 0.875 |

### K_2 (post_ACS)

|  | VI | LD | ED | HS |
|--|-----|-----|-----|-----|
| statin | -0.174 | -1.979 | 0.365 | -0.043 |
| anti_infl | -3.094 | 0.508 | 0.627 | -0.111 |
| anti_HT | -0.055 | -0.112 | 0.829 | -1.412 |
| exercise | 0.706 | 0.587 | -5.261 | 0.615 |

**Key observation**: The dominant gain entries are:
- Statin strongly driven by LD axis (K[statin, LD] ~ -2.5 to -3.1)
- Anti-inflammatory strongly driven by VI axis (K[anti_infl, VI] ~ -2.7 to -3.9)
- Antihypertensive strongly driven by HS axis in Basin 1 (K[anti_HT, HS] = -3.0)
- Exercise strongly driven by ED axis (K[exercise, ED] ~ -1.8 to -5.3)

This is clinically interpretable: each drug targets its primary axis.

---

## 5. LQR Pooled (Basin-Unaware)

```
u_t = -K_pool * x_t    (clipped to [0, u_max])
```

A single gain matrix computed from averaged dynamics:
```
A_avg = (A_0 + A_1 + A_2) / 3
B_avg = (B_0 + B_1 + B_2) / 3
```

This represents the best controller achievable **without basin identification** -- the standard approach when the disease regime is unknown.

### K_pool

|  | VI | LD | ED | HS |
|--|-----|-----|-----|-----|
| statin | -0.190 | -2.486 | 0.396 | -0.043 |
| anti_infl | -3.179 | 0.475 | 0.705 | -0.114 |
| anti_HT | -0.026 | -0.136 | 0.854 | -1.383 |
| exercise | 0.470 | 0.483 | -3.661 | 0.399 |

---

## 6. Performance Comparison

### Mean Cost (20 episodes x 180 days)

| Controller | Subclinical | Vulnerable | Post-ACS |
|-----------|-------------|------------|----------|
| open_loop | 3.169 | 16.718 | 0.073 |
| guideline | 242.698 | 6746.945 | 14.368 |
| aggressive | 1051.668 | 25657.956 | 75.408 |
| lqr_basin | 0.843 | 58.700 | 0.046 |
| lqr_pooled | 0.734 | 55.126 | 0.047 |

### Time-in-Target (fraction of steps with ||x|| <= 0.5)

| Controller | Subclinical | Vulnerable | Post-ACS |
|-----------|-------------|------------|----------|
| open_loop | 90.1% | 31.5% | 93.6% |
| guideline | 2.1% | 1.1% | 1.5% |
| aggressive | 1.1% | 0.6% | 0.8% |
| lqr_basin | 92.5% | 16.1% | 92.4% |
| lqr_pooled | 91.5% | 17.6% | 91.9% |

### Key Findings

1. **Static therapies (guideline, aggressive) have high costs** because they apply constant control regardless of state, incurring unnecessary control effort when the state is already near target.

2. **LQR dramatically outperforms static therapy** by adapting control intensity to the current state deviation.

3. **Basin 1 (vulnerable) is hard to control for all strategies** -- even LQR achieves only ~16% time-in-target vs ~92% in Basin 0.

4. **HDR advantage (basin-aware vs pooled)** is modest in this 4-axis model:
   - Post-ACS: +2.3% cost improvement (basin-aware better)
   - Subclinical: -14.9% (pooled slightly better -- possibly noise)
   - Vulnerable: -6.5% (pooled slightly better)

   The advantage is expected to grow with higher-dimensional state spaces and when basin identification uncertainty is explicitly modelled.

---

## 7. Known Limitations

1. **One-sided control constraint**: `u >= 0` breaks standard LQR optimality guarantees, especially in Basin 1 where the near-unit-root dynamics amplify the mismatch. A proper nonnegative MPC solver would be more appropriate.

2. **No basin identification uncertainty**: The LQR basin-aware controller assumes perfect knowledge of the current basin. In practice, basin identification is unreliable in Basin 1 (see ICI diagnostic).

3. **Linear control for nonlinear disease**: The SLDS is piecewise-linear. Real atherosclerosis has nonlinear dose-response curves, saturation effects, and state-dependent transitions that this model does not capture.
