# Testing Guide

## Running Tests

```bash
# Install in dev mode (includes pytest)
pip install -e ".[dev]"

# Run all tests with verbose output
python -m pytest tests/ -v

# Run a specific test class
python -m pytest tests/test_athero.py::TestStability -v

# Run a single test
python -m pytest tests/test_athero.py::TestSimulation::test_lqr_reduces_cost_vs_open_loop -v

# Run with coverage (if pytest-cov installed)
python -m pytest tests/ -v --cov=hdr_athero --cov-report=term-missing
```

---

## Test Suite Architecture

The test suite (`tests/test_athero.py`) is organised into 6 test classes covering the full verification stack from parameter correctness through to end-to-end simulation integrity.

### Test Class 1: TestParameterShapes (6 tests)

Validates that all parameter matrices have correct dimensions:

| Test | Validates | Expected Shape |
|------|----------|----------------|
| `test_tau_shape` | TAU decay constants | (3, 4) |
| `test_A_shapes` | A_k dynamics matrices | 3 x (4, 4) |
| `test_B_shapes` | B_k intervention maps | 3 x (4, 4) |
| `test_J_shapes` | J^(k) coupling matrices | 3 x (4, 4) |
| `test_Q_shapes` | Q_k noise covariances | 3 x (4, 4) |
| `test_cost_matrices` | Q_COST, R_COST | (4,4), (4,4) |

**Why these matter**: Dimension mismatches would cause silent broadcasting errors in numpy, producing incorrect but finite results.

### Test Class 2: TestParameterProperties (7 tests)

Validates biological constraints that must hold for the model to be physically meaningful:

| Test | Constraint | Biological Reason |
|------|-----------|-------------------|
| `test_tau_positive` | TAU > 0 | Decay constants must be positive (finite recovery time) |
| `test_J_zero_diagonal` | diag(J) = 0 | No self-coupling on diagonal (handled by D matrix) |
| `test_J_nonnegative_offdiag` | J off-diag >= 0 | Disease coupling is self-reinforcing |
| `test_B_therapeutic_sign` | B <= 0 | Interventions drive state toward zero |
| `test_Q_spd` | Q is SPD | Noise covariance must be symmetric positive definite |
| `test_P_trans_rows_valid` | P_TRANS rows sum <= 1 | Valid probability distribution |
| `test_basin1_tau_larger` | TAU[1] >= TAU[0] | Vulnerable basin recovers more slowly |

### Test Class 3: TestStability (3 tests)

Validates spectral properties ensuring the system is stable (no unbounded growth):

| Test | Validates | Criterion |
|------|----------|-----------|
| `test_all_basins_stable` | All rho(A_k) < 1 | Global stability |
| `test_spectral_ordering` | rho(post_ACS) < rho(subclinical) < rho(vulnerable) | Basin severity ordering |
| `test_vulnerable_near_unit_root` | rho(Basin 1) in [0.88, 0.998] | ICI regime check |

### Test Class 4: TestControllability (2 tests)

Validates that the system is controllable in every basin:

| Test | Validates | Criterion |
|------|----------|-----------|
| `test_all_basins_controllable` | Controllability rank = 4 | Full state dimension reachable |
| `test_gramian_positive_definite` | Gramian min eigenvalue > 0 | All directions reachable |

### Test Class 5: TestSimulation (6 tests)

End-to-end integration tests:

| Test | What It Tests | Duration |
|------|-------------|----------|
| `test_open_loop_runs` | Open-loop produces finite states/costs | T=30 |
| `test_lqr_basin_runs` | LQR controller produces finite states | T=30 |
| `test_lqr_reduces_cost_vs_open_loop` | LQR mean cost < open-loop (10 episodes, Basin 0) | T=90 |
| `test_basin_transitions_occur` | Multiple basins visited | T=5000 |
| `test_controls_within_bounds` | All controllers: u in [0, u_max] | T=30 |
| `test_recovery_surrogate_amplification` | Basin 1 tau_tilde > 3x Basin 0 | N/A |

---

## What the Tests Prove

Collectively, the 24 tests establish:

1. **Dimensional correctness**: All matrix operations will produce correctly-shaped outputs
2. **Physical validity**: Parameters satisfy the biological constraints assumed by HDR theory
3. **Stability**: The system will not diverge in any basin (bounded trajectories)
4. **Controllability**: Every state-space direction is reachable from the control inputs
5. **Simulation integrity**: The complete simulation pipeline produces finite, bounded results
6. **Controller effectiveness**: LQR demonstrably outperforms open-loop in the stable regime
7. **Recovery surrogate**: The ICI effect (Basin 1 amplification) is quantitatively present

---

## Adding New Tests

When extending the model, add tests that verify:

1. **New parameters**: Shape, sign, and magnitude constraints
2. **New basins**: Stability (rho < 1), controllability (rank = n), and simulation viability
3. **New controllers**: Control bounds respected, finite costs, improvement over baseline
4. **New axes**: Update all shape tests and cross-coupling property tests

Example test for a new control channel:

```python
def test_new_control_bounds(self):
    """New controller respects u in [0, u_max]."""
    model = AtheroSLDS(seed=42)
    x0 = np.array([2.0, 2.0, 2.0, 2.0])
    episode = model.simulate(x0, 0, new_controller, T=30)
    assert np.all(episode.controls >= -1e-10)
    assert np.all(episode.controls <= 1.0 + 1e-10)
```
