# HDR-Athero

**Literature-calibrated 4-axis atherosclerosis model for the Homeodynamic Remediation Framework (HDR v5.4)**

## Overview

This repository implements a reduced-dimension switched linear dynamical system (SLDS) for atherosclerotic cardiovascular disease, parameterised from published clinical trial data and meta-analyses.  It is a computational companion to the HDR Framework v5.4, instantiating the framework's core objects — basin dynamics matrices $A_k$, coupling matrices $J^{(k)}$, and intervention input matrices $B_k$ — for a 4-axis model spanning:

| Axis | Symbol | Primary Biomarker |
|------|--------|-------------------|
| Vascular inflammation | Δx_VI | hsCRP, IL-6 |
| Lipid dysregulation | Δx_LD | LDL-C, ApoB |
| Endothelial dysfunction | Δx_ED | FMD (%) |
| Haemodynamic stress | Δx_HS | Systolic BP |

Three basins model clinically distinct disease regimes:

- **Basin 0**: Subclinical compensated atherosclerosis (ρ ≈ 0.75)
- **Basin 1**: Vulnerable progressive disease (ρ ≈ 0.94) — the maladaptive basin
- **Basin 2**: Post-ACS stabilisation (ρ ≈ 0.62)

Four control channels map to the standard therapeutic arsenal:

| Channel | Intervention |
|---------|-------------|
| u₁ | Lipid-lowering (statin intensity / PCSK9i) |
| u₂ | Anti-inflammatory (colchicine / canakinumab) |
| u₃ | Antihypertensive (ACEi/ARB + CCB) |
| u₄ | Exercise (MET-hours/week, standardised) |

## What This Is

A **literature-calibrated parameter set** — not system identification from patient-level data.  Every numerical entry in the matrices is derived from published aggregate effect sizes (trial results, meta-analyses, Mendelian randomisation studies) and converted into the HDR state-space formulation.  Source annotations and conversion arithmetic are documented in `hdr_athero/model/sources.py`.

## What This Is Not

- Not validated for any clinical use
- Not identified from individual patient data
- Not a substitute for formal system identification (HDR Phase 2)

## Quick Start

```bash
# Install (with dev dependencies for testing)
pip install -e ".[dev]"

# Run the self-contained simulation
python -m hdr_athero.simulation.run

# Run tests
python -m pytest tests/ -v

# Run spectral analysis only
python -c "from hdr_athero.analysis.spectral import print_spectral_report; print_spectral_report()"
```

## Results Summary

### Test Suite: 24/24 Passed

| Category | Tests | Status |
|----------|-------|--------|
| Parameter shapes | 6 | All passed |
| Parameter properties | 7 | All passed |
| Spectral stability | 3 | All passed |
| Controllability | 2 | All passed |
| Simulation integration | 6 | All passed |

### Spectral Analysis

| Basin | ρ(A_k) | Recovery Amplification | Stability |
|-------|--------|----------------------|-----------|
| Subclinical | 0.9101 | 5.82× | Stable |
| Vulnerable | 0.9909 | 55.44× | Stable (near-unit-root) |
| Post-ACS | 0.8780 | 4.36× | Stable |

### Controller Comparison (Mean Cost, 20 episodes × 180 days)

| Controller | Subclinical | Vulnerable | Post-ACS |
|-----------|-------------|------------|----------|
| Open-loop | 3.17 | 16.72 | 0.07 |
| Guideline | 242.70 | 6746.95 | 14.37 |
| LQR (basin-aware) | 0.84 | 58.70 | 0.05 |
| LQR (pooled) | 0.73 | 55.13 | 0.05 |

Full results in [`results/`](results/) and [`docs/`](docs/).

## Repository Structure

```
HDR-Athero/
├── hdr_athero/
│   ├── model/
│   │   ├── parameters.py    # All numerical matrices with source annotations
│   │   ├── sources.py       # Literature sources and conversion arithmetic
│   │   └── slds.py          # SLDS model class
│   ├── simulation/
│   │   ├── run.py           # Main simulation entry point
│   │   └── controller.py    # 5 controllers: open-loop, guideline, aggressive, LQR
│   └── analysis/
│       └── spectral.py      # Spectral radii, controllability Gramian
├── tests/
│   └── test_athero.py       # 24 tests: shapes, properties, stability, simulation
├── docs/
│   ├── model_specification.md  # Full mathematical specification
│   ├── controllers.md          # Controller documentation and gains
│   ├── literature_sources.md   # Complete parameter provenance
│   └── testing_guide.md        # Testing architecture and guide
├── results/
│   ├── test_report.txt         # Detailed test results
│   ├── simulation_output.txt   # Full simulation output
│   └── spectral_summary.txt    # Spectral and controllability analysis
├── CLAUDE.md               # Claude Code instructions
├── pyproject.toml
└── README.md
```

## Literature Sources

Key data sources for parameterisation:

- **CTT Collaboration (2010)**: LDL-C dose–response and CHD risk reduction
- **CANTOS / Ridker et al. (2017)**: Canakinumab hsCRP reduction independent of LDL
- **Ashor et al. (2015)**: Exercise dose–response on FMD (2 MET → +1% FMD)
- **CRP CHD Genetics Collaboration (2011)**: MR evidence on CRP causality
- **GLAGOV / Nicholls et al. (2016)**: PCSK9i plaque regression quantification

Full source documentation in `hdr_athero/model/sources.py` and [`docs/literature_sources.md`](docs/literature_sources.md).

## Disclaimer

All constructions are theoretical and computational.  No clinical validity is asserted.  
The Thalion Initiative, Boston, Massachusetts, 2025.
