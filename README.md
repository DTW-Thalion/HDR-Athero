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
# Install
pip install -e .

# Run the self-contained simulation
python -m hdr_athero.simulation.run

# Run tests
pytest tests/ -v
```

## Repository Structure

```
HDR-Athero/
├── hdr_athero/
│   ├── model/
│   │   ├── parameters.py    # All numerical matrices with source annotations
│   │   ├── sources.py       # Literature sources and conversion arithmetic
│   │   ├── basins.py        # Basin definitions and A_k construction
│   │   └── slds.py          # SLDS model class
│   ├── simulation/
│   │   ├── run.py           # Main simulation entry point
│   │   ├── controller.py    # Mode A LQR controller
│   │   └── observer.py      # Simplified IMM observer
│   └── analysis/
│       ├── spectral.py      # Spectral analysis of A_k
│       ├── controllability.py # Controllability Gramian analysis
│       └── report.py        # Generate summary report
├── tests/
│   ├── test_parameters.py   # Matrix property verification
│   ├── test_stability.py    # ρ(A_k) < 1 for all basins
│   └── test_simulation.py   # End-to-end smoke test
├── results/                 # Auto-generated outputs
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

Full source documentation in `hdr_athero/model/sources.py`.

## Disclaimer

All constructions are theoretical and computational.  No clinical validity is asserted.  
The Thalion Initiative, Boston, Massachusetts, 2025.
