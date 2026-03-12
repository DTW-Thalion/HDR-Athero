# CLAUDE.md — HDR-Athero

## Project Overview

Literature-calibrated 4-axis atherosclerosis SLDS model for HDR v5.4.
This is a computational companion — not a clinical tool.

## Quick Commands

```bash
# Install in dev mode
pip install -e ".[dev]"

# Run all tests
pytest tests/ -v

# Run the full simulation suite
python -m hdr_athero.simulation.run

# Run just the spectral analysis
python -c "from hdr_athero.analysis.spectral import print_spectral_report; print_spectral_report()"

# Run the latent space recovery experiment
python -c "from hdr_athero.analysis.latent_recovery import run_experiment; run_experiment()"
```

## Structure

- `hdr_athero/model/parameters.py` — All numerical matrices (A_k, B_k, J, Q, tau)
- `hdr_athero/model/sources.py` — Literature provenance for every parameter
- `hdr_athero/model/slds.py` — SLDS model class with step/simulate
- `hdr_athero/simulation/controller.py` — Open-loop, guideline, LQR controllers
- `hdr_athero/simulation/run.py` — Main simulation entry point
- `hdr_athero/analysis/spectral.py` — Spectral radii, controllability Gramian
- `hdr_athero/analysis/observation_model.py` — Observation matrix C, noise R, missingness generator
- `hdr_athero/analysis/latent_recovery.py` — PCA / AE / VAE latent space recovery experiment

## Key Design Decisions

- n=4 state dim (VI, LD, ED, HS); reduced from 8-axis parent HDR
- m=4 control dim (statin, anti-inflammatory, antihypertensive, exercise)
- K=3 basins (subclinical, vulnerable, post-ACS)
- Δt = 1 day
- All B_k entries ≤ 0 (therapeutic sign convention)
- All J off-diagonal entries ≥ 0 (disease is self-reinforcing)
- Parameters from published trial data — see sources.py

## Dependencies

Core: numpy, scipy
Dev: pytest, matplotlib (optional)
Python ≥ 3.10
