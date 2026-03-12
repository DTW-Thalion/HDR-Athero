#!/bin/bash
# ─────────────────────────────────────────────────────
# HDR-Athero: Git + GitHub Setup Script
#
# Run this after extracting the tarball to create the
# GitHub repository and push the initial commit.
# ─────────────────────────────────────────────────────

set -e

echo "═══════════════════════════════════════════════════"
echo "  HDR-Athero: Repository Setup"
echo "═══════════════════════════════════════════════════"
echo ""

# ── Step 1: Create GitHub repo ──
echo "Step 1: Creating GitHub repository..."
echo "  Run this in your terminal (requires gh CLI):"
echo ""
echo "    gh repo create HDR-Athero --public --description \\"
echo "      'Literature-calibrated 4-axis atherosclerosis SLDS for HDR v5.4'"
echo ""

# ── Step 2: Initialize and push ──
echo "Step 2: Initialize Git and push..."
echo "  Run these commands from inside the HDR-Athero directory:"
echo ""
echo "    cd HDR-Athero"
echo "    git init"
echo "    git add ."
echo "    git commit -m 'Initial commit: 4-axis atherosclerosis SLDS'"
echo "    git branch -M main"
echo "    git remote add origin git@github.com:<YOUR_USERNAME>/HDR-Athero.git"
echo "    git push -u origin main"
echo ""

# ── Step 3: Verify ──
echo "Step 3: Verify the installation:"
echo ""
echo "    pip install -e '.[dev]'"
echo "    pytest tests/ -v"
echo "    python -m hdr_athero.simulation.run"
echo ""

# ── Claude Code usage ──
echo "═══════════════════════════════════════════════════"
echo "  Claude Code: Recommended Commands"
echo "═══════════════════════════════════════════════════"
echo ""
echo "  # Open the project in Claude Code"
echo "  claude --project HDR-Athero"
echo ""
echo "  # Or from within the project directory:"
echo "  cd HDR-Athero && claude"
echo ""
echo "  # Claude Code will read CLAUDE.md automatically"
echo "  # for project context and conventions."
echo ""
echo "  # Example Claude Code prompts:"
echo '  #   "Run the simulation and explain the spectral analysis"'
echo '  #   "Add a new control channel for SGLT2 inhibitors"'
echo '  #   "Increase Basin 1 τ_VI to 28 days and re-run spectral check"'
echo '  #   "Add a Mode B committor calculation"'
echo ""
