"""
Literature sources and conversion arithmetic for all numerical parameters.

Every entry in the J and B_k matrices is traced to a published source.
This module documents the provenance and the unit-conversion steps.

IMPORTANT: These are aggregate population-level effect sizes converted to
state-space parameters.  They are NOT identified from individual patient data.
All conversions involve assumptions documented inline.

Axes (n=4, standardised deviations from youthful normative reference):
    0: Δx_VI  — Vascular inflammation    (hsCRP proxy, standardised)
    1: Δx_LD  — Lipid dysregulation      (LDL-C proxy, standardised)
    2: Δx_ED  — Endothelial dysfunction  (FMD proxy, sign-inverted: higher = worse)
    3: Δx_HS  — Haemodynamic stress      (SBP proxy, standardised)

Control channels (m=4):
    0: u_1 — Lipid-lowering intensity    (standardised 0–1 scale)
    1: u_2 — Anti-inflammatory therapy   (standardised 0–1 scale)
    2: u_3 — Antihypertensive intensity  (standardised 0–1 scale)
    3: u_4 — Exercise dose               (standardised 0–1 scale)
"""

# ═══════════════════════════════════════════════════════════════════════════
# DECAY TIME CONSTANTS (τ_i in days)
# ═══════════════════════════════════════════════════════════════════════════

TAU_SOURCES = {
    "VI": {
        "basin_0": {
            "value_days": 5.0,
            "source": "hsCRP half-life ~19h (Pepys & Hirschfield, J Clin Invest 2003); "
                      "but clinical recovery of hsCRP after acute perturbation "
                      "(e.g. infection) takes 3-7 days.  We use 5 days as midpoint.",
            "confidence": "moderate",
        },
        "basin_1": {
            "value_days": 21.0,
            "source": "In chronic inflammatory states (RA, metabolic syndrome), "
                      "hsCRP remains elevated for weeks after perturbation.  "
                      "CANTOS showed sustained elevation >2 mg/L in placebo arm "
                      "over 3.7y median follow-up (Ridker et al., NEJM 2017). "
                      "We estimate 3 weeks = 21 days for the effective return "
                      "time constant in the maladaptive regime.",
            "confidence": "low — this is the hardest axis to constrain",
        },
        "basin_2": {
            "value_days": 7.0,
            "source": "Post-ACS, high-intensity statin reduces hsCRP within "
                      "days to weeks.  PROVE-IT TIMI 22 showed significant CRP "
                      "reduction by 30 days (Ridker et al., NEJM 2005). "
                      "We estimate 7 days reflecting aggressive pharmacotherapy.",
            "confidence": "moderate",
        },
    },
    "LD": {
        "basin_0": {
            "value_days": 10.0,
            "source": "LDL-C half-life in plasma ~2.5 days (Grundy, J Lipid Res 2004).  "
                      "But the HDR axis is the dysregulation state, not the particle.  "
                      "After statin initiation, steady-state LDL-C is reached in "
                      "~2 weeks.  We use 10 days as the effective return constant.",
            "confidence": "moderate-high (well-characterised pharmacokinetics)",
        },
        "basin_1": {
            "value_days": 21.0,
            "source": "In statin-resistant or undertreated states, LDL remains "
                      "persistently elevated.  Effective τ reflects metabolic "
                      "inertia and hepatic LDLR downregulation.",
            "confidence": "moderate",
        },
        "basin_2": {
            "value_days": 5.0,
            "source": "Post-ACS high-intensity statin produces rapid LDL-C "
                      "reduction (atorvastatin 80mg: ~50% reduction in 2 weeks). "
                      "IMPROVE-IT showed additive effect with ezetimibe.",
            "confidence": "moderate-high",
        },
    },
    "ED": {
        "basin_0": {
            "value_days": 2.0,
            "source": "FMD responds to acute perturbations within hours. "
                      "Atorvastatin improves endothelial function within 24h "
                      "(Laufs et al., Circulation 2001).  We use 2 days as the "
                      "effective return constant for the intact endothelium.",
            "confidence": "moderate",
        },
        "basin_1": {
            "value_days": 14.0,
            "source": "Chronically dysfunctional endothelium (low NO bioavailability, "
                      "NOS uncoupling) recovers slowly.  Exercise training studies "
                      "show FMD improvement emerges after 2-4 weeks of training "
                      "(Ashor et al., Sports Med 2015).  We use 14 days.",
            "confidence": "moderate",
        },
        "basin_2": {
            "value_days": 5.0,
            "source": "Post-ACS endothelial function is acutely impaired but "
                      "responds to pharmacotherapy within days to a week.",
            "confidence": "low — limited post-ACS FMD time-series data",
        },
    },
    "HS": {
        "basin_0": {
            "value_days": 1.0,
            "source": "Blood pressure is the most rapidly responsive axis. "
                      "Acute BP perturbations (stress, salt load) normalise "
                      "within hours via baroreflex.  We use 1 day.",
            "confidence": "high",
        },
        "basin_1": {
            "value_days": 5.0,
            "source": "In resistant hypertension with arterial stiffening, "
                      "BP recovery from perturbation is slower.  We estimate "
                      "5 days reflecting reduced vascular compliance.",
            "confidence": "moderate",
        },
        "basin_2": {
            "value_days": 2.0,
            "source": "Post-ACS patients on aggressive antihypertensive therapy.  "
                      "Beta-blockade + ACEi provides rapid haemodynamic control.",
            "confidence": "moderate",
        },
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# COUPLING MATRIX J SOURCES
# ═══════════════════════════════════════════════════════════════════════════

J_SOURCES = {
    # Row: VI (inflammation), Column: LD (lipids)
    # Inflammation ← Lipids (ox-LDL drives NF-κB, NLRP3)
    "J_VI_LD": {
        "sign": "+",
        "basin_0_value": 0.04,
        "basin_1_value": 0.10,
        "source": "ox-LDL activates endothelial NF-κB and NLRP3 inflammasome. "
                  "CANTOS showed that lipid-independent inflammation pathway "
                  "accounts for ~15% MACE reduction (Ridker et al., NEJM 2017). "
                  "CTT meta-analysis: 22% risk reduction per 1 mmol/L LDL-C "
                  "lowering over 5 years (CTT 2010).  The ratio of inflammation-"
                  "mediated vs total risk suggests J_VI_LD is moderate.",
        "conversion": "Standardised coupling coefficient.  In Basin 1, "
                      "self-reinforcing inflammatory cascade approximately "
                      "doubles the effective coupling (Libby, Nature 2021).",
        "confidence": "moderate — sign and relative magnitude well-supported; "
                      "absolute value requires longitudinal cross-lagged data",
    },

    # Row: VI, Column: ED (endothelial dysfunction)
    "J_VI_ED": {
        "sign": "+",
        "basin_0_value": 0.03,
        "basin_1_value": 0.07,
        "source": "Dysfunctional endothelium upregulates VCAM-1 and ICAM-1, "
                  "increasing monocyte transmigration (Gimbrone & García-Cardeña, "
                  "Circ Res 2016).  Effect is secondary to the lipid pathway.",
        "confidence": "moderate — mechanism clear, magnitude uncertain",
    },

    # Row: LD, Column: VI (lipids ← inflammation)
    "J_LD_VI": {
        "sign": "+",
        "basin_0_value": 0.025,
        "basin_1_value": 0.06,
        "source": "IL-6/STAT3 pathway upregulates hepatic PCSK9 expression, "
                  "reducing LDL receptor density.  MR evidence (CRP CHD "
                  "Genetics Collaboration, BMJ 2011) shows CRP itself is "
                  "NOT causal for CHD, but IL-6 pathway IS.  The coupling "
                  "is through IL-6 → PCSK9, not CRP directly. "
                  "Canakinumab did not lower LDL in CANTOS, suggesting this "
                  "coupling is modest at the IL-1β level.  We set it lower "
                  "than J_VI_LD to reflect this asymmetry.",
        "confidence": "moderate — CANTOS null LDL result constrains magnitude",
    },

    # Row: ED, Column: VI (endothelial ← inflammation)
    "J_ED_VI": {
        "sign": "+",
        "basin_0_value": 0.035,
        "basin_1_value": 0.08,
        "source": "Inflammatory cytokines reduce NO bioavailability via "
                  "eNOS uncoupling and superoxide generation (Deanfield, "
                  "Halcox & Rabelink, Circulation 2007).  This is the "
                  "primary pathway from inflammation to endothelial damage.",
        "confidence": "moderate-high — mechanism well-characterised",
    },

    # Row: ED, Column: HS (endothelial ← haemodynamic)
    "J_ED_HS": {
        "sign": "+",
        "basin_0_value": 0.03,
        "basin_1_value": 0.06,
        "source": "Elevated BP causes mechanical endothelial injury, "
                  "especially at bifurcations with disturbed flow. "
                  "BP lowering improves FMD: meta-analysis shows ~1-2% "
                  "FMD improvement per 10 mmHg SBP reduction "
                  "(Modena et al., JACC 2002).",
        "confidence": "moderate",
    },

    # Row: HS, Column: ED (haemodynamic ← endothelial)
    "J_HS_ED": {
        "sign": "+",
        "basin_0_value": 0.04,
        "basin_1_value": 0.07,
        "source": "Endothelial dysfunction reduces NO-mediated vasodilation, "
                  "increasing peripheral resistance.  Well-established in "
                  "hypertension pathophysiology (Vanhoutte et al., Br J "
                  "Pharmacol 2017).",
        "confidence": "moderate",
    },

    # Row: HS, Column: VI (haemodynamic ← inflammation)
    "J_HS_VI": {
        "sign": "+",
        "basin_0_value": 0.015,
        "basin_1_value": 0.035,
        "source": "Systemic inflammation promotes arterial stiffness and "
                  "increases vascular resistance via cytokine-mediated "
                  "smooth muscle contraction.  Effect is secondary.",
        "confidence": "low — indirect pathway, poorly quantified",
    },

    # Row: VI, Column: HS (inflammation ← haemodynamic)
    "J_VI_HS": {
        "sign": "+",
        "basin_0_value": 0.01,
        "basin_1_value": 0.025,
        "source": "Mechanical stress activates endothelial mechanosensors, "
                  "promoting NF-κB-mediated inflammatory signalling at "
                  "sites of disturbed flow (Gimbrone & García-Cardeña, "
                  "Circ Res 2016).  Weaker than direct metabolic couplings.",
        "confidence": "low — primarily a local (plaque-level) effect",
    },

    # Remaining off-diagonal entries are weak or negligible
    "J_LD_ED": {
        "sign": "~0",
        "basin_0_value": 0.005,
        "basin_1_value": 0.01,
        "source": "Weak indirect pathway.  Endothelial dysfunction does not "
                  "directly alter hepatic lipoprotein metabolism.",
        "confidence": "low",
    },

    "J_LD_HS": {
        "sign": "~0",
        "basin_0_value": 0.0,
        "basin_1_value": 0.0,
        "source": "No known direct coupling from BP to lipid metabolism.",
        "confidence": "high (absence of mechanism)",
    },

    "J_HS_LD": {
        "sign": "~0",
        "basin_0_value": 0.0,
        "basin_1_value": 0.0,
        "source": "No known direct coupling from lipids to BP.",
        "confidence": "high (absence of mechanism)",
    },

    "J_ED_LD": {
        "sign": "+",
        "basin_0_value": 0.02,
        "basin_1_value": 0.04,
        "source": "ox-LDL directly impairs endothelial function by "
                  "depleting tetrahydrobiopterin (BH4), uncoupling eNOS "
                  "(Steinberg, PNAS 1997).  Moderate effect.",
        "confidence": "moderate",
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# INTERVENTION MAP B_k SOURCES
# ═══════════════════════════════════════════════════════════════════════════

BK_SOURCES = {
    # Statin → LDL-C  (strongest, best-characterised entry)
    "B_LD_statin": {
        "value_basin_0": -0.35,
        "value_basin_1": -0.30,
        "value_basin_2": -0.40,
        "source": "CTT Collaboration (Lancet 2010): each 1 mmol/L LDL-C "
                  "reduction → 22% proportional reduction in major vascular "
                  "events over 5 years.  High-intensity statin (atorvastatin "
                  "80mg / rosuvastatin 40mg) reduces LDL-C by ~50% from "
                  "baseline (~1.8 mmol/L absolute).  'Rule of 6%': each "
                  "doubling of statin dose → additional ~6% LDL-C reduction. "
                  "Standardised to u₁ ∈ [0,1] where 1 = max tolerated dose.",
        "conversion": "At u₁=1 (max statin), Δ(LDL-C) ≈ -50%. "
                      "In standardised units (σ_LD ≈ 1 mmol/L for population), "
                      "this is approximately -1.8/σ ≈ -1.8 SD over ~14 days. "
                      "Per time step (Δt=1 day): -1.8/14 ≈ -0.13 SD/day at "
                      "full dose, but this is the steady-state effect already "
                      "captured by τ.  B_k represents the incremental control "
                      "gain beyond the decay, so we scale to ~-0.35 (see "
                      "parameters.py for full derivation).",
        "confidence": "high (best-characterised entry in the entire matrix)",
    },

    # Statin → Inflammation (pleiotropic)
    "B_VI_statin": {
        "value_basin_0": -0.08,
        "value_basin_1": -0.06,
        "value_basin_2": -0.10,
        "source": "Meta-analysis (Kinlay, Curr Cardiol Rep 2007): statins "
                  "reduce hsCRP by ~15-25% independent of LDL lowering. "
                  "JUPITER: rosuvastatin reduced hsCRP by 37% (Ridker et al., "
                  "NEJM 2008).  But JACC meta-regression (Robinson et al., "
                  "JACC 2005) found pleiotropic effects do not contribute "
                  "additional CHD risk reduction beyond LDL lowering. "
                  "We set this at ~25% of the lipid effect to reflect the "
                  "observed CRP reduction while acknowledging the uncertain "
                  "causal contribution.",
        "confidence": "moderate — CRP reduction is real, causal impact debated",
    },

    # Anti-inflammatory → Inflammation
    "B_VI_antiinfl": {
        "value_basin_0": -0.25,
        "value_basin_1": -0.18,
        "value_basin_2": -0.22,
        "source": "CANTOS (Ridker et al., NEJM 2017): canakinumab 150mg "
                  "reduced hsCRP by ~40% from baseline (median 4.2 → ~2.5 "
                  "mg/L at 48 months), producing 15% MACE reduction.  "
                  "Those achieving hsCRP <2 mg/L had 25% MACE reduction "
                  "(Ridker et al., Lancet 2018).  LoDoCo2 (Nidorf et al., "
                  "NEJM 2020): colchicine 0.5mg daily reduced hsCRP by ~30% "
                  "and MACE by 31%.  Basin 1 value is lower reflecting "
                  "attenuated drug effect against self-sustaining cascade.",
        "confidence": "moderate-high (CANTOS is the cleanest data source)",
    },

    # Antihypertensive → Haemodynamic stress
    "B_HS_antihyp": {
        "value_basin_0": -0.30,
        "value_basin_1": -0.25,
        "value_basin_2": -0.35,
        "source": "Meta-analysis (Law et al., BMJ 2009): one standard-dose "
                  "antihypertensive reduces SBP by ~9.1 mmHg; two drugs "
                  "reduce by ~13.4 mmHg.  SPRINT: intensive target (<120 "
                  "mmHg) vs standard (<140) reduced MACE by 25%.  "
                  "Population σ_SBP ≈ 15 mmHg.  At max dose combination "
                  "(u₃=1): expected ~20 mmHg reduction = ~1.3 SD.",
        "confidence": "high (extensive trial data)",
    },

    # Exercise → Endothelial function
    "B_ED_exercise": {
        "value_basin_0": -0.30,
        "value_basin_1": -0.20,
        "value_basin_2": -0.15,
        "source": "Ashor et al. (Sports Med 2015): 2 MET increase in exercise "
                  "intensity → 1 percentage-point improvement in FMD. "
                  "All modalities improve FMD significantly: aerobic WMD "
                  "2.79% (95% CI 2.12–3.45), resistance 2.52%, combined "
                  "2.07%.  Population σ_FMD ≈ 3-4%.  At u₄=1 (max exercise "
                  "programme ~20 MET-h/wk): expected ~3% FMD improvement "
                  "≈ 0.8-1.0 SD.  Basin 2 value lower due to post-ACS "
                  "exercise restrictions.",
        "confidence": "moderate-high (meta-analytic dose–response available)",
    },

    # Exercise → Haemodynamic stress
    "B_HS_exercise": {
        "value_basin_0": -0.12,
        "value_basin_1": -0.10,
        "value_basin_2": -0.08,
        "source": "Meta-analysis (Cornelissen & Smart, JACC 2013): aerobic "
                  "exercise reduces SBP by ~3-4 mmHg on average. Resistance "
                  "exercise ~2-3 mmHg.  Smaller effect than pharmacotherapy.",
        "confidence": "moderate",
    },

    # Exercise → Inflammation
    "B_VI_exercise": {
        "value_basin_0": -0.10,
        "value_basin_1": -0.07,
        "value_basin_2": -0.05,
        "source": "Meta-analysis (Fedewa et al., Med Sci Sports Exerc 2017): "
                  "exercise training reduces CRP by ~0.5 mg/L.  Effect is "
                  "smaller than pharmacological anti-inflammatory therapy.",
        "confidence": "moderate",
    },

    # Antihypertensive → Endothelial function
    "B_ED_antihyp": {
        "value_basin_0": -0.10,
        "value_basin_1": -0.08,
        "value_basin_2": -0.10,
        "source": "ACEi/ARB improve endothelial function through bradykinin-"
                  "mediated NO release.  Meta-analysis (Shahin et al., "
                  "Vascular 2011): RAS inhibitors improve FMD by ~2% absolute.",
        "confidence": "moderate",
    },
}
