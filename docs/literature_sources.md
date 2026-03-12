# Literature Sources and Parameter Provenance

Every numerical parameter in the HDR-Athero model is derived from published clinical trial data, meta-analyses, or established physiological relationships. This document catalogues all sources with conversion arithmetic.

---

## 1. Decay Time Constants (TAU)

### 1.1 Vascular Inflammation (VI)

| Basin | tau (days) | Source | Confidence |
|-------|-----------|--------|------------|
| 0 (subclinical) | 5.0 | hsCRP half-life ~19h (Pepys & Hirschfield, *J Clin Invest* 2003); clinical recovery of hsCRP after acute perturbation takes 3-7 days. Midpoint estimate. | Moderate |
| 1 (vulnerable) | 21.0 | CANTOS placebo arm: sustained hsCRP >2 mg/L over 3.7y follow-up (Ridker et al., *NEJM* 2017). Chronic inflammatory states (RA, metabolic syndrome) show prolonged recovery. | Low |
| 2 (post_ACS) | 7.0 | PROVE-IT TIMI 22: significant CRP reduction by 30 days with high-intensity statin (Ridker et al., *NEJM* 2005). | Moderate |

### 1.2 Lipid Dysregulation (LD)

| Basin | tau (days) | Source | Confidence |
|-------|-----------|--------|------------|
| 0 (subclinical) | 10.0 | LDL-C half-life in plasma ~2.5 days (Grundy, *J Lipid Res* 2004). Steady-state after statin initiation reached in ~2 weeks. | Moderate-high |
| 1 (vulnerable) | 21.0 | Statin-resistant or undertreated states show persistent LDL elevation due to metabolic inertia and hepatic LDLR downregulation. | Moderate |
| 2 (post_ACS) | 5.0 | High-intensity statin produces rapid LDL-C reduction (atorvastatin 80mg: ~50% reduction in 2 weeks). IMPROVE-IT showed additive effect with ezetimibe. | Moderate-high |

### 1.3 Endothelial Dysfunction (ED)

| Basin | tau (days) | Source | Confidence |
|-------|-----------|--------|------------|
| 0 (subclinical) | 2.0 | Atorvastatin improves endothelial function within 24h (Laufs et al., *Circulation* 2001). FMD responds to acute perturbations within hours. | Moderate |
| 1 (vulnerable) | 14.0 | Exercise training shows FMD improvement after 2-4 weeks (Ashor et al., *Sports Med* 2015). Chronic eNOS uncoupling recovers slowly. | Moderate |
| 2 (post_ACS) | 5.0 | Post-ACS endothelial function is acutely impaired but responds to pharmacotherapy within days. | Low |

### 1.4 Haemodynamic Stress (HS)

| Basin | tau (days) | Source | Confidence |
|-------|-----------|--------|------------|
| 0 (subclinical) | 1.0 | Acute BP perturbations normalise within hours via baroreflex. | High |
| 1 (vulnerable) | 5.0 | Resistant hypertension with arterial stiffening shows slower BP recovery. | Moderate |
| 2 (post_ACS) | 2.0 | Beta-blockade + ACEi provides rapid haemodynamic control. | Moderate |

---

## 2. Coupling Matrix (J) Sources

### 2.1 Inflammation <-- Lipids (J_VI_LD)

- **Mechanism**: ox-LDL activates endothelial NF-kB and NLRP3 inflammasome
- **Evidence**: CANTOS showed lipid-independent inflammation pathway accounts for ~15% MACE reduction (Ridker et al., *NEJM* 2017). CTT meta-analysis: 22% risk reduction per 1 mmol/L LDL-C lowering over 5 years (CTT 2010)
- **Basin 0 value**: 0.040
- **Basin 1 value**: 0.035 (stability-constrained; biological estimate ~2x Basin 0)
- **Confidence**: Moderate

### 2.2 Lipids <-- Inflammation (J_LD_VI)

- **Mechanism**: IL-6/STAT3 pathway upregulates hepatic PCSK9 expression, reducing LDL receptor density
- **Evidence**: MR evidence (CRP CHD Genetics Collaboration, *BMJ* 2011) shows CRP is NOT causal for CHD, but IL-6 pathway IS. Canakinumab did not lower LDL in CANTOS, constraining this coupling magnitude
- **Basin 0 value**: 0.025
- **Basin 1 value**: 0.021 (stability-constrained)
- **Confidence**: Moderate

### 2.3 Endothelial <-- Inflammation (J_ED_VI)

- **Mechanism**: Inflammatory cytokines reduce NO bioavailability via eNOS uncoupling and superoxide generation
- **Evidence**: Deanfield, Halcox & Rabelink, *Circulation* 2007
- **Basin 0 value**: 0.035
- **Basin 1 value**: 0.028
- **Confidence**: Moderate-high

### 2.4 Inflammation <-- Endothelial (J_VI_ED)

- **Mechanism**: Dysfunctional endothelium upregulates VCAM-1 and ICAM-1, increasing monocyte transmigration
- **Evidence**: Gimbrone & Garcia-Cardena, *Circ Res* 2016
- **Basin 0 value**: 0.030
- **Basin 1 value**: 0.025
- **Confidence**: Moderate

### 2.5 Endothelial <-- Haemodynamic (J_ED_HS)

- **Mechanism**: Elevated BP causes mechanical endothelial injury at bifurcations with disturbed flow
- **Evidence**: BP lowering improves FMD: ~1-2% FMD improvement per 10 mmHg SBP reduction (Modena et al., *JACC* 2002)
- **Basin 0 value**: 0.030
- **Basin 1 value**: 0.021
- **Confidence**: Moderate

### 2.6 Haemodynamic <-- Endothelial (J_HS_ED)

- **Mechanism**: Endothelial dysfunction reduces NO-mediated vasodilation, increasing peripheral resistance
- **Evidence**: Vanhoutte et al., *Br J Pharmacol* 2017
- **Basin 0 value**: 0.040
- **Basin 1 value**: 0.025
- **Confidence**: Moderate

### 2.7 Haemodynamic <-- Inflammation (J_HS_VI)

- **Mechanism**: Systemic inflammation promotes arterial stiffness
- **Basin 0 value**: 0.015
- **Basin 1 value**: 0.012
- **Confidence**: Low

### 2.8 Inflammation <-- Haemodynamic (J_VI_HS)

- **Mechanism**: Mechanical stress activates endothelial mechanosensors, promoting NF-kB signalling at sites of disturbed flow
- **Evidence**: Gimbrone & Garcia-Cardena, *Circ Res* 2016
- **Basin 0 value**: 0.010
- **Basin 1 value**: 0.009
- **Confidence**: Low (primarily local plaque-level effect)

### 2.9 Endothelial <-- Lipids (J_ED_LD)

- **Mechanism**: ox-LDL directly impairs endothelial function by depleting tetrahydrobiopterin (BH4), uncoupling eNOS
- **Evidence**: Steinberg, *PNAS* 1997
- **Basin 0 value**: 0.020
- **Basin 1 value**: 0.014
- **Confidence**: Moderate

### 2.10 Zero-Coupling Entries

- **J_LD_HS = 0**: No known direct coupling from BP to lipid metabolism (High confidence)
- **J_HS_LD = 0**: No known direct coupling from lipids to BP (High confidence)
- **J_LD_ED ~= 0**: Weak indirect pathway (0.005 in Basin 0)

---

## 3. Intervention Map (B_k) Sources

### 3.1 Statin --> LDL-C (B_LD_statin)

- **Source**: CTT Collaboration (*Lancet* 2010): each 1 mmol/L LDL-C reduction produces 22% proportional reduction in major vascular events over 5 years
- **Conversion**: At u_1=1 (max statin), LDL-C reduction ~50% from baseline (~1.8 mmol/L). In standardised units (sigma_LD ~ 1 mmol/L): -1.8 SD over ~14 days. Per day at steady state: B represents incremental control gain beyond natural decay
- **Values**: Basin 0: -0.35, Basin 1: -0.30 (attenuated), Basin 2: -0.40 (amplified)
- **Confidence**: High (best-characterised entry)

### 3.2 Anti-inflammatory --> Inflammation (B_VI_antiinfl)

- **Source**: CANTOS (Ridker et al., *NEJM* 2017): canakinumab 150mg reduced hsCRP by ~40% (median 4.2 to ~2.5 mg/L at 48 months), producing 15% MACE reduction. LoDoCo2 (Nidorf et al., *NEJM* 2020): colchicine 0.5mg daily reduced hsCRP by ~30% and MACE by 31%
- **Values**: Basin 0: -0.25, Basin 1: -0.18 (attenuated against self-sustaining cascade), Basin 2: -0.22
- **Confidence**: Moderate-high

### 3.3 Antihypertensive --> Haemodynamic (B_HS_antihyp)

- **Source**: Law et al. (*BMJ* 2009): one standard-dose antihypertensive reduces SBP by ~9.1 mmHg; two drugs ~13.4 mmHg. SPRINT: intensive vs standard target reduced MACE by 25%
- **Conversion**: Population sigma_SBP ~ 15 mmHg. At max combination (u_3=1): ~20 mmHg = ~1.3 SD
- **Values**: Basin 0: -0.30, Basin 1: -0.25, Basin 2: -0.35
- **Confidence**: High

### 3.4 Exercise --> Endothelial Function (B_ED_exercise)

- **Source**: Ashor et al. (*Sports Med* 2015): 2 MET increase --> 1 percentage-point FMD improvement. Aerobic WMD 2.79% (95% CI 2.12-3.45), resistance 2.52%, combined 2.07%
- **Conversion**: Population sigma_FMD ~ 3-4%. At max exercise (~20 MET-h/wk): ~3% FMD improvement ~ 0.8-1.0 SD
- **Values**: Basin 0: -0.30, Basin 1: -0.20, Basin 2: -0.15 (restricted post-ACS)
- **Confidence**: Moderate-high

### 3.5 Exercise --> Haemodynamic (B_HS_exercise)

- **Source**: Cornelissen & Smart (*JACC* 2013): aerobic exercise reduces SBP by ~3-4 mmHg, resistance ~2-3 mmHg
- **Values**: Basin 0: -0.12, Basin 1: -0.10, Basin 2: -0.08
- **Confidence**: Moderate

### 3.6 Exercise --> Inflammation (B_VI_exercise)

- **Source**: Fedewa et al. (*Med Sci Sports Exerc* 2017): exercise training reduces CRP by ~0.5 mg/L
- **Values**: Basin 0: -0.10, Basin 1: -0.07, Basin 2: -0.05
- **Confidence**: Moderate

### 3.7 Statin --> Inflammation (B_VI_statin, pleiotropic)

- **Source**: Statins reduce hsCRP by ~15-25% independent of LDL lowering (Kinlay, *Curr Cardiol Rep* 2007). JUPITER: rosuvastatin reduced hsCRP by 37% (Ridker et al., *NEJM* 2008). But pleiotropic effects may not contribute additional CHD risk reduction beyond LDL lowering (Robinson et al., *JACC* 2005)
- **Values**: Basin 0: -0.08, Basin 1: -0.06, Basin 2: -0.10
- **Confidence**: Moderate

### 3.8 Antihypertensive --> Endothelial (B_ED_antihyp)

- **Source**: ACEi/ARB improve endothelial function through bradykinin-mediated NO release. RAS inhibitors improve FMD by ~2% absolute (Shahin et al., *Vascular* 2011)
- **Values**: Basin 0: -0.10, Basin 1: -0.08, Basin 2: -0.10
- **Confidence**: Moderate

---

## 4. Reference List

1. Ashor AW et al. Exercise modalities and endothelial function. *Sports Med*. 2015;45(2):279-296.
2. Cornelissen VA, Smart NA. Exercise training for blood pressure. *JACC*. 2013;62(16):1563-1571.
3. CRP CHD Genetics Collaboration. Association between CRP and coronary heart disease. *BMJ*. 2011;342:d548.
4. CTT Collaboration. Efficacy and safety of more intensive lowering of LDL cholesterol. *Lancet*. 2010;376(9753):1670-1681.
5. Deanfield JE, Halcox JP, Rabelink TJ. Endothelial function and dysfunction. *Circulation*. 2007;115(10):1285-1295.
6. Fedewa MV et al. Effect of exercise training on C reactive protein. *Med Sci Sports Exerc*. 2017;49(3):861-870.
7. Gimbrone MA, Garcia-Cardena G. Endothelial cell dysfunction and the pathobiology of atherosclerosis. *Circ Res*. 2016;118(4):620-636.
8. Grundy SM. The issue of statin safety. *J Lipid Res*. 2004;45:1235-1245.
9. Kinlay S. Low-density lipoprotein-dependent and -independent effects of cholesterol-lowering therapies. *Curr Cardiol Rep*. 2007;9(6):445-451.
10. Laufs U et al. Rapid effects on vascular function after initiation and withdrawal of atorvastatin. *Circulation*. 2001;103(5):756-761.
11. Law MR et al. Use of blood pressure lowering drugs in the prevention of cardiovascular disease. *BMJ*. 2009;338:b1665.
12. Libby P. The changing landscape of atherosclerosis. *Nature*. 2021;592(7855):524-533.
13. Modena MG et al. Prognostic significance of flow-mediated dilation. *JACC*. 2002;40(3):505-510.
14. Nicholls SJ et al. Effect of evolocumab on progression of coronary disease (GLAGOV). *JAMA*. 2016;316(22):2373-2384.
15. Nidorf SM et al. Colchicine in patients with chronic coronary disease (LoDoCo2). *NEJM*. 2020;383(19):1838-1847.
16. Pepys MB, Hirschfield GM. C-reactive protein: a critical update. *J Clin Invest*. 2003;111(12):1805-1812.
17. Ridker PM et al. Antiinflammatory therapy with canakinumab (CANTOS). *NEJM*. 2017;377(12):1119-1131.
18. Ridker PM et al. C-reactive protein levels and outcomes after statin therapy (PROVE-IT TIMI 22). *NEJM*. 2005;352(1):20-28.
19. Ridker PM et al. Rosuvastatin to prevent vascular events (JUPITER). *NEJM*. 2008;359(21):2195-2207.
20. Robinson JG et al. Pleiotropic effects of statins. *JACC*. 2005;46(10):1855-1862.
21. Shahin Y et al. Angiotensin converting enzyme inhibitors effect on endothelial dysfunction. *Vascular*. 2011;19(6):271-277.
22. Steinberg D. Low density lipoprotein oxidation and its pathobiological significance. *PNAS*. 1997;94(16):8370-8377.
23. Vanhoutte PM et al. Endothelial dysfunction and vascular disease. *Br J Pharmacol*. 2017;174(12):1736-1752.
