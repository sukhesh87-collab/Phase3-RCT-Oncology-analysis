# Phase 3 RCT Analysis: Darbepoetin Alfa in Extensive-Stage Small Cell Lung Cancer

**Secondary analysis of a real Phase 3, randomized, double-blind, placebo-controlled oncology trial demonstrating clinical trial biostatistics methods consistent with pharma/CRO industry standards.**

**Author**: Sukhesh Sudan, MPH, BDS
**Status**: Complete — June 2026
**Trial**: NCT00119613 / Amgen Study 20010145
**Data Source**: Project Data Sphere (projectdatasphere.org)

---

## Overview

This project demonstrates an end-to-end Phase 3 clinical trial biostatistics workflow using real, de-identified individual patient data from a completed Amgen oncology trial. The drug under study — darbepoetin alfa (Aranesp/NESP), an erythropoiesis-stimulating agent (ESA) — was evaluated for its effect on overall survival in 479 patients with previously untreated extensive-stage small cell lung cancer (ES-SCLC) receiving platinum/etoposide chemotherapy.

The analysis covers the full clinical trial analytics pipeline: CDISC-consistent dataset construction (ADSL, ADTTE), ITT/Per-Protocol population derivation, Kaplan-Meier survival analysis, Cox proportional hazards regression, patient-level safety summarization (TEAEs), and longitudinal quality of life assessment (FACT-F, FACT-G, HSI, EQ-5D VAS).

**Skills demonstrated:**
- CDISC ADaM-consistent dataset construction (ADSL, ADTTE) from pre-ADaM CRT-style source data
- ITT, Safety, and Per-Protocol population derivation per ICH E9 guidance
- Kaplan-Meier estimation, log-rank test, Cox proportional hazards model in R
- Pre-specified sensitivity analysis (Per-Protocol population)
- Pre-specified subgroup analysis with forest plot
- Patient-level TEAE safety summarization (preferred term, body system class, CV subcategory)
- Longitudinal QoL analysis from wide-format PRO datasets (FACT-F, FACT-G, HSI, EQ-5D VAS)
- Clinical Study Results Report in CSR-style format
- Pre-specified Study Protocol and Statistical Analysis Plan (SAP) with TFL shells

---

## Key Results

### Primary Efficacy: Overall Survival

| Analysis | N | Events | HR (95% CI) | p-value |
|---|---|---|---|---|
| **Primary (ITT, Cox PH)** | **479** | **397** | **1.122 (0.921–1.367)** | **0.2534** |
| Sensitivity (Per-Protocol) | 426 | 349 | 1.050 (0.851–1.296) | 0.6485 |

Median OS: **9.1 months** (darbepoetin) vs **9.0 months** (placebo). No statistically significant OS benefit was demonstrated.

### Key Safety Finding

| Cardiac/CV Subcategory | Darbepoetin (n=239) | Placebo (n=240) |
|---|---|---|
| Any CV TEAE | 55 (23.0%) | 34 (14.2%) |
| Embolism / Thrombosis | 23 (9.6%) | 11 (4.6%) |
| Cerebrovascular Accident | 11 (4.6%) | 6 (2.5%) |
| MI / Coronary Artery Disorders | 8 (3.3%) | 3 (1.3%) |

Despite near 3-fold reduction in anemia incidence (12.1% vs 30.4%), confirming pharmacodynamic activity, darbepoetin was associated with a clinically meaningful excess of thromboembolic and cardiovascular events — consistent with the ESA class effect that informed the FDA's 2008 black box warning for ESAs in oncology.

### Quality of Life: Null Finding

FACT-F fatigue scores were similar between arms at all timepoints (end-of-study mean: 32.3 vs 31.4). No between-arm difference reached the MCID of 3–4 points for FACT-F. Hemoglobin Symptom Index (HSI) scores were identical (0.7 vs 0.6–0.7 throughout), despite confirmed hemoglobin correction — demonstrating dissociation between surrogate pharmacological activity and patient-reported outcomes.

---

## Clinical Context

This trial was part of the evidence base that informed the **FDA's 2008 black box warning** for erythropoiesis-stimulating agents (ESAs) in oncology. Multiple trials, including this one, showed that ESAs raised hemoglobin (a surrogate endpoint) without improving survival and with a consistent thromboembolic safety signal — a regulatory landmark case in surrogate endpoint failure. Being able to discuss this clinical and regulatory context, alongside the statistical methods used, is the intended portfolio signal from this project.

---

## Repository Contents

### Documents
| File | Description |
|---|---|
| `Project2_Study_Protocol.docx` | Pre-specified study protocol: PICO, design, dataset structure, ADSL/ADTTE build plans, analysis populations |
| `Project2_SAP.docx` | Statistical Analysis Plan: methods, censoring rules, TFL shells, multiplicity considerations |
| `Project2_Results_Report.docx` | Clinical study results report (CSR-style): Abstract, Methods, Results, Discussion, Conclusions |

### Analysis Scripts (R) — Run in Order
| Order | File | Description |
|---|---|---|
| 1 | `Project2_01_ADSL.R` | ADSL construction: load all SAS datasets via haven, derive ITT/Safety/PP population flags, merge baseline demographics and OS endpoint variables into one-row-per-patient subject-level dataset |
| 2 | `Project2_02_ADTTE_OS_Analysis.R` | ADTTE construction, Kaplan-Meier estimation, log-rank test, Cox PH primary OS analysis, PH assumption check (Schoenfeld residuals), Table 1 baseline characteristics, Table 2 OS results, pre-specified subgroup forest plot |
| 3 | `Project2_03_PPFL_Sensitivity.R` | Per-Protocol flag derivation (≥2 chemotherapy cycles completed + no major protocol deviation per C_DISP) and PP sensitivity Cox PH analysis |
| 4 | `Project2_04_AE_Final.R` | TEAE safety analysis: AE summary by severity and SAE status, most frequent TEAEs (≥5% in either arm), body system organ class summary, cardiac/CV subcategory detail |
| 5 | `Project2_05_QOL_Hgb.R` | QoL analysis from wide-format A_QOL dataset using pivot_longer (FACT-F, FACT-G, HSI, EQ-5D VAS, change from baseline); hemoglobin over time from C_HEMAT |

### Output Figures
| File | Description |
|---|---|
| `Figure1_KM_Overall_Survival.png` | Kaplan-Meier curves by treatment arm with number-at-risk table (ITT) |
| `Figure2_Forest_Subgroup.png` | Forest plot: OS subgroup analysis (ECOG, age, sex) |
| `Figure4_FACTF_Fatigue.png` | FACT-F fatigue score over time by treatment arm |
| `Figure5_EQ5D_VAS.png` | EQ-5D VAS score over time by treatment arm |

### Output Tables (Aggregate Summary Statistics)
| File | Description |
|---|---|
| `Table1_Baseline.csv` | Baseline characteristics by arm (ITT, N=479) |
| `Table2_OS_Results.csv` | Primary OS efficacy results |
| `Table3a_AE_Summary.csv` | TEAE summary by category |
| `Table3b_AE_Frequent.csv` | Most frequent TEAEs (≥5% in either arm) |
| `Table3c_AE_BodySystem.csv` | TEAEs by body system organ class |
| `Table3d_AE_Cardiac.csv` | Cardiac/CV subcategory detail |
| `Table5_Hemoglobin.csv` | Mean hemoglobin over time by arm |
| `Table6a_FACTF.csv` | FACT-F scores by arm and timepoint |
| `Table6b_FACTG.csv` | FACT-G scores by arm and timepoint |
| `Table6c_HSI.csv` | Hemoglobin Symptom Index by arm and timepoint |
| `Table6d_VAS.csv` | EQ-5D VAS by arm and timepoint |
| `Table6e_EQ5D.csv` | EQ-5D utility (baseline vs end of study) |

---

## Data Access

**Raw data and patient-level derived datasets are not included in this repository.**

This analysis used de-identified individual patient data from Project Data Sphere (projectdatasphere.org), obtained under a data use agreement (DUA) for non-commercial methodological research. Per the DUA, raw datasets and patient-level derived files (ADSL, ADTTE) cannot be redistributed publicly.

**To replicate this analysis:**
1. Register at projectdatasphere.org (free, requires professional/institutional affiliation)
2. Request access to study: Amgen Study 20010145 / NCT00119613
3. Download SAS (.sas7bdat) datasets
4. Run scripts in order: `Project2_01_ADSL.R` → `Project2_02_ADTTE_OS_Analysis.R` → `Project2_03_PPFL_Sensitivity.R` → `Project2_04_AE_Final.R` → `Project2_05_QOL_Hgb.R`
5. Set your working directory to the folder containing downloaded SAS files before running any script

All scripts use relative paths and require R v4.3.3+ with packages: `haven`, `dplyr`, `tidyr`, `survival`, `survminer`, `ggplot2`, `gtsummary`, `janitor`, `lubridate`

---

## Methods Summary

| Element | Specification |
|---|---|
| Trial design | Phase 3, randomized, double-blind, placebo-controlled |
| Drug | Darbepoetin alfa (NESP) vs Placebo |
| Chemotherapy | Platinum + etoposide (both arms) |
| Primary endpoint | Overall Survival (OS) — time from randomization to death |
| Analysis population | Intent-to-Treat (ITT): N=479 (all patients with available data) |
| Safety population | N=239 darb, N=240 placebo (≥1 dose received) |
| Per-Protocol population | N=426 (≥2 chemotherapy cycles, no major protocol deviation) |
| Primary analysis | Kaplan-Meier + log-rank test + Cox PH (two-sided α=0.05) |
| Sensitivity analysis | Per-Protocol Cox PH |
| Subgroup analysis | Pre-specified: ECOG 0-1 vs 2+; age <65 vs ≥65; sex |
| Safety analysis | Patient-level TEAE incidence by severity, SAE, drug relationship, body system |
| QoL analysis | Descriptive means by arm and timepoint (FACT-F, FACT-G, HSI, EQ-5D VAS) |
| Software | R v4.3.3; set.seed(2026) for reproducibility |
| Standards followed | ICH E9 (population definitions), CDISC ADaM (dataset conventions), CONSORT (reporting) |

---

## Disclaimer

This is a secondary re-analysis conducted for methodological demonstration purposes. Results replicate the direction of findings from the original trial but may not identically reproduce original Amgen analyses due to partial data availability (N=479 of 600 enrolled). No novel clinical conclusions should be drawn from this analysis. Raw data were not generated or modified; only analytical derivations and summaries were produced.

---

## Related Project

**Project 1**: [Comparative Effectiveness Study — SGLT2i vs DPP4i and Heart Failure Hospitalization](https://github.com/sukhesh87-collab/rwe-comparative-effectiveness-study)
A new-user, active-comparator RWE study using synthetic patient data demonstrating propensity score matching and pharmacoepidemiology methods.
