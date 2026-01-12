# MPA Analysis Scripts Guide

## Overview
This guide explains the differences between the MPA analysis scripts (08, 08a, 08b, 08c) and when to use each one.

---

## Script Summaries

### **08_MPA_POD.R** - Within-Patient Comparison
**Comparison:** 0R samples vs 2R+ samples **within the same patients**

**Patient Group:** Only patients who develop 2R at some point

**Samples:**
- 0R samples from these patients
- 2R+ samples from these same patients

**Purpose:** Tests **rejection STATE effects** (what happens when rejection occurs)

**Key Question:** "Within patients who reject, how do MPA levels change from quiescence (0R) to rejection (2R+)?"

**Output Directory:** `Results2/Feedback_Analysis/Within_Patient_POD_Stratified/`

---

### **08a_MPA_POD.R** - Patient Phenotype (Broader Inclusion)
**Comparison:** ALL samples from 0R/1R-only patients vs ALL samples from 2R patients

**Patient Groups:**
- 0R/1R-only patients: Never develop 2R+ (n=XX patients)
- 2R patients: Develop 2R at some point (n=18 patients)

**Samples:**
- ALL samples from 0R/1R-only patients (includes both 0R and 1R)
- ALL samples from 2R patients (includes 0R, 1R, AND 2R+ samples)

**Purpose:** Tests **inherent patient-level differences** with broader inclusion criteria

**Key Question:** "Do patients who develop severe rejection (2R) inherently differ from patients who only have mild/no rejection (0R/1R)?"

**Output Directory:** `Results2/Feedback_Analysis/0R1R_only_vs_2R/`

**Note:** Similar to 08b but includes 1R patients in the "non-severe" group

**Note**
- Crude Overall Summary - 

Rejection state confounding: The 2R patient group includes:

0R samples (quiescent)
1R samples (mild rejection)
2R+ samples (severe rejection)
So you're mixing different biological states in one group

Time confounding: If 2R patients tend to have samples at different POD than 0R/1R-only patients, that could drive differences

Treatment confounding: Patients who develop 2R might get different treatments/dose adjustments
---

### **08b_MPA_POD.R** - Patient Phenotype (Strictest Inclusion)
**Comparison:** ALL samples from 0R-only patients vs ALL samples from 2R patients

**Patient Groups:**
- 0R-only patients: Never have ANY rejection (n=5 patients, strictest non-rejectors)
- 2R patients: Develop 2R at some point (n=18 patients)

**Samples:**
- ALL samples from 0R-only patients (all are 0R)
- ALL samples from 2R patients (includes 0R, 1R, AND 2R+ samples)

**Purpose:** Tests **inherent patient-level differences** with strictest non-rejector definition

**Key Question:** "Do patients who NEVER reject inherently differ from patients prone to rejection?"

**Output Directory:** `Results2/Feedback_Analysis/0R_only_vs_2R/`

**Note:** Similar to 08a but uses strictest definition of non-rejectors (0R-only)

---

### **08c_MPA_patient_phenotype.R** - Patient Phenotype (State-Controlled)
**Comparison 1:** 0R samples from 0R-only patients vs 0R samples from 2R patients

**Comparison 2:** 0R/1R samples from 0R/1R-only patients vs 0R/1R samples from 2R patients

**Patient Groups:**
- 0R-only or 0R/1R-only patients (depending on comparison)
- 2R patients

**Samples:**
- **Part 1-4:** Only 0R samples (same rejection state in both groups)
- **Part 5:** Only 0R/1R samples (same rejection state range in both groups)

**Purpose:** Tests **pure patient phenotype** while **controlling for rejection state**

**Key Questions:** 
1. "Do 0R samples from never-rejectors look different from 0R samples from future-rejectors?"
2. "Do 0R/1R samples from non-severe rejectors differ from 0R/1R samples from severe rejectors?"

**Output Directory:** `Results2/Feedback_Analysis/Patient_Phenotype/`

**Note:** This is the CLEANEST patient-level comparison because rejection state is controlled

---

## Comparison Matrix

| Script | Patient Group 1 | Patient Group 2 | Sample Selection | Tests |
|--------|----------------|-----------------|------------------|-------|
| **08** | 2R patients | 2R patients (same) | 0R vs 2R+ (within) | Rejection STATE |
| **08a** | 0R/1R-only | 2R patients | ALL vs ALL | Patient PHENOTYPE (broad) |
| **08b** | 0R-only | 2R patients | ALL vs ALL | Patient PHENOTYPE (strict) |
| **08c Part 1-4** | 0R-only | 2R patients | 0R vs 0R | Patient PHENOTYPE (state-controlled) |
| **08c Part 5** | 0R/1R-only | 2R patients | 0R/1R vs 0R/1R | Patient PHENOTYPE (state-controlled) |

---

## Clinical Interpretation

### Within-Patient (08)
- **Significant result:** MPA levels change during rejection episodes
- **Clinical use:** Understand drug behavior during acute rejection
- **Limitation:** Only includes patients who reject

### Between-Patient, All Samples (08a, 08b)
- **Significant result:** Rejection-prone patients have inherently different MPA levels
- **Clinical use:** Identify patients at risk based on baseline drug levels
- **Limitation:** Mixes rejection states (confounds phenotype with state)
- **08a vs 08b:** 08a has more power (larger n), 08b has stricter definition

### Between-Patient, State-Controlled (08c)
- **Significant result:** Patient phenotype differences independent of rejection state
- **Clinical use:** Strongest evidence for baseline risk stratification
- **Advantage:** Cleanest comparison, controls for confounding
- **0R comparison:** Strictest quiescence
- **0R/1R comparison:** More power, clinically relevant (1R often untreated)

---

## Recommended Analysis Strategy

1. **Start with 08c** (state-controlled) → Cleanest patient phenotype test
2. **Then run 08a and 08b** → See if including all samples changes conclusions
3. **Finally run 08** (within-patient) → Understand rejection state effects

### Interpretation Patterns

**If 08c significant but 08a/08b not:**
- Rejection state dilutes the signal
- Patient differences exist but only detectable in controlled comparison

**If 08a/08b significant but 08c not:**
- Differences driven by mixing rejection states
- 2R+ samples drive the effect, not inherent patient phenotype

**If all are significant:**
- Strong evidence for inherent patient-level differences
- MPA levels consistently lower/higher in rejection-prone patients

**If only 08 is significant:**
- Rejection is a STATE phenomenon, not a PHENOTYPE
- Focus on monitoring during at-risk periods

---

## POD Stratifications

All scripts include:
1. **Overall** (all POD)
2. **Early vs Late** (≤30 days vs >30 days)
3. **Excluding POD<10** (removes immediate post-transplant)
4. **POD 10-45 only** (acute window)

**Clinical Rationale:**
- Early period: Drug stabilization, protocol-driven dosing
- Late period: Individualized dosing, steady state
- POD<10: High variability, not clinically actionable
- POD 10-45: High-risk window for rejection

---

## Output Structure

```
Results2/Feedback_Analysis/
├── Within_Patient_POD_Stratified/          # 08
├── 0R1R_only_vs_2R/                        # 08a
├── 0R_only_vs_2R/                          # 08b
└── Patient_Phenotype/                      # 08c
```

Each directory contains:
- `MMF_Overall.png`
- `MMF_Early.png`
- `MMF_Late.png`
- `MMF_Early_POD10plus.png`
- `MMF_Late_POD10plus.png`
- `MMF_Early_POD10_45.png`
- `MMF_Late_POD10_45.png`

08c also includes:
- `MMF_0R1R_Samples_*.png` (Part 5 analysis)

---

## Statistical Notes

- **Test:** Wilcoxon rank-sum (Mann-Whitney U)
- **Visualization:** Boxplots with jittered points, red diamonds = median
- **Adjustment:** No multiple testing correction (exploratory analysis)
- **Sample size:** Displayed on each plot for transparency

---

## Questions? 

Compare results across scripts to understand:
1. Is rejection a STATE or PHENOTYPE phenomenon?
2. Does including 1R patients change conclusions (08a vs 08b)?
3. Are differences driven by rejection state or patient baseline (08a/08b vs 08c)?
4. Do MPA levels change during rejection episodes (08)?
