# Quick Reference: R Analysis Scripts

## Script Quick Selector

| **Analysis Goal** | **Use This Script** | **Sample Size** | **Design** |
|-------------------|---------------------|-----------------|------------|
| Find biomarkers (discovery) | `03_tests` | n=133 samples | Independent |
| Validate at patient level | `05a_PatientLevel` | n=18 patients | Independent |
| Within-patient changes (rejection) | `06_Paired_Pre2R_vs_Active2R` | Variable | Paired (2-timepoint) |
| Recovery/reversibility | `06b_Paired_Pre_Active_Post` | Variable | Paired (3-timepoint) |
| Longitudinal trajectories | `07` | All patients | Mixed-effects (NEW) |

---

## One-Line Script Descriptions

| **Script** | **What It Does** |
|------------|------------------|
| `00_source` | Loads everything you need to start |
| `02_setup` | Prepares data for analysis |
| `03_tests` | Compares samples with vs without rejection |
| `05a_PatientLevel` | Compares patients with median aggregation |
| `05b_PatientLevel_SameAs03Tests` | Same as 03 but at patient level |
| `06_Paired_Pre2R_vs_Active2R` | Within-patient pre→active rejection changes |
| `06b_Paired_Pre_Active_Post` | Within-patient pre→active→post recovery patterns |
| `07` | Mixed-effects models for longitudinal trajectories (NEW) |

---

## Clinical Questions Answered

| **Question** | **Script** |
|--------------|------------|
| "What metabolites differ in rejection samples?" | `03_tests` |
| "Do patients with rejection have different metabolite profiles?" | `05a_PatientLevel` |
| "How do metabolites change when a patient develops rejection?" | `06_Paired_Pre2R_vs_Active2R` |
| "Do metabolite changes reverse after treatment?" | `06b_Paired_Pre_Active_Post` |
| "How do metabolite trajectories differ over time with rejection?" | `07` (NEW) |
| "Do rejection patients have different metabolite slopes?" | `07` (ACR×Time interaction) |

---

## Statistical Tests by Script

| **Script** | **Primary Test** | **Secondary Test** |
|------------|------------------|-------------------|
| `03_tests` | Wilcoxon rank-sum | Welch's t-test (log2) |
| `05a_PatientLevel` | Wilcoxon (non-log2) | Welch's t-test (log2) |
| `06_Paired_Pre2R_vs_Active2R` | Wilcoxon signed-rank (paired) | Paired t-test |
| `06b_Paired_Pre_Active_Post` | Friedman test | Pairwise Wilcoxon |
| `07` | Linear mixed-effects (LME) | Model comparison (AIC) |

---

## Data Objects Created by 02_setup

| **Object** | **Description** | **N Samples** |
|------------|-----------------|---------------|
| `patients_with_2R` | All samples from patients with ≥1 2R episode | 133 |
| `pre_2R` | Samples before any 2R episode | Variable |
| `active_2R` | Samples during a 2R episode | 22 |
| `post_2R` | Samples immediately after 2R episode | Variable |
| `acr_0r` | All 0R samples | 43 |
| `acr_1r` | All 1R samples | 68 |
| `acr_0_1r` | Combined 0R+1R samples | 111 |
| `acr_2r` | All 2R+ samples | 22 |

---

## Typical Workflow

```r
# 1. Setup (run once per session)
source("Scripts/00_source")

# 2. Discovery & validation
source("Scripts/03_tests")            # Sample-level
source("Scripts/05a_PatientLevel")    # Patient-level

# 3. Mechanistic investigation
source("Scripts/06_Paired_Pre2R_vs_Active2R")  # Paired 2-timepoint

# 4. Extended analysis (optional, if data permits)
source("Scripts/06b_Paired_Pre_Active_Post")   # Paired 3-timepoint
source("Scripts/07")                           # Mixed-effects models (NEW)
```

## Advanced Workflow (Maximum Power)

```r
# 1. Setup
source("Scripts/00_source")

# 2. Comprehensive paired + mixed-effects analysis
source("Scripts/06_Paired_Pre2R_vs_Active2R")  # Simple paired (interpretable)
source("Scripts/07")                           # Mixed-effects (sophisticated)
# Script 07 automatically compares with script 06 results

# 3. Report metabolites significant in BOTH (most robust)
```

---

## Output Interpretation

### FDR < 0.05
✅ **High confidence finding** - Report in main results

### p < 0.05 (unadjusted)
⚠️ **Nominal significance** - Use cautiously, validate

### p < 0.1 (unadjusted)
🔍 **Exploratory/suggestive** - Hypothesis generation only

---

## Decision Tree

```
START
  ↓
Do you want to discover biomarkers?
  YES → Run 03_tests, then validate with 05a_PatientLevel
  NO  → Continue
  ↓
Do you want to understand within-patient changes?
  YES → Run 06_Paired_Pre2R_vs_Active2R
  NO  → Continue
  ↓
Do you have patients with all 3 timepoints (pre, active, post)?
  YES → Run 06b_Paired_Pre_Active_Post (recovery analysis)
  NO  → Analysis complete - review results
  ↓
Do you want to analyze longitudinal trajectories?
  YES → Run 07 (mixed-effects analysis)
  NO  → Analysis complete - review results
```

---

## Common Pitfalls

❌ **Don't:** Use 03_tests results alone (inflated type I error from repeated measures)  
✅ **Do:** Validate with 05a_PatientLevel

❌ **Don't:** Interpret post-2R as "biological recovery" (confounded by treatment)  
✅ **Do:** Recognize treatment effects in post-2R samples

❌ **Don't:** Expect large sample size for 06b (requires all 3 timepoints)  
✅ **Do:** Use 06 as primary paired analysis, 06b as hypothesis-generating

❌ **Don't:** Ignore exploratory plots (p < 0.1 may reveal biology)  
✅ **Do:** Use exploratory findings to guide follow-up studies

❌ **Don't:** Overinterpret mixed-effects results (focus on fixed effects)  
✅ **Do:** Use mixed-effects models to explore individual trajectories

---

## File Locations

**Scripts:** `Immunosuppressive metabolites feature table/Scripts/`  
**Data:** `Immunosuppressive metabolites feature table/IS_metabolites_ACR.csv`  
**Results:** `Results2/` (auto-created by scripts)  
**Guide:** `Immunosuppressive metabolites feature table/ANALYSIS_GUIDE.md` (detailed documentation)

---

## Status: ✅ ALL SCRIPTS READY TO USE

Last validated: All scripts error-free and outputs confirmed.

---

**For detailed information, see `ANALYSIS_GUIDE.md`**
