# Immunosuppressive Metabolites ACR Analysis Guide

## Overview
This document provides a comprehensive guide to the R analysis scripts for immunosuppressive metabolite data in acute cellular rejection (ACR) studies. All scripts are up-to-date, error-free, and ready to use.

---

## Project Structure

```
Scripts/
├── 00_source                    # Sources all utilities and setup
├── 02_setup                     # Data import and preprocessing
├── 03_tests                     # Sample-level independent group comparisons
├── 05a_PatientLevel             # Patient-level analysis with median/one-per aggregation
├── 05b_PatientLevel_SameAs03Tests  # Patient-level analysis mirroring 03_tests
├── 06_Paired_Pre2R_vs_Active2R  # 2-timepoint paired analysis
├── 06b_Paired_Pre_Active_Post   # 3-timepoint paired analysis
└── 07                          # Mixed-effects models for longitudinal data (NEW)

Utilities/
├── plot_enhanced_bars           # Enhanced bar plot function
├── plot_conc_corr              # Correlation plots
├── plot_RM                     # Repeated measures plots
└── plot_ttest_bars             # T-test bar plots
```

---

## Script Descriptions

### 00_source
**Purpose:** Sources all required packages, utilities, and runs data setup.

**Key Features:**
- Loads tidyverse, ggplot2, and statistical packages
- Sources all utility functions
- Runs 02_setup to prepare data objects

**When to use:** Run this at the start of any analysis session to ensure all dependencies are loaded.

---

### 02_setup
**Purpose:** Import raw data and create all analysis-ready datasets.

**Key Features:**
- Imports `IS_metabolites_ACR.csv`
- Creates log2 and non-log2 versions for statistics vs visualization
- Identifies patients with 2R+ episodes (n=18 patients, 133 total samples)
- Creates ACR category tables (0R, 1R, 0R+1R combined, 2R+)
- **Defines pre-2R, active-2R, and post-2R timepoint classifications**

**Output Objects:**
- `patients_with_2R`: All samples from patients who had at least one 2R episode
- `pre_2R`: Samples before any 2R episode
- `active_2R`: Samples during a 2R episode
- `post_2R`: Samples immediately after a 2R episode (before next 2R if any)
- `acr_0r`, `acr_1r`, `acr_0_1r`, `acr_2r`: Samples by ACR grade

**Important Notes:**
- Log2-transformed data used for statistical tests (normality assumption)
- Non-log2 data used for visualization and interpretation
- Post-2R classification is conservative (only immediate follow-up samples)

---

### 03_tests
**Purpose:** Sample-level independent group comparisons (0R+1R vs 2R+).

**Clinical Question:** "Are metabolite concentrations different in samples with 2R+ compared to samples without rejection?"

**Key Features:**
- Independent (unpaired) tests comparing two groups of samples
- Wilcoxon rank-sum test (non-parametric, primary method)
- Welch's t-test on log2 data (parametric, secondary)
- FDR correction for multiple testing
- Creates bar plots for significant metabolites (FDR < 0.05)
- Exploratory plots for suggestive findings (p < 0.1)

**Strengths:**
- Maximum statistical power (uses all samples)
- Simple, straightforward interpretation
- Standard approach for biomarker discovery

**Limitations:**
- Does NOT account for repeated samples from same patient
- Cannot control for patient-level confounders
- Treats each sample as independent

**When to use:**
- Initial exploratory analysis
- Biomarker discovery
- When patient-level matching is not feasible

---

### 05a_PatientLevel
**Purpose:** Patient-level analysis using median or one-representative-sample aggregation per patient.

**Clinical Question:** "Are metabolite levels different between patients with and without 2R episodes?"

**Key Features:**
- Aggregates multiple samples per patient (median or one representative)
- Wilcoxon test on non-log2 aggregated data (primary)
- Welch's t-test on log2 aggregated data (secondary)
- FDR correction
- Generates bar plots for significant results (FDR < 0.05)
- Exploratory plots for trending results (p < 0.1)

**Strengths:**
- Each patient contributes equally (addresses repeated measures)
- More conservative than sample-level analysis
- Patient-level interpretation (clinically relevant unit)

**Limitations:**
- Reduced statistical power (fewer data points)
- Information loss from aggregation
- Still an unpaired design (doesn't use within-patient changes)

**When to use:**
- When repeated measures from patients are common
- When patient-level differences are of primary interest
- As a more conservative complement to 03_tests

---

### 05b_PatientLevel_SameAs03Tests
**Purpose:** Patient-level analysis that exactly mirrors 03_tests methodology but at patient level.

**Key Features:**
- Same test sequence as 03_tests (Wilcoxon, Welch's t-test)
- Patient aggregation for fair comparison
- Directly comparable results to 03_tests

**When to use:**
- To compare sample-level vs patient-level results
- To validate findings from 03_tests at patient level

---

### 06_Paired_Pre2R_vs_Active2R
**Purpose:** Two-timepoint paired within-patient analysis comparing pre-2R vs active-2R samples.

**Clinical Question:** "Within the same patients, how do metabolite levels change when they develop acute rejection?"

**Key Features:**
- Includes only patients with BOTH pre-2R AND active-2R samples
- Paired statistical tests (paired t-test, Wilcoxon signed-rank)
- Aggregates multiple samples per timepoint using median
- Controls for patient-level confounders by design
- Creates paired difference plots and trajectory visualizations
- Median-focused plots with p-value annotations

**Strengths:**
- Controls for inter-patient variability (strongest design)
- Directly measures within-patient metabolic changes
- Addresses "what changes during rejection?" question
- Most clinically relevant for understanding rejection biology

**Limitations:**
- Smaller sample size (requires both timepoints per patient)
- Post-2R not included (focused on rejection onset only)
- Cannot assess recovery/reversibility

**When to use:**
- Primary analysis for understanding rejection-associated changes
- When sufficient patients have both pre and active samples
- For causal inference about rejection effects

**Output:**
- Paired test results with FDR correction
- Paired difference plots
- Individual patient trajectory plots
- Summary statistics for within-patient changes

---

### 06b_Paired_Pre_Active_Post
**Purpose:** Three-timepoint paired within-patient analysis including recovery/post-treatment samples.

**Clinical Question:** "Do metabolic changes during rejection reverse after treatment? Are there persistent alterations?"

**Key Features:**
- Requires pre-2R, active-2R, AND post-2R samples from same patients
- Friedman test (non-parametric repeated measures ANOVA for 3 groups)
- Pairwise comparisons with FDR correction:
  - Pre vs Active: Rejection onset
  - Active vs Post: Recovery/treatment response
  - Pre vs Post: Complete normalization?
- Three-timepoint trajectory visualizations
- Reversibility pattern classification:
  - **Full recovery:** Post returns close to pre baseline
  - **Partial recovery:** Post improves but doesn't normalize
  - **No recovery:** Post remains as altered as active
  - **No change:** Stable across all timepoints

**Strengths:**
- Addresses recovery/reversibility questions
- Identifies persistent vs transient changes
- May reveal subclinical ongoing processes
- Complete picture of rejection and recovery biology

**Limitations:**
- Smallest sample size (requires all 3 timepoints per patient)
- Post-2R confounded by treatment (steroids, immunosuppression changes)
- Post-2R timing variable (may not be standardized)
- Harder to interpret clinically (treatment effects vs biology)

**When to use:**
- When you have sufficient patients with all 3 timepoints (typically n≥3)
- To assess whether changes are reversible vs persistent
- To identify markers of incomplete recovery or ongoing injury
- For hypothesis generation about long-term consequences

**Important Considerations:**
- Post-2R samples include treatment effects, not just biological recovery
- Timing of post-2R samples may vary between patients
- Smaller sample size reduces statistical power
- Use in conjunction with 06_Paired_Pre2R_vs_Active2R, not as replacement

**Output:**
- `Results2/Paired_3Timepoint/3Timepoint_Statistical_Results.csv`: All test results
- `Results2/Paired_3Timepoint/Reversibility_Pattern_Summary.csv`: Recovery patterns
- Trajectory plots for significant metabolites
- Friedman test results (overall 3-way comparison)
- Pairwise comparison results

---

### 07
**Purpose:** Mixed-effects models for longitudinal data analysis.

**Clinical Question:** "How do metabolite levels change over time in relation to acute rejection and recovery?"

**Key Features:**
- Linear mixed-effects models to account for repeated measures
- Fixed effects for timepoints (pre, active, post) and group (2R+ vs 0R/1R)
- Random intercepts and slopes for patients
- Model selection based on AIC/BIC
- Visualization of fitted trajectories

**Strengths:**
- Accounts for all available data (balanced or unbalanced designs)
- Controls for patient-specific variability
- Can model complex correlation structures

**Limitations:**
- Requires more complex modeling and interpretation
- Assumes linearity and normality of residuals
- Computationally intensive

**When to use:**
- When analyzing data from 06b_Paired_Pre_Active_Post
- For a more nuanced understanding of longitudinal changes
- When interested in both fixed effects (population-level) and random effects (individual-level)

**Output:**
- Model summaries with fixed and random effects
- AIC/BIC values for model comparison
- Plots of observed vs fitted values, residuals

---

## When to Use Which Script?

### For Biomarker Discovery:
**Start with:** `03_tests` (sample-level, maximum power)  
**Validate with:** `05a_PatientLevel` (patient-level, more conservative)

### For Understanding Within-Patient Changes:
**2 timepoints:** `06_Paired_Pre2R_vs_Active2R` (paired 2-timepoint)  
**3 timepoints:** `06b_Paired_Pre_Active_Post` (if sufficient data)  
**Multiple timepoints:** `07` (mixed-effects models - most powerful)

### For Longitudinal Trajectories:
**Primary:** `07` (mixed-effects models with all timepoints)  
**Validate with:** `06` or `06b` (paired analysis as sanity check)

### For Patient-Level Characterization:
**Use:** `05a_PatientLevel` or `05b_PatientLevel_SameAs03Tests`

### For Recovery/Reversibility Questions:
**Simple approach:** `06b_Paired_Pre_Active_Post` (requires all 3 timepoints)  
**Advanced approach:** `07` (interaction model: ACR×Time)

### For Maximum Statistical Power with Repeated Measures:
**Use:** `07` (mixed-effects models use ALL data efficiently)

### Decision Tree:

```
Do you have repeated measures (>1 sample per patient)?
  NO  → Use 03_tests (sample-level) or 05a (patient-level aggregation)
  YES → Continue
       ↓
How many timepoints per patient?
  Exactly 2  → Use 06_Paired_Pre2R_vs_Active2R
  Exactly 3  → Use 06b_Paired_Pre_Active_Post
  Variable or >3 → Use 07 (mixed-effects models)
       ↓
Do you want to test time effects or trajectories?
  YES → Use 07 (mixed-effects with ACR×Time interaction)
  NO  → Use 06 or 06b (simpler paired tests)
       ↓
Are you interested in individual recovery patterns?
  YES → Use 06b (reversibility classification)
  NO  → Use 07 (population-level trajectories)
```

---

## Statistical Considerations

### Multiple Testing Correction
All scripts use FDR (false discovery rate) correction via `p.adjust(method = "fdr")`:
- **Primary threshold:** FDR < 0.05 (high confidence)
- **Exploratory threshold:** Unadjusted p < 0.1 (hypothesis generation)

### Log2 Transformation
- **For statistics:** Log2-transformed data (assumes log-normality)
- **For visualization:** Non-log2 data (interpretable concentrations)
- **Why:** Normalizes distributions, makes fold changes symmetric

### Sample Size Requirements
- **03_tests:** n=133 samples (43 0R, 68 1R, 22 2R+)
- **05a_PatientLevel:** n=18 patients with 2R episodes
- **06_Paired_Pre2R_vs_Active2R:** Depends on patients with both pre and active
- **06b_Paired_Pre_Active_Post:** Depends on patients with all 3 timepoints (typically smallest)
- **07:** Depends on number of observations per patient (balanced or unbalanced)

---

## Output Structure

```
Results2/
├── Sample_Level_Tests/              # From 03_tests
│   ├── Statistical_Results_*.csv
│   ├── Plots_FDR0.05/
│   └── Plots_p0.1_exploratory/
├── Patient_Level/                   # From 05a, 05b
│   ├── *_Statistical_Results.csv
│   └── Plots/
├── Paired_Pre2R_vs_Active2R/       # From 06
│   ├── paired_test_results.csv
│   └── paired_plots/
└── Paired_3Timepoint/              # From 06b
    ├── 3Timepoint_Statistical_Results.csv
    ├── Reversibility_Pattern_Summary.csv
    ├── Plots_Friedman_FDR0.05/
    └── Plots_Friedman_p0.1_unadjusted/
```

---

## Clinical Interpretation Guide

### 03_tests: Sample-Level Findings
**Interpretation:** "Samples collected during 2R+ episodes have different metabolite concentrations compared to samples without rejection."

**Clinical Meaning:** Potential biomarkers for detecting active rejection.

**Limitations:** Cannot tell if differences are patient-specific or universal.

---

### 05a_PatientLevel: Patient-Level Findings
**Interpretation:** "Patients who experience 2R+ episodes have different baseline/overall metabolite levels compared to those who don't."

**Clinical Meaning:** May identify patients at risk or with different metabolic phenotypes.

**Limitations:** Still doesn't show within-patient changes.

---

### 06_Paired_Pre2R_vs_Active2R: Within-Patient Changes
**Interpretation:** "Metabolite levels change within individual patients when they develop acute rejection."

**Clinical Meaning:**
- **Increase during rejection:** May reflect inflammatory response, tissue injury, or immune activation
- **Decrease during rejection:** May indicate consumption, altered clearance, or metabolic reprogramming

**Strongest Evidence:** This is the gold standard for identifying rejection-associated changes.

---

### 06b_Paired_Pre_Active_Post: Recovery/Reversibility
**Interpretation:** "Metabolite changes during rejection may or may not normalize after treatment."

**Clinical Meaning:**

**Reversible Changes (active≠pre, post≈pre):**
- Acute functional changes
- Treatment-responsive
- Likely reflects acute rejection biology

**Persistent Changes (active≠pre, post≈active):**
- May indicate permanent injury
- Incomplete recovery
- Subclinical ongoing processes
- Could predict long-term outcomes

**Late Changes (active≈pre, post≠active):**
- Treatment effects
- Post-rejection remodeling
- May not be directly related to rejection

**Important:** Post-2R samples are confounded by treatment (steroids, immunosuppression adjustments), so "persistent" doesn't necessarily mean "irreversible biology" - could reflect ongoing treatment effects.

---

## Workflow Recommendations

### Standard Analysis Pipeline

1. **Data Preparation**
   ```r
   source("Scripts/00_source")  # This also runs 02_setup
   ```

2. **Discovery Phase**
   ```r
   source("Scripts/03_tests")            # Sample-level screening
   source("Scripts/05a_PatientLevel")     # Patient-level validation
   ```

3. **Mechanistic Investigation**
   ```r
   source("Scripts/06_Paired_Pre2R_vs_Active2R")  # Within-patient changes
   ```

4. **Extended Analysis (if data permits)**
   ```r
   source("Scripts/06b_Paired_Pre_Active_Post")   # Recovery patterns
   source("Scripts/07")                           # Mixed-effects models (NEW)
   ```

### Advanced Workflow (with Mixed-Effects Models)

1. **Setup**
   ```r
   source("Scripts/00_source")
   ```

2. **Comprehensive Analysis**
   ```r
   # Discovery
   source("Scripts/03_tests")
   
   # Paired analysis (simple, interpretable)
   source("Scripts/06_Paired_Pre2R_vs_Active2R")
   
   # Mixed-effects (sophisticated, uses all data)
   source("Scripts/07")
   
   # Compare results from 06 and 07 (built into script 07)
   ```

3. **Extended Recovery Analysis (if sufficient 3-timepoint data)**
   ```r
   source("Scripts/06b_Paired_Pre_Active_Post")
   ```

### Recommended Combinations

**For Publication (Most Robust):**
```r
source("Scripts/03_tests")                      # Discovery
source("Scripts/06_Paired_Pre2R_vs_Active2R")  # Paired validation
source("Scripts/07")                           # Mixed-effects confirmation
# Report metabolites significant in ALL THREE
```

**For Maximum Power (Limited Sample Size):**
```r
source("Scripts/07")                           # Mixed-effects (uses all data)
source("Scripts/06_Paired_Pre2R_vs_Active2R")  # Validate key findings
```

**For Recovery Questions:**
```r
source("Scripts/06b_Paired_Pre_Active_Post")   # Individual recovery patterns
source("Scripts/07")                           # Population-level trajectories (ACR×Time)
```
```
