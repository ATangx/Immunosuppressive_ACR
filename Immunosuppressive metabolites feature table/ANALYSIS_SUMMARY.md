# Immunosuppressive Metabolite Analysis - Complete Script Documentation

## Project Overview
This analysis examines immunosuppressive metabolite levels across different acute cellular rejection (ACR) categories in renal transplant patients, with particular focus on post-transplant day (POD) stratification and patient-level analyses.

## Data Structure
- **Sample Categories**: pre-2R (0R/1R), active-2R (2Ra/2Rb), post-2R (subsequent samples)
- **POD Stratification**: Early (≤30d), Mid (31-90d), Late (>90d) - clinically meaningful groups
- **Patient-level vs Sample-level**: Analyses control for multiple samples per patient

---

## Analysis Scripts Description

### Core Setup and Infrastructure

#### `00_source`
**Purpose**: Central sourcing script that loads all dependencies and data preparation
- Loads required libraries with minimal output suppression
- Sources data preparation scripts (01_libraries, 02_setup) 
- Creates standardized `df_clean` dataset with consistent variable names
- Defines global variables like `metabolite_columns`

#### `01_libraries` 
**Purpose**: Library loading and package management
- Loads all required R packages for analysis and visualization
- Includes: tidyverse, ggplot2, patchwork, etc.

#### `02_setup`
**Purpose**: Data import and cleaning
- Reads the main CSV data file
- Creates standardized variable names (Patient_ID, Sample_Category, POD_group, etc.)
- Implements clinically meaningful POD stratification (Early ≤30d, Mid 31-90d, Late >90d)
- Defines metabolite columns for analysis

---

### Primary Statistical Analyses

#### `03_tests` - **Sample-Level Group Comparisons**
**Purpose**: Independent group comparisons between rejection categories
- **Analysis Type**: Sample-level (treats each sample independently)
- **Comparisons**: pre-2R vs active-2R, pre-2R vs post-2R, active-2R vs post-2R
- **Statistical Tests**: Wilcoxon rank-sum tests with FDR correction
- **Output**: Statistical results, significance tables, volcano plots
- **Use Case**: Initial screening for metabolite differences between rejection states

#### `04_PODdist` - **POD Distribution Analysis**
**Purpose**: Examine metabolite distributions across post-transplant day groups
- **Analysis Type**: Sample-level with POD stratification  
- **Groups**: Early (≤30d), Mid (31-90d), Late (>90d)
- **Statistical Tests**: ANOVA and pairwise comparisons
- **Output**: Distribution plots, statistical summaries by POD group
- **Use Case**: Understanding time-dependent changes in metabolite levels

---

### Patient-Level Analyses (Control for Multiple Samples)

#### `05a_PatientLevel` - **Basic Patient-Level Comparisons**
**Purpose**: Patient-level analysis aggregating multiple samples per patient
- **Analysis Type**: Patient-level (one value per patient per category)
- **Aggregation**: Median metabolite levels per patient per rejection category
- **Statistical Tests**: Wilcoxon tests on patient-aggregated data
- **Output**: Patient-level statistical results and comparisons
- **Use Case**: Control for patient-specific factors and multiple sampling

#### `05b_PatientLevel_SameAs03Tests` - **Patient-Level Validation**
**Purpose**: Validate that patient-level analysis gives similar results to sample-level
- **Analysis Type**: Patient-level using same logic as 03_tests
- **Comparison**: Results should be similar to 03_tests but with patient-level control
- **Output**: Comparative statistical results
- **Use Case**: Method validation and sensitivity analysis

#### `05c_PatientLevel_BaselineOnly` - **Baseline Analysis**
**Purpose**: Patient-level analysis using only baseline (first) samples per category
- **Analysis Type**: Patient-level using temporal selection
- **Selection**: First sample per patient per rejection category
- **Statistical Tests**: Wilcoxon tests on baseline samples only
- **Output**: Baseline-focused statistical results
- **Use Case**: Control for temporal effects and sampling bias

#### `05d_PatientLevel_PODstrat` - **Patient-Level POD Stratified**
**Purpose**: Combine patient-level control with POD stratification
- **Analysis Type**: Patient-level within POD strata
- **Groups**: Early (≤30d), Mid (31-90d), Late (>90d) 
- **Statistical Tests**: Separate analyses within each POD group
- **Output**: POD-stratified results with patient-level control and trajectory plots
- **Use Case**: Most comprehensive analysis controlling for both patient factors and time

---

### Specialized Stratified Analyses

#### `05_PODstrat` - **POD Stratified Sample-Level Analysis**
**Purpose**: Sample-level analysis stratified by post-transplant day groups
- **Analysis Type**: Sample-level within POD strata
- **Stratification**: Early, Mid, Late POD groups with clinically meaningful cutoffs
- **Statistical Tests**: Separate comparisons within each POD stratum
- **Output**: POD-specific statistical results
- **Use Case**: Identify time-dependent metabolite changes while maintaining sample-level granularity

#### `05_PODstrat_plots` - **POD Stratified Visualization**
**Purpose**: Comprehensive plotting for POD-stratified analyses
- **Output**: Box plots, violin plots, trajectory plots for each POD stratum
- **Visualization**: Enhanced plots with darker lines and better readability
- **Use Case**: Visual interpretation of POD-stratified results

---

### Advanced Paired Analysis

#### `06_Paired_Pre2R_vs_Active2R` - **Within-Patient Paired Analysis**
**Purpose**: Paired analysis comparing pre-2R vs active-2R within same patients
- **Analysis Type**: Matched pairs (strongest control for confounders)
- **Design**: Only includes patients with BOTH pre-2R AND active-2R samples
- **Statistical Tests**: Paired t-tests and Wilcoxon signed-rank tests
- **Output**: 
  - Paired statistical results with effect sizes
  - Before-after trajectory plots for individual patients
  - Paired difference distributions
  - Volcano plots showing effect sizes vs significance
- **Use Case**: Strongest causal inference for rejection-related metabolite changes

---

## Key Analysis Distinctions

### Sample-Level vs Patient-Level
- **Sample-Level** (03_tests, 04_PODdist, 05_PODstrat): Treats each sample as independent
- **Patient-Level** (05a-05d): Aggregates multiple samples per patient to control for patient-specific factors

### Independent vs Paired Comparisons  
- **Independent** (03_tests, 05a-05d): Compares different groups of samples/patients
- **Paired** (06_Paired): Compares changes within the same patients over time

### POD Stratification
- **Clinically Meaningful**: Early ≤30d, Mid 31-90d, Late >90d (replaced quartile-based stratification)
- **Used in**: 04_PODdist, 05d_PatientLevel_PODstrat, 05_PODstrat

---

## Output Organization

### Results2/ Directory Structure
```
Results2/
├── Sample_Level_Tests/          # 03_tests results
├── POD_Distribution/            # 04_PODdist results  
├── Patient_Level_Basic/         # 05a results
├── Patient_Level_Validation/    # 05b results
├── Patient_Level_Baseline/      # 05c results
├── Patient_Level_PODstrat/      # 05d results
├── POD_Stratified_Sample/       # 05_PODstrat results
├── POD_Stratified_Plots/        # 05_PODstrat_plots results
└── Paired_Pre2R_vs_Active2R/    # 06_Paired results
```

### Key Output Files
- **Statistical Results**: CSV files with p-values, FDR corrections, effect sizes
- **Plots**: High-resolution PNG files for publication
- **Summary Tables**: Analysis-specific summary statistics
- **Trajectory Plots**: Patient-level change visualizations

---

## Recommended Analysis Workflow

1. **Start with**: `03_tests` for initial sample-level screening
2. **Validate with**: `05a_PatientLevel` to ensure patient-level consistency  
3. **Examine timing**: `04_PODdist` or `05_PODstrat` for time-dependent effects
4. **Best control**: `05d_PatientLevel_PODstrat` for comprehensive patient+time control
5. **Strongest inference**: `06_Paired_Pre2R_vs_Active2R` for causal relationships

---

## Technical Improvements Made

### Output Suppression
- Replaced "nuclear" suppression with minimal, standard R methods
- Uses `suppressMessages()`, `capture.output()`, and temporary `cat` override
- Maintains essential output while reducing console noise

### POD Stratification  
- Removed quartile-based stratification (data-driven, potentially arbitrary)
- Implemented clinically meaningful groups: Early ≤30d, Mid 31-90d, Late >90d
- Added median split as sensitivity analysis where appropriate

### Code Standardization
- All scripts use consistent variable names and data sources
- Standardized output directories (`Results2/`)
- Consistent statistical testing approaches across analyses
- Improved plot aesthetics and readability

### Statistical Rigor
- Added effect size calculations where appropriate
- Consistent FDR multiple testing correction
- Proper paired analysis implementation
- Patient-level aggregation to control confounders

---

## Clinical Interpretation Guide

### Effect Priorities
1. **Paired Analysis** (06_Paired): Strongest evidence for rejection-related changes
2. **Patient-Level POD Stratified** (05d): Best control for confounders
3. **Patient-Level Basic** (05a): Good control, simpler interpretation
4. **Sample-Level** (03_tests): Initial screening, larger sample sizes

### Time Dependency
- Use POD-stratified analyses (04, 05d, 05_PODstrat) to understand temporal patterns
- Early POD effects may reflect acute immunological processes
- Late POD effects may reflect chronic adaptations

### Statistical Significance
- FDR < 0.05: Strong evidence after multiple testing correction
- Effect sizes: Consider clinical/biological significance beyond statistical significance
- Paired analysis effect sizes most interpretable (within-patient changes)

---

This analysis framework provides a comprehensive, statistically rigorous approach to examining immunosuppressive metabolite changes in acute cellular rejection, with appropriate controls for patient-level factors, temporal effects, and multiple testing.
