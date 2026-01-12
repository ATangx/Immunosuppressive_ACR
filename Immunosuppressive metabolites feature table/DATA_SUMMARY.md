# COMPREHENSIVE DATA SUMMARY

**Generated:** January 12, 2026  
**Dataset:** Immunosuppressive Metabolites in Heart Transplant Recipients

---

## 1. OVERALL DATASET OVERVIEW

| Metric | Count |
|--------|-------|
| **Total Samples** | 354 |
| **Total Patients** | 57 |
| **Average Samples/Patient** | 6.2 |

---

## 2. SAMPLE DISTRIBUTION BY ACR STATUS (ALL SAMPLES)

| ACR Grade | N Samples | N Patients | % of Total Samples |
|-----------|-----------|------------|--------------------|
| **0R** | 174 | 50 | 49.2% |
| **1R** | 158 | 51 | 44.6% |
| **2R+** | 22 | 18 | 6.2% |
| **TOTAL** | 354 | 57 | 100% |

---

## 3. ACR GROUP CLASSIFICATION

Patients are classified by their **WORST rejection episode**:

| ACR Groups | N Patients | N Samples | Samples/Patient |
|-----------|------------|-----------|-----------------|
| **0R** | 5 | 25 | 5.0 |
| **1R** | 34 | 196 | 5.8 |
| **0R/1R** (combined) | 39 | 221 | 5.7 |
| **2R+** | 18 | 133 | 7.4 |

**Notes:**
- 0R/1R-only group **INCLUDES** the 0R-only patients
- 1R group **EXCLUDES** the 0R-only patients
- 2R+ patients have more samples per patient (surveillance bias - more frequent biopsies)

---

## 4. SAMPLE BREAKDOWN WITHIN 2R+ PATIENTS

Among the **18 patients** who ever develop 2R+ rejection:

| ACR Grade | N Samples | N Patients |
|-----------|-----------|------------|
| 0R | 43 | 15 |
| 1R | 68 | 17 |
| 2R+ | 22 | 18 |
| **TOTAL** | **133** | **18** |

**Key insight:** 2R+ patients contribute samples across ALL ACR grades (0R, 1R, 2R+), reflecting longitudinal sampling before, during, and sometimes after rejection episodes.

---

## 5. POST-OPERATIVE DAY (POD) DISTRIBUTION

### Overall POD Statistics (All Samples)
| Statistic | Days |
|-----------|------|
| Min | 7 |
| Q1 | 21 |
| Median | 42 |
| Mean | 56.3 |
| Q3 | 87 |
| Max | 180 |

### POD Distribution by ACR Status
| ACR Grade | N | Min | Q1 | Median | Q3 | Max |
|-----------|---|-----|-----|--------|-----|-----|
| 0R | 174 | 7 | 28 | 54.5 | 106 | 180 |
| 1R | 158 | 7 | 17.2 | 29 | 57 | 179 |
| 2R+ | 22 | 7 | 10.8 | 20 | 59.5 | 174 |

**Key insight:** 2R+ samples tend to occur earlier (median POD = 20) compared to 0R samples (median POD = 54.5), consistent with early acute rejection.

---

## 6. METABOLITES IN DATASET

**Total metabolites:** 14

### MPA (Mycophenolic Acid) Metabolites (n=3)
1. `MPA (C18)` - Active metabolite measured in serum, C18 chromatography
2. `MPA (HILIC)` - Active metabolite measured in serum, HILIC chromatography
3. `MPAG (C18)` - **Inactive metabolite (glucuronide), renally cleared**

### Corticosteroid Metabolites (n=11)
4. `11-Dehydrocorticosterone (HILIC)`
5. `Prednisone (HILIC)`
6. `Cortisone (HILIC)`
7. `Tetrahydrocorticosterone (C18)`
8. `Tetrahydrocortisone (C18)`
9. `Prednisone (C18)`
10. `Cortisone acetate (HILIC)`
11. `Dexamethasone acetate anhydrous (C18)`
12. `Fludrocortisone acetate (HILIC)`
13. `Fludrocortisone (C18)`
14. `Hydrocortisone sodium succinate (HILIC)`

**Note:** Some metabolites measured using two different chromatography methods (C18 vs HILIC) for technical validation.

---

## 7. CLINICAL VARIABLES AVAILABLE

**Source:** `OHT_Clinical.csv`  
**Clinical data available for:** 61 patients (57 match metabolomics data)

| Variable | Type | Details |
|----------|------|---------|
| **H** | ID | Patient identifier |
| **Age** | Continuous | Patient age at transplant |
| **Race** | Categorical | Asian (n=3), Black (n=33), Hispanic (n=2), White (n=23) |
| **Sex** | Binary | Female (n=21), Male (n=40) |
| **BMI** | Continuous | Body mass index |
| **GFR** | Continuous | Glomerular filtration rate (Min=13, Median=57, Max=140, Missing=1) |

---

## 8. KEY PATIENT GROUPS FOR ANALYSES

### Primary Comparison Groups (defined in `02_setup`)

#### 1. **patients_with_0R_only**
- **Definition:** Patients who NEVER develop rejection
- **N patients:** 5
- **N samples:** 25
- **Use case:** "True controls" - never rejected

#### 2. **patients_with_0R_1R_only**
- **Definition:** Patients who never develop 2R+ (includes 0R-only)
- **N patients:** 39
- **N samples:** 221
- **Use case:** Mild/no rejection group

#### 3. **patients_with_2R**
- **Definition:** Patients who develop at least one 2R+ episode
- **N patients:** 18
- **N samples:** 133 (includes 0R, 1R, and 2R+ samples from these patients)
- **Use case:** Clinically significant rejection group

---

## 9. COMMON ANALYSIS COMPARISONS

### Between-Patient Comparisons
| Comparison | Script | Design |
|------------|--------|--------|
| 0R-only patients vs 2R+ patients | `08b_MPA_POD.R` | All samples from 0R-only vs all samples from 2R+ |
| 0R/1R-only patients vs 2R+ patients | `08a_MPA_POD.R` | All samples from 0R/1R-only vs all samples from 2R+ |

### Within-Patient Comparisons
| Comparison | Script | Design |
|------------|--------|--------|
| 0R samples vs 2R+ samples | `08_MPA_POD.R` | Within 2R+ patients only |

### State-Controlled Comparisons
| Comparison | Script | Design |
|------------|--------|--------|
| 0R from 0R-only vs 0R from 2R+ | `08c_MPA_patient_phenotype.R` | Same ACR state, different patient phenotypes |
| 0R/1R from 0R/1R-only vs 0R/1R from 2R+ | `08c_MPA_patient_phenotype.R` | Same ACR state, different patient phenotypes |

### Stratified Analyses
| Stratification | Script(s) | Variables |
|----------------|-----------|-----------|
| Sex | `11_sex_stratified.R` | Male vs Female |
| Race | `09_race_stratified.R` | Black vs White (others too small) |
| Age | `10_age_stratified.R` | Age quartiles |
| GFR | `12_GFR_stratified.R`, `12a_MPAG_GFR_stratified.R` | 4 CKD stages |
| POD timing | `08_MPA_POD.R` (multiple sections) | Early (≤30d) vs Late (>30d) |

---

## 10. STATISTICAL POWER CONSIDERATIONS

### Sample Size Limitations
- **0R-only patients:** Only 5 patients - **severely underpowered** for between-patient analyses
- **2R+ samples:** Only 22 samples - **limited power** for sample-level comparisons
- **2R+ patients:** 18 patients - **moderate power** for patient-level analyses

### Recommended Approaches
1. **Primary:** Within-patient comparisons (paired designs) maximize power
2. **Secondary:** Patient-level aggregation reduces noise
3. **Exploratory:** Stratified analyses (sex, GFR, race) require caution due to small subgroups

---

## 11. DATA CHARACTERISTICS

### Longitudinal Nature
- **Design:** Repeated measures (median 6.2 samples/patient)
- **Follow-up:** Up to 180 days post-transplant
- **Timing:** Surveillance biopsies at protocol intervals + clinical indication

### Missing Data
- **Metabolomics:** Complete (all metabolites measured in all samples)
- **Clinical:** GFR missing for 1 patient

### Data Format
- **Metabolite values:** Log2-transformed in `IS_metabolites_ACR.csv`
- **For visualization:** Convert back using `2^value`
- **For statistics:** Use log2 values (normalizes distributions, symmetric fold changes)

---

## 12. HOW TO USE THIS DATASET

### Quick Start
```r
# Load everything
source("Scripts/00_source")

# View this summary
source("Scripts/00_data_summary.R")

# Run basic comparisons
source("Scripts/03_tests")  # Sample-level discovery
```

### For MPA-Specific Analyses
```r
# See MPA_ANALYSIS_GUIDE.md for detailed roadmap of:
# - 08_MPA_POD.R (within-patient, 0R vs 2R+)
# - 08a_MPA_POD.R (between-patient, 0R/1R-only vs 2R+)
# - 08b_MPA_POD.R (between-patient, 0R-only vs 2R+)
# - 08c_MPA_patient_phenotype.R (state-controlled)
# - 11_sex_stratified.R (sex stratification)
# - 12_GFR_stratified.R (GFR stratification)
# - 12a_MPAG_GFR_stratified.R (MPAG + GFR)
```

---

## SUMMARY

This is a **moderate-sized, longitudinal metabolomics dataset** with:
- ✅ **Good:** Longitudinal design, clinical metadata, multiple timepoints
- ⚠️ **Limitation:** Small n for 2R+ (n=18 patients, 22 samples)
- ⚠️ **Limitation:** Very small 0R-only group (n=5)
- 💡 **Strength:** Within-patient comparisons (paired designs) maximize power
- 💡 **Strength:** 57 patients with comprehensive clinical annotations

**Recommended analysis strategy:** Focus on within-patient changes and patient-level aggregations. Use stratified analyses cautiously due to sample size.

---

**For questions or updates to this summary, see:**
- `Scripts/00_data_summary.R` (script that generates this)
- `Scripts/02_setup` (data processing pipeline)
- `QUICK_REFERENCE.md` (analysis workflow guide)
- `MPA_ANALYSIS_GUIDE.md` (MPA-specific analyses)
