# TERMINOLOGY UPDATE SUMMARY

**Date:** January 12, 2026  
**Changes:** Updated all references from MMF/Mycophenolate to MPA/MPAG

---

## Rationale for Changes

### Why MPA (not MMF)?
- **MMF (Mycophenolate Mofetil)** = pro-drug, oral medication
- **MPA (Mycophenolic Acid)** = active metabolite, what's actually measured in serum
- **MPAG (Mycophenolic Acid Glucuronide)** = inactive metabolite, renally cleared

**Bottom line:** When measuring drug levels in blood/serum, we're measuring **MPA** (the active form), not MMF (the pro-drug). This is standard clinical pharmacology terminology.

---

## Files Renamed

| Old Filename | New Filename |
|--------------|--------------|
| `08_MMF_POD.R` | `08_MPA_POD.R` |
| `08a_MMF_POD.R` | `08a_MPA_POD.R` |
| `08b_MMF_POD.R` | `08b_MPA_POD.R` |
| `08c_MMF_patient_phenotype.R` | `08c_MPA_patient_phenotype.R` |
| `MMF_ANALYSIS_GUIDE.md` | `MPA_ANALYSIS_GUIDE.md` |

---

## Terminology Changes Throughout All Files

### Text Changes
- `MMF` → `MPA` (in titles, labels, text)
- `Mycophenolate (C18)` → `MPA (C18)`
- `Mycophenolate (HILIC)` → `MPA (HILIC)`
- `Mycophenolic acid O-acyl-glucuronide` → `MPAG`
- `MPA-glucuronide` → `MPAG`

### Variable Names
- `mmf_data` → `mpa_data`
- `mmf_long` → `mpa_long`
- `mmf_0r_only` → `mpa_0r_only`
- `mmf_2r_patients` → `mpa_2r_patients`
- `mmf_0r1r_only` → `mpa_0r1r_only`
- `mmf_metabolites` → `mpa_metabolites`

### Plot Filenames
- `MMF_Overall.png` → `MPA_Overall.png`
- `MMF_Early.png` → `MPA_Early.png`
- `MMF_Late.png` → `MPA_Late.png`
- All other MMF_* plot files → MPA_*

---

## Files Updated

### Analysis Scripts (Primary)
- ✅ `Scripts/08_MPA_POD.R`
- ✅ `Scripts/08a_MPA_POD.R`
- ✅ `Scripts/08b_MPA_POD.R`
- ✅ `Scripts/08c_MPA_patient_phenotype.R`

### Stratified Analysis Scripts
- ✅ `Scripts/09_race_stratified.R`
- ✅ `Scripts/10_age_stratified.R`
- ✅ `Scripts/11_sex_stratified.R`
- ✅ `Scripts/12_GFR_stratified.R`
- ✅ `Scripts/12a_MPAG_GFR_stratified.R`
- ✅ `Scripts/13_GFR_ACR_association.R`

### Utility Scripts
- ✅ `Scripts/00_data_summary.R`

### Documentation
- ✅ `Scripts/MPA_ANALYSIS_GUIDE.md`
- ✅ `DATA_SUMMARY.md`
- ✅ `Rmd/ACR_Immunosuppression.Rmd`

---

## Column Names in Data Files

**NOTE:** The raw data file (`IS_metabolites_ACR.csv`) still has the original column names:
- `Mycophenolate (C18)`
- `Mycophenolate (HILIC)`  
- `Mycophenolic acid O-acyl-glucuronide (C18)`

**This is intentional** - we don't modify raw data files. Instead:
- R scripts reference these columns using backticks: `` `Mycophenolate..C18.` ``
- After reading the data, they're referred to as "MPA" in variable names and plots
- Output labels show "MPA (C18)", "MPA (HILIC)", and "MPAG"

---

## What This Means for Your Analyses

### Correct Interpretation
✅ **"MPA levels"** = serum levels of mycophenolic acid (active metabolite)  
✅ **"MPAG levels"** = serum levels of glucuronide metabolite (inactive, renally cleared)

### Clinical Context
- Patients take **MMF** (oral pro-drug)
- Liver converts MMF → **MPA** (active immunosuppressant)
- Kidneys clear **MPAG** (inactive metabolite)

### Expected Findings
- **MPA:** May differ by rejection status (drug exposure/metabolism)
- **MPAG:** Should correlate with GFR (renal clearance)

---

## Backward Compatibility

### Old References Still Work
The actual column names in the raw data haven't changed, so:
- Old code referencing `` `Mycophenolate..C18.` `` still works
- But new plots/outputs will show "MPA (C18)"

### Plot Files
**Important:** If you've generated plots before this update:
- Old plots are named `MMF_*.png`
- New plots will be named `MPA_*.png`
- You may have duplicates - consider deleting old MMF_* plots

---

## Next Steps

1. **Re-run analyses** to generate plots with updated naming
2. **Update any manuscripts/presentations** to use "MPA" terminology
3. **Delete old MMF_*.png plots** to avoid confusion

---

## Clinical Pharmacology Note

**Proper terminology for publications:**
- "Serum MPA levels" or "MPA concentrations"
- "MPAG levels" or "glucuronide metabolite"
- Can mention "in patients receiving MMF therapy" in methods

**Avoid:**
- ❌ "MMF levels" (MMF is the pro-drug, not measured in serum)
- ❌ "Mycophenolate levels" (ambiguous - is it the pro-drug or active form?)

---

**Questions?** See:
- `MPA_ANALYSIS_GUIDE.md` for analysis roadmap
- `DATA_SUMMARY.md` for dataset overview
- `Scripts/00_data_summary.R` to regenerate data summary
