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

### What Changed ✅
- **File names:** `*_MMF_*` → `*_MPA_*`
- **Titles/Labels:** `"MMF"` → `"MPA"` (in plot titles, output text)
- **Variable names:** `mmf_data` → `mpa_data`, `mmf_long` → `mpa_long`, etc.
- **Plot filenames:** `MMF_Overall.png` → `MPA_Overall.png`, etc.
- **Output text:** References to "MMF" in printed output → "MPA"

### What Did NOT Change ❌
- **Column names in data:** Still `Mycophenolate..C18.`, `Mycophenolate..HILIC.`, `Mycophenolic.acid.O.acyl.glucuronide..C18.`
- **Column references in code:** Still use backticks with original names: `` `Mycophenolate..C18.` ``
- **The actual data file:** `IS_metabolites_ACR.csv` unchanged

### Summary
- **In the code:** Use original column names (`` `Mycophenolate..C18.` ``)
- **In the output:** Display as "MPA (C18)", "MPA (HILIC)", "MPAG"
- **In filenames/variables:** Use `mpa_*` naming

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

**IMPORTANT:** The raw data file (`IS_metabolites_ACR.csv`) column names are **UNCHANGED**:
- `Mycophenolate (C18)` → becomes `Mycophenolate..C18.` in R
- `Mycophenolate (HILIC)` → becomes `Mycophenolate..HILIC.` in R
- `Mycophenolic acid O-acyl-glucuronide (C18)` → becomes `Mycophenolic.acid.O.acyl.glucuronide..C18.` in R

**In the R code:**
- ✅ Use original column names: `` `Mycophenolate..C18.` ``
- ✅ Use original column names: `` `Mycophenolic.acid.O.acyl.glucuronide..C18.` ``
- ❌ Do NOT use: `` `MPA..C18.` `` (this will cause errors!)

**In the output/labels:**
- When displaying to user, show as "MPA (C18)", "MPA (HILIC)", "MPAG"
- In plot titles, axis labels, printed text: use "MPA" terminology
- The data transformation happens via `gsub("\\.\\.", " ", Metabolite)` after pivoting

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
