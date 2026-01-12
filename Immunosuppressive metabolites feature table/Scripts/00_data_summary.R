# =============================================================================
# COMPREHENSIVE DATA SUMMARY
# =============================================================================
# This script generates a complete overview of the dataset including:
# - Total samples and patients
# - Sample distribution by ACR status
# - Patient phenotype classification
# - Metabolite list
# - Clinical variables available
#
# Run this after sourcing 00_source
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("                    COMPREHENSIVE DATA SUMMARY                                 \n")
cat("================================================================================\n\n")

# =============================================================================
# SECTION 1: OVERALL DATASET OVERVIEW
# =============================================================================

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("1. OVERALL DATASET OVERVIEW\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Load data
data_path <- "/Users/ailintang/Desktop/Immunosuppressive_ACR/Immunosuppressive metabolites feature table"
metabolite_data <- as_tibble(read.csv(file.path(data_path, "IS_metabolites_ACR.csv")))
nonlog2_data <- metabolite_data %>% 
    mutate(across(6:ncol(.), ~ 2^.x))

# Basic counts
n_total_samples <- nrow(nonlog2_data)
n_total_patients <- n_distinct(nonlog2_data$H)

cat("Total Samples:  ", n_total_samples, "\n")
cat("Total Patients: ", n_total_patients, "\n\n")

# =============================================================================
# SECTION 2: SAMPLE DISTRIBUTION BY ACR STATUS (ALL SAMPLES)
# =============================================================================

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("2. SAMPLE DISTRIBUTION BY ACR STATUS (ALL SAMPLES)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

acr_summary <- nonlog2_data %>%
    group_by(ACR) %>%
    summarise(
        n_samples = n(),
        n_patients = n_distinct(H),
        .groups = "drop"
    ) %>%
    arrange(ACR)

print(acr_summary)

# Calculate 2R+ total
n_2r_samples <- sum(grepl("^2R", nonlog2_data$ACR, ignore.case = TRUE))
n_2r_patients <- nonlog2_data %>%
    filter(grepl("^2R", ACR, ignore.case = TRUE)) %>%
    summarise(n = n_distinct(H)) %>%
    pull(n)

cat("\nCombined 2R+ samples: ", n_2r_samples, " (from ", n_2r_patients, " patients)\n\n")

# =============================================================================
# SECTION 3: PATIENT PHENOTYPE CLASSIFICATION
# =============================================================================

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("3. PATIENT PHENOTYPE CLASSIFICATION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("Patients are classified by their WORST rejection episode:\n\n")

# 0R-only patients
patients_with_0R_only <- nonlog2_data %>%
    group_by(H) %>%
    filter(all(ACR == "0R")) %>%
    ungroup()

n_patients_0R_only <- n_distinct(patients_with_0R_only$H)
n_samples_0R_only <- nrow(patients_with_0R_only)

# 0R/1R-only patients (includes patients with ONLY 0R)
patients_with_0R_1R_only <- nonlog2_data %>%
    group_by(H) %>%
    filter(all(ACR %in% c("0R", "1R"))) %>%
    ungroup()

n_patients_0R_1R_only <- n_distinct(patients_with_0R_1R_only$H)
n_samples_0R_1R_only <- nrow(patients_with_0R_1R_only)

# 1R-only patients (max rejection is 1R, excludes patients with only 0R)
patients_with_1R_max <- nonlog2_data %>%
    group_by(H) %>%
    filter(any(ACR == "1R") & all(ACR %in% c("0R", "1R"))) %>%
    ungroup()

n_patients_1R_max <- n_distinct(patients_with_1R_max$H)
n_samples_1R_max <- nrow(patients_with_1R_max)

# 2R+ patients
patients_with_2R <- nonlog2_data %>%
    group_by(H) %>%
    filter(any(grepl("^2R", ACR, ignore.case = TRUE))) %>%
    ungroup()

n_patients_2R <- n_distinct(patients_with_2R$H)
n_samples_2R <- nrow(patients_with_2R)

# Create summary table
phenotype_summary <- tibble(
    Phenotype = c("0R-only", "1R-max (excludes 0R-only)", "0R/1R-only (combined)", "2R+ (ever)"),
    N_Patients = c(n_patients_0R_only, n_patients_1R_max, n_patients_0R_1R_only, n_patients_2R),
    N_Samples = c(n_samples_0R_only, n_samples_1R_max, n_samples_0R_1R_only, n_samples_2R),
    Samples_per_Patient = round(c(n_samples_0R_only, n_samples_1R_max, n_samples_0R_1R_only, n_samples_2R) / 
                                 c(n_patients_0R_only, n_patients_1R_max, n_patients_0R_1R_only, n_patients_2R), 1)
)

print(phenotype_summary)
cat("\nNOTE: 0R/1R-only group INCLUDES the 0R-only patients\n")
cat("      1R-max group EXCLUDES the 0R-only patients\n\n")

# =============================================================================
# SECTION 4: WITHIN 2R+ PATIENTS - SAMPLE BREAKDOWN BY ACR STATUS
# =============================================================================

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("4. SAMPLE BREAKDOWN WITHIN 2R+ PATIENTS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("Among the ", n_patients_2R, " patients who ever develop 2R+ rejection:\n\n")

within_2r_summary <- patients_with_2R %>%
    group_by(ACR) %>%
    summarise(
        n_samples = n(),
        n_patients = n_distinct(H),
        .groups = "drop"
    ) %>%
    arrange(ACR)

print(within_2r_summary)

cat("\nThese patients have ", n_samples_2R, " total samples across all ACR grades\n\n")

# =============================================================================
# SECTION 5: POD (POST-OPERATIVE DAY) DISTRIBUTION
# =============================================================================

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("5. POST-OPERATIVE DAY (POD) DISTRIBUTION\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("Overall POD statistics:\n")
print(summary(nonlog2_data$POD))

cat("\nPOD distribution by ACR status:\n")
pod_by_acr <- nonlog2_data %>%
    group_by(ACR) %>%
    summarise(
        n = n(),
        Min = min(POD, na.rm = TRUE),
        Q1 = quantile(POD, 0.25, na.rm = TRUE),
        Median = median(POD, na.rm = TRUE),
        Q3 = quantile(POD, 0.75, na.rm = TRUE),
        Max = max(POD, na.rm = TRUE),
        .groups = "drop"
    )
print(pod_by_acr)
cat("\n")

# =============================================================================
# SECTION 6: METABOLITES IN DATASET
# =============================================================================

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("6. METABOLITES IN DATASET\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

metabolite_names <- names(metabolite_data)[6:ncol(metabolite_data)]
n_metabolites <- length(metabolite_names)

cat("Total metabolites: ", n_metabolites, "\n\n")

cat("Metabolite list:\n")
for (i in 1:length(metabolite_names)) {
    cat(sprintf("%2d. %s\n", i, metabolite_names[i]))
}

# Group metabolites by drug class
cat("\n\nMetabolites by drug class:\n\n")

mpa_metabolites <- metabolite_names[grepl("Mycophenol", metabolite_names, ignore.case = TRUE)]
cat("MPA (Mycophenolic Acid) metabolites (n=", length(mpa_metabolites), "):\n", sep = "")
for (m in mpa_metabolites) cat("  - ", m, "\n", sep = "")

cat("\n")

steroid_metabolites <- metabolite_names[!grepl("Mycophenol", metabolite_names, ignore.case = TRUE)]
cat("Corticosteroid metabolites (n=", length(steroid_metabolites), "):\n", sep = "")
for (m in steroid_metabolites) cat("  - ", m, "\n", sep = "")

cat("\n")

# =============================================================================
# SECTION 7: CLINICAL VARIABLES AVAILABLE
# =============================================================================

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("7. CLINICAL VARIABLES AVAILABLE\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Check if clinical data exists
clinical_file <- file.path(data_path, "..", "OHT_Clinical.csv")
if (file.exists(clinical_file)) {
    clinical_data <- read_csv(clinical_file, show_col_types = FALSE)
    
    cat("Clinical variables (from OHT_Clinical.csv):\n\n")
    
    clinical_vars <- names(clinical_data)
    for (i in 1:length(clinical_vars)) {
        var <- clinical_vars[i]
        cat(sprintf("%2d. %s\n", i, var))
        
        # Show summary for key variables
        if (var == "GFR") {
            cat("    Summary: ")
            cat("Min=", round(min(clinical_data$GFR, na.rm = TRUE), 1),
                ", Median=", round(median(clinical_data$GFR, na.rm = TRUE), 1),
                ", Max=", round(max(clinical_data$GFR, na.rm = TRUE), 1),
                ", Missing=", sum(is.na(clinical_data$GFR)), "\n", sep = "")
        } else if (var %in% c("Sex", "Race", "Age_Group")) {
            tab <- table(clinical_data[[var]], useNA = "ifany")
            cat("    Distribution: ", paste(names(tab), "=", tab, collapse = ", "), "\n")
        }
    }
    
    cat("\nClinical data available for ", nrow(clinical_data), " patients\n", sep = "")
    cat("Match with metabolomics: ", sum(clinical_data$H %in% nonlog2_data$H), " patients\n\n", sep = "")
} else {
    cat("Clinical data file not found at expected location\n\n")
}

# =============================================================================
# SECTION 8: KEY PATIENT GROUPS FOR ANALYSES
# =============================================================================

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("8. KEY PATIENT GROUPS FOR ANALYSES\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("Primary comparison groups defined in 02_setup:\n\n")

cat("1. patients_with_0R_only\n")
cat("   - Patients who NEVER develop rejection\n")
cat("   - N patients: ", n_patients_0R_only, "\n", sep = "")
cat("   - N samples:  ", n_samples_0R_only, "\n\n", sep = "")

cat("2. patients_with_0R_1R_only\n")
cat("   - Patients who never develop 2R+ (includes 0R-only)\n")
cat("   - N patients: ", n_patients_0R_1R_only, "\n", sep = "")
cat("   - N samples:  ", n_samples_0R_1R_only, "\n\n", sep = "")

cat("3. patients_with_2R\n")
cat("   - Patients who develop at least one 2R+ episode\n")
cat("   - N patients: ", n_patients_2R, "\n", sep = "")
cat("   - N samples:  ", n_samples_2R, " (includes 0R, 1R, and 2R+ samples)\n\n", sep = "")

cat("Common comparisons in analyses:\n")
cat("  - 0R-only patients vs 2R+ patients (08b_MPA_POD.R)\n")
cat("  - 0R/1R-only patients vs 2R+ patients (08a_MPA_POD.R)\n")
cat("  - Within 2R+ patients: 0R samples vs 2R+ samples (08_MPA_POD.R)\n")
cat("  - State-controlled: 0R from 0R-only vs 0R from 2R+ (08c_MPA_patient_phenotype.R)\n\n")

# =============================================================================
# SUMMARY FOOTER
# =============================================================================

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("                           END OF SUMMARY                                     \n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("To run this summary:\n")
cat("  source('Immunosuppressive metabolites feature table/Scripts/00_source')\n")
cat("  source('Immunosuppressive metabolites feature table/Scripts/00_data_summary.R')\n\n")
