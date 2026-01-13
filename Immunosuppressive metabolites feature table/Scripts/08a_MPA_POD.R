#! remove

# MPA Analysis with POD Stratification - 0R/1R-only Patients vs 2R Patients
# Comparing ALL samples from patients who have ONLY 0R/1R vs ALL samples from patients who develop 2R
# This is a patient-level phenotype comparison (similar to 08b, but using 0R/1R-only instead of 0R-only)

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# NOTE: patients_with_0R_1R_only and patients_with_2R are already created in 02_setup

# ============================================================================
# VERIFY PATIENT GROUPS
# ============================================================================

cat("\n=== PATIENT GROUPS (from 02_setup) ===\n\n")

cat("Patients with ONLY 0R/1R (non-severe rejectors):", 
    length(unique(patients_with_0R_1R_only$H)), "patients,", 
    nrow(patients_with_0R_1R_only), "samples\n")
cat("Patients who develop 2R at any point:", 
    length(unique(patients_with_2R$H)), "patients,", 
    nrow(patients_with_2R), "samples\n\n")

# ============================================================================
# PART 1: Overall MPA comparison (0R/1R-only patients vs 2R patients - ALL samples)
# ============================================================================

cat("\n=== OVERALL MPA ANALYSIS (0R/1R-only patients vs 2R patients) ===\n")
cat("Comparing: All samples from 0R/1R-only patients vs ALL samples from patients who develop 2R\n")
cat("NOTE: This tests inherent patient-level differences, not just rejection state\n\n")

# Get all samples from 0R/1R-only patients
mpa_0r1r_only <- patients_with_0R_1R_only %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "0R/1R-only patients")

# Get ALL samples from patients who develop 2R (including their 0R, 1R, 2R+ samples)
mpa_2r_patients <- patients_with_2R %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "2R patients (all samples)")

# Combine for comparison
mpa_data <- bind_rows(mpa_0r1r_only, mpa_2r_patients)

cat("Total samples in comparison:\n")
cat("  0R/1R-only patients:", nrow(mpa_0r1r_only), "samples from", 
    length(unique(mpa_0r1r_only$H)), "patients\n")
cat("    - 0R samples:", sum(mpa_0r1r_only$ACR == "0R"), "\n")
cat("    - 1R samples:", sum(mpa_0r1r_only$ACR == "1R"), "\n")
cat("  2R patients:", nrow(mpa_2r_patients), "samples from", 
    length(unique(mpa_2r_patients$H)), "patients (includes 0R, 1R, 2R+ samples)\n")
cat("    - 0R samples:", sum(mpa_2r_patients$ACR == "0R"), "\n")
cat("    - 1R samples:", sum(mpa_2r_patients$ACR == "1R"), "\n")
cat("    - 2R+ samples:", sum(grepl("^2R", mpa_2r_patients$ACR, ignore.case = TRUE)), "\n\n")

# Run Wilcoxon tests
cat("MPA (C18):\n")
wilcox_c18 <- wilcox.test(
    `Mycophenolate..C18.` ~ Patient_Group, 
    data = mpa_data, 
    exact = FALSE
)
print(wilcox_c18)

cat("\nMPA (HILIC):\n")
wilcox_hilic <- wilcox.test(
    `Mycophenolate..HILIC.` ~ Patient_Group, 
    data = mpa_data, 
    exact = FALSE
)
print(wilcox_hilic)

# Plot overall comparison
mpa_long <- mpa_data %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))

# Get sample sizes for x-axis labels
n_0r1r_only <- nrow(mpa_0r1r_only)
n_2r <- nrow(mpa_2r_patients)

p_overall <- ggplot(mpa_long, aes(x = Patient_Group, y = Level, fill = Patient_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    scale_fill_manual(values = c("0R/1R-only patients" = "lightgreen", "2R patients (all samples)" = "lightcoral")) +
    scale_x_discrete(labels = c("0R/1R-only patients" = paste0("0R/1R-only\n(n=", n_0r1r_only, ")"),
                                "2R patients (all samples)" = paste0("2R patients\n(n=", n_2r, ")"))) +
    labs(title = "MPA Levels: 0R/1R-only Patients vs 2R Patients (All Samples)",
         subtitle = paste0("Patient-level comparison | C18: p=", round(wilcox_c18$p.value, 3), 
                          "; HILIC: p=", round(wilcox_hilic$p.value, 3)),
         x = "Patient Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "none")

print(p_overall)

# Create output directory
dir.create("Results2/Feedback_Analysis/0R1R_only_vs_2R", showWarnings = FALSE, recursive = TRUE)
ggsave("Results2/Feedback_Analysis/0R1R_only_vs_2R/MPA_Overall.png", p_overall, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 2: POD-Stratified Analysis (Early ≤30 days vs Late >30 days)
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA ANALYSIS ===\n\n")

# Add POD stratum
mpa_data <- mpa_data %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (≤30 days)", "Late (>30 days)"))

# Function to run tests for one stratum
analyze_stratum <- function(data, stratum) {
    cat("--- ", stratum, " ---\n")
    
    stratum_data <- data %>% filter(POD_Stratum == stratum)
    n_0r1r <- sum(stratum_data$Patient_Group == "0R/1R-only patients")
    n_2r <- sum(stratum_data$Patient_Group == "2R patients (all samples)")
    
    cat("Sample size:", nrow(stratum_data), 
        "(0R/1R-only:", n_0r1r, ", 2R patients:", n_2r, ")\n\n")
    
    # Wilcoxon tests
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ Patient_Group, 
                         data = stratum_data, exact = FALSE)
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ Patient_Group, 
                           data = stratum_data, exact = FALSE)
    
    cat("C18:   p =", round(w_c18$p.value, 3), "\n")
    cat("HILIC: p =", round(w_hilic$p.value, 3), "\n\n")
    
    # Create plot with sample sizes in x-axis labels
    stratum_long <- stratum_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite", values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(stratum_long, aes(x = Patient_Group, y = Level, fill = Patient_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R/1R-only patients" = "lightgreen", "2R patients (all samples)" = "lightcoral")) +
        scale_x_discrete(labels = c("0R/1R-only patients" = paste0("0R/1R-only\n(n=", n_0r1r, ")"),
                                    "2R patients (all samples)" = paste0("2R patients\n(n=", n_2r, ")"))) +
        labs(title = paste("MPA Levels:", stratum, "(0R/1R-only vs 2R patients)"),
             subtitle = paste0("C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "Patient Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(plot = p, p_c18 = w_c18$p.value, p_hilic = w_hilic$p.value))
}

# Run for each stratum
early_results <- analyze_stratum(mpa_data, "Early (≤30 days)")
print(early_results$plot)
ggsave("Results2/Feedback_Analysis/0R1R_only_vs_2R/MPA_Early.png", early_results$plot, 
       width = 10, height = 6, dpi = 300)

late_results <- analyze_stratum(mpa_data, "Late (>30 days)")
print(late_results$plot)
ggsave("Results2/Feedback_Analysis/0R1R_only_vs_2R/MPA_Late.png", late_results$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 3: POD-Stratified Analysis EXCLUDING POD<10
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA ANALYSIS (EXCLUDING POD<10) ===\n")
cat("NOTE: Removing very early post-transplant samples (POD<10)\n\n")

# Filter out POD<10 samples
mmf_0r1r_pod10plus <- patients_with_0R_1R_only %>%
    filter(POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "0R/1R-only patients")

mmf_2r_pod10plus <- patients_with_2R %>%
    filter(POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "2R patients (all samples)")

mmf_data_pod10plus <- bind_rows(mmf_0r1r_pod10plus, mmf_2r_pod10plus)

cat("Total samples with POD >= 10:", nrow(mmf_data_pod10plus), "\n")
cat("0R/1R-only patients:", nrow(mmf_0r1r_pod10plus), "\n")
cat("2R patients (all samples):", nrow(mmf_2r_pod10plus), "\n\n")

# Add POD stratum
mmf_data_pod10plus <- mmf_data_pod10plus %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (>30 days)"))

# Run for each stratum
early_results_pod10 <- analyze_stratum(mmf_data_pod10plus, "Early (10-30 days)")
print(early_results_pod10$plot)
ggsave("Results2/Feedback_Analysis/0R1R_only_vs_2R/MPA_Early_POD10plus.png", early_results_pod10$plot, 
       width = 10, height = 6, dpi = 300)

late_results_pod10 <- analyze_stratum(mmf_data_pod10plus, "Late (>30 days)")
print(late_results_pod10$plot)
ggsave("Results2/Feedback_Analysis/0R1R_only_vs_2R/MPA_Late_POD10plus.png", late_results_pod10$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 4: POD-Stratified Analysis EXCLUDING POD<10 AND POD>45
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA ANALYSIS (POD 10-45 ONLY) ===\n")
cat("NOTE: Removing very early (POD<10) and late (POD>45) samples\n")
cat("      to focus on the acute post-transplant window\n\n")

# Filter to POD 10-45 only
mmf_0r1r_pod10_45 <- patients_with_0R_1R_only %>%
    filter(POD >= 10 & POD <= 45) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "0R/1R-only patients")

mmf_2r_pod10_45 <- patients_with_2R %>%
    filter(POD >= 10 & POD <= 45) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "2R patients (all samples)")

mmf_data_pod10_45 <- bind_rows(mmf_0r1r_pod10_45, mmf_2r_pod10_45)

cat("Total samples with POD 10-45:", nrow(mmf_data_pod10_45), "\n")
cat("0R/1R-only patients:", nrow(mmf_0r1r_pod10_45), "\n")
cat("2R patients (all samples):", nrow(mmf_2r_pod10_45), "\n\n")

# Add POD stratum
mmf_data_pod10_45 <- mmf_data_pod10_45 %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (31-45 days)"))

# Run for each stratum
early_results_pod10_45 <- analyze_stratum(mmf_data_pod10_45, "Early (10-30 days)")
print(early_results_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/0R1R_only_vs_2R/MPA_Early_POD10_45.png", early_results_pod10_45$plot, 
       width = 10, height = 6, dpi = 300)

late_results_pod10_45 <- analyze_stratum(mmf_data_pod10_45, "Late (31-45 days)")
print(late_results_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/0R1R_only_vs_2R/MPA_Late_POD10_45.png", late_results_pod10_45$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== SUMMARY (0R/1R-ONLY vs 2R PATIENTS COMPARISON) ===\n")
cat("Overall:      C18 p=", round(wilcox_c18$p.value, 3), 
    ", HILIC p=", round(wilcox_hilic$p.value, 3), "\n")
cat("Early ≤30d:   C18 p=", round(early_results$p_c18, 3), 
    ", HILIC p=", round(early_results$p_hilic, 3), "\n")
cat("Late >30d:    C18 p=", round(late_results$p_c18, 3), 
    ", HILIC p=", round(late_results$p_hilic, 3), "\n")
cat("Early 10-30d: C18 p=", round(early_results_pod10$p_c18, 3), 
    ", HILIC p=", round(early_results_pod10$p_hilic, 3), "\n")
cat("Late >30d:    C18 p=", round(late_results_pod10$p_c18, 3), 
    ", HILIC p=", round(late_results_pod10$p_hilic, 3), "\n")
cat("Early 10-30d: C18 p=", round(early_results_pod10_45$p_c18, 3), 
    ", HILIC p=", round(early_results_pod10_45$p_hilic, 3), " (POD 10-45 only)\n")
cat("Late 31-45d:  C18 p=", round(late_results_pod10_45$p_c18, 3), 
    ", HILIC p=", round(late_results_pod10_45$p_hilic, 3), " (POD 10-45 only)\n")
cat("\nPlots saved to: Results2/Feedback_Analysis/0R1R_only_vs_2R/\n\n")

cat("\n=== INTERPRETATION ===\n")
cat("This analysis compares PATIENT-LEVEL differences:\n")
cat("  - 0R/1R-only patients: Never develop severe (2R+) rejection (n=", 
    length(unique(mpa_0r1r_only$H)), " patients, all samples)\n", sep = "")
cat("  - 2R patients: Develop 2R at some point (n=", 
    length(unique(mpa_2r_patients$H)), " patients, ALL samples including 0R/1R/2R+)\n\n", sep = "")

cat("Key Insights:\n")
cat("  - Tests inherent patient phenotype differences (similar to 08b)\n")
cat("  - Broader patient group than 08b (includes 1R patients)\n")
cat("  - More statistical power than 0R-only comparisons\n")
cat("  - Mixes rejection states within the 2R patient group\n\n")

cat("Comparison with other analyses:\n")
cat("  08_MPA_POD.R:  Within-patient (0R vs 2R+ samples in patients who reject)\n")
cat("                 → Tests rejection STATE effects\n\n")
cat("  08a_MPA_POD.R: 0R/1R-only patients vs 2R patients (all samples) [THIS SCRIPT]\n")
cat("                 → Patient phenotype with broader inclusion (0R/1R-only)\n\n")
cat("  08b_MPA_POD.R: 0R-only patients vs 2R patients (all samples)\n")
cat("                 → Patient phenotype with strictest non-rejectors (0R-only)\n\n")
cat("  08c_MMF_POD.R: Patient phenotype controlled for rejection state\n")
cat("                 → 0R vs 0R; 0R/1R vs 0R/1R (cleanest patient comparison)\n\n")

cat("Clinical Implications:\n")
cat("  - If significant: Patients prone to severe rejection differ from non-severe rejectors\n")
cat("  - Compare with 08b: Does including 1R patients change the conclusions?\n")
cat("  - Compare with 08c: Are differences driven by rejection state or patient phenotype?\n\n")





