# MPA Analysis: Patient Phenotype Comparison (Rejection State Controlled)
# Comparing 0R samples from 0R-only patients vs 0R samples from 2R patients
# This isolates PATIENT-LEVEL differences while controlling for rejection state

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# NOTE: patients_with_0R_only and patients_with_2R are already created in 02_setup

# ============================================================================
# VERIFY PATIENT GROUPS
# ============================================================================

cat("\n=== PATIENT PHENOTYPE COMPARISON (0R samples only) ===\n\n")
cat("Research Question: Do 0R samples differ between patients who never reject\n")
cat("                   vs patients who will develop 2R?\n\n")

cat("Patients with ONLY 0R (never reject):", 
    length(unique(patients_with_0R_only$H)), "patients,", 
    nrow(patients_with_0R_only), "samples (all 0R)\n")
cat("Patients who develop 2R at any point:", 
    length(unique(patients_with_2R$H)), "patients,", 
    nrow(patients_with_2R), "total samples\n")
cat("  - 0R samples from 2R patients:", 
    sum(patients_with_2R$ACR == "0R"), "samples\n\n")

# ============================================================================
#- PART 1: Overall 0R Sample Comparison
# ============================================================================

cat("\n=== OVERALL 0R SAMPLE COMPARISON ===\n")
cat("Comparing: 0R samples from 0R-only patients vs 0R samples from 2R patients\n")
cat("NOTE: Same rejection state (0R), different patient phenotypes\n\n")

# Get 0R samples from 0R-only patients (all their samples are 0R)
mpa_0r_only <- patients_with_0R_only %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "0R-only patients")

# Get ONLY 0R samples from patients who develop 2R
mpa_0r_from_2r <- patients_with_2R %>%
    filter(ACR == "0R") %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "2R patients (0R samples)")

# Combine for comparison
mpa_data <- bind_rows(mpa_0r_only, mpa_0r_from_2r)

cat("Sample sizes:\n")
cat("  0R samples from 0R-only patients:", nrow(mpa_0r_only), "samples from", 
    length(unique(mpa_0r_only$H)), "patients\n")
cat("  0R samples from 2R patients:", nrow(mpa_0r_from_2r), "samples from", 
    length(unique(mpa_0r_from_2r$H)), "patients\n\n")

# Run Wilcoxon tests
cat("MPA (C18):\n")
wilcox_c18 <- wilcox.test(
    `Mycophenolate..C18.` ~ Patient_Phenotype, 
    data = mpa_data, 
    exact = FALSE
)
print(wilcox_c18)

cat("\nMPA (HILIC):\n")
wilcox_hilic <- wilcox.test(
    `Mycophenolate..HILIC.` ~ Patient_Phenotype, 
    data = mpa_data, 
    exact = FALSE
)
print(wilcox_hilic)

# Get medians for interpretation
medians_c18 <- mpa_data %>%
    group_by(Patient_Phenotype) %>%
    summarise(
        median_C18 = median(`Mycophenolate..C18.`, na.rm = TRUE),
        median_HILIC = median(`Mycophenolate..HILIC.`, na.rm = TRUE)
    )
print(medians_c18)

# Plot overall comparison
mpa_long <- mpa_data %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))

# Get sample sizes for x-axis labels
n_0r_only <- nrow(mpa_0r_only)
n_0r_from_2r <- nrow(mpa_0r_from_2r)

p_overall <- ggplot(mpa_long, aes(x = Patient_Phenotype, y = Level, fill = Patient_Phenotype)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    scale_fill_manual(values = c("0R-only patients" = "lightblue", 
                                  "2R patients (0R samples)" = "lightcoral")) +
    scale_x_discrete(labels = c("0R-only patients" = paste0("0R-only\npatients\n(n=", n_0r_only, ")"),
                                "2R patients (0R samples)" = paste0("2R patients\n(0R samples)\n(n=", n_0r_from_2r, ")"))) +
    labs(title = "MPA Levels in 0R Samples: Patient Phenotype Comparison",
         subtitle = paste0("0R-only patients vs 0R samples from 2R patients | C18: p=", 
                          round(wilcox_c18$p.value, 3), "; HILIC: p=", round(wilcox_hilic$p.value, 3)),
         x = "Patient Phenotype", y = "Level") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 9))

print(p_overall)

# Create output directory
dir.create("Results2/Feedback_Analysis/Patient_Phenotype", showWarnings = FALSE, recursive = TRUE)
ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R_Samples_Overall.png", p_overall, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 2: POD-Stratified Analysis (Early ≤30 days vs Late >30 days)
# ============================================================================

cat("\n\n=== POD-STRATIFIED 0R SAMPLE COMPARISON ===\n\n")

# Add POD stratum
mpa_data <- mpa_data %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (≤30 days)", "Late (>30 days)"))

# Function to run tests for one stratum
analyze_stratum <- function(data, stratum) {
    cat("--- ", stratum, " ---\n")
    
    stratum_data <- data %>% filter(POD_Stratum == stratum)
    n_0r_only <- sum(stratum_data$Patient_Phenotype == "0R-only patients")
    n_0r_from_2r <- sum(stratum_data$Patient_Phenotype == "2R patients (0R samples)")
    
    cat("Sample size:", nrow(stratum_data), 
        "(0R-only:", n_0r_only, ", 0R from 2R patients:", n_0r_from_2r, ")\n")
    
    # Check if we have enough samples
    if (n_0r_only < 2 || n_0r_from_2r < 2) {
        cat("Insufficient samples for comparison\n\n")
        return(NULL)
    }
    
    # Wilcoxon tests
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ Patient_Phenotype, 
                         data = stratum_data, exact = FALSE)
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ Patient_Phenotype, 
                           data = stratum_data, exact = FALSE)
    
    cat("C18:   p =", round(w_c18$p.value, 3), "\n")
    cat("HILIC: p =", round(w_hilic$p.value, 3), "\n\n")
    
    # Get medians
    meds <- stratum_data %>%
        group_by(Patient_Phenotype) %>%
        summarise(
            median_C18 = median(`Mycophenolate..C18.`, na.rm = TRUE),
            median_HILIC = median(`Mycophenolate..HILIC.`, na.rm = TRUE)
        )
    print(meds)
    cat("\n")
    
    # Create plot with sample sizes in x-axis labels
    stratum_long <- stratum_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite", values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(stratum_long, aes(x = Patient_Phenotype, y = Level, fill = Patient_Phenotype)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R-only patients" = "lightblue", 
                                      "2R patients (0R samples)" = "lightcoral")) +
        scale_x_discrete(labels = c("0R-only patients" = paste0("0R-only\npatients\n(n=", n_0r_only, ")"),
                                    "2R patients (0R samples)" = paste0("2R patients\n(0R samples)\n(n=", n_0r_from_2r, ")"))) +
        labs(title = paste("MPA Levels in 0R Samples:", stratum),
             subtitle = paste0("Baseline state comparison | C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "Patient Phenotype", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 9))
    
    return(list(plot = p, p_c18 = w_c18$p.value, p_hilic = w_hilic$p.value))
}

# Run for each stratum
early_results <- analyze_stratum(mpa_data, "Early (≤30 days)")
if (!is.null(early_results)) {
    print(early_results$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R_Samples_Early.png", 
           early_results$plot, width = 10, height = 6, dpi = 300)
}

late_results <- analyze_stratum(mpa_data, "Late (>30 days)")
if (!is.null(late_results)) {
    print(late_results$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R_Samples_Late.png", 
           late_results$plot, width = 10, height = 6, dpi = 300)
}

# ============================================================================
# PART 3: POD-Stratified Analysis EXCLUDING POD<10
# ============================================================================

cat("\n\n=== POD-STRATIFIED 0R SAMPLE COMPARISON (EXCLUDING POD<10) ===\n")
cat("NOTE: Removing very early post-transplant samples (POD<10)\n\n")

# Filter out POD<10 samples
mmf_0r_only_pod10plus <- patients_with_0R_only %>%
    filter(POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "0R-only patients")

mmf_0r_from_2r_pod10plus <- patients_with_2R %>%
    filter(ACR == "0R", POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "2R patients (0R samples)")

mmf_data_pod10plus <- bind_rows(mmf_0r_only_pod10plus, mmf_0r_from_2r_pod10plus)

cat("Sample sizes with POD >= 10:\n")
cat("  0R samples from 0R-only patients:", nrow(mmf_0r_only_pod10plus), "\n")
cat("  0R samples from 2R patients:", nrow(mmf_0r_from_2r_pod10plus), "\n\n")

# Add POD stratum
mmf_data_pod10plus <- mmf_data_pod10plus %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (>30 days)"))

# Run for each stratum
early_results_pod10 <- analyze_stratum(mmf_data_pod10plus, "Early (10-30 days)")
if (!is.null(early_results_pod10)) {
    print(early_results_pod10$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R_Samples_Early_POD10plus.png", 
           early_results_pod10$plot, width = 10, height = 6, dpi = 300)
}

late_results_pod10 <- analyze_stratum(mmf_data_pod10plus, "Late (>30 days)")
if (!is.null(late_results_pod10)) {
    print(late_results_pod10$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R_Samples_Late_POD10plus.png", 
           late_results_pod10$plot, width = 10, height = 6, dpi = 300)
}

# ============================================================================
# PART 4: POD-Stratified Analysis EXCLUDING POD<10 AND POD>45
# ============================================================================

cat("\n\n=== POD-STRATIFIED 0R SAMPLE COMPARISON (POD 10-45 ONLY) ===\n")
cat("NOTE: Focusing on acute post-transplant window (POD 10-45)\n\n")

# Filter to POD 10-45 only
mmf_0r_only_pod10_45 <- patients_with_0R_only %>%
    filter(POD >= 10 & POD <= 45) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "0R-only patients")

mmf_0r_from_2r_pod10_45 <- patients_with_2R %>%
    filter(ACR == "0R", POD >= 10, POD <= 45) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "2R patients (0R samples)")

mmf_data_pod10_45 <- bind_rows(mmf_0r_only_pod10_45, mmf_0r_from_2r_pod10_45)

cat("Sample sizes with POD 10-45:\n")
cat("  0R samples from 0R-only patients:", nrow(mmf_0r_only_pod10_45), "\n")
cat("  0R samples from 2R patients:", nrow(mmf_0r_from_2r_pod10_45), "\n\n")

# Add POD stratum
mmf_data_pod10_45 <- mmf_data_pod10_45 %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (31-45 days)"))

# Run for each stratum
early_results_pod10_45 <- analyze_stratum(mmf_data_pod10_45, "Early (10-30 days)")
if (!is.null(early_results_pod10_45)) {
    print(early_results_pod10_45$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R_Samples_Early_POD10_45.png", 
           early_results_pod10_45$plot, width = 10, height = 6, dpi = 300)
}

late_results_pod10_45 <- analyze_stratum(mmf_data_pod10_45, "Late (31-45 days)")
if (!is.null(late_results_pod10_45)) {
    print(late_results_pod10_45$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R_Samples_Late_POD10_45.png", 
           late_results_pod10_45$plot, width = 10, height = 6, dpi = 300)
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== SUMMARY (0R SAMPLES: PATIENT PHENOTYPE COMPARISON) ===\n")
cat("Overall:      C18 p=", round(wilcox_c18$p.value, 3), 
    ", HILIC p=", round(wilcox_hilic$p.value, 3), "\n")

if (!is.null(early_results)) {
    cat("Early ≤30d:   C18 p=", round(early_results$p_c18, 3), 
        ", HILIC p=", round(early_results$p_hilic, 3), "\n")
}
if (!is.null(late_results)) {
    cat("Late >30d:    C18 p=", round(late_results$p_c18, 3), 
        ", HILIC p=", round(late_results$p_hilic, 3), "\n")
}
if (!is.null(early_results_pod10)) {
    cat("Early 10-30d: C18 p=", round(early_results_pod10$p_c18, 3), 
        ", HILIC p=", round(early_results_pod10$p_hilic, 3), "\n")
}
if (!is.null(late_results_pod10)) {
    cat("Late >30d:    C18 p=", round(late_results_pod10$p_c18, 3), 
        ", HILIC p=", round(late_results_pod10$p_hilic, 3), "\n")
}
if (!is.null(early_results_pod10_45)) {
    cat("Early 10-30d: C18 p=", round(early_results_pod10_45$p_c18, 3), 
        ", HILIC p=", round(early_results_pod10_45$p_hilic, 3), " (POD 10-45 only)\n")
}
if (!is.null(late_results_pod10_45)) {
    cat("Late 31-45d:  C18 p=", round(late_results_pod10_45$p_c18, 3), 
        ", HILIC p=", round(late_results_pod10_45$p_hilic, 3), " (POD 10-45 only)\n")
}

cat("\nPlots saved to: Results2/Feedback_Analysis/Patient_Phenotype/\n\n")

# ============================================================================
# INTERPRETATION
# ============================================================================

cat("\n=== CLINICAL INTERPRETATION ===\n\n")
cat("This is the CLEANEST patient-level comparison:\n")
cat("  - Same rejection state (0R) in both groups\n")
cat("  - Controls for acute inflammatory effects of rejection\n")
cat("  - Tests inherent patient-level differences in MPA during quiescence\n\n")

cat("Patient Groups:\n")
cat("  - 0R-only patients: Never have ANY rejection (n=", 
    length(unique(mpa_0r_only$H)), " patients)\n", sep = "")
cat("  - 2R patients (0R samples): Samples taken when NOT rejecting from patients who DO develop 2R (n=", 
    length(unique(mpa_0r_from_2r$H)), " patients)\n\n", sep = "")

cat("Key Questions Addressed:\n")
cat("  1. Do rejection-prone patients have lower MPA levels even when NOT rejecting?\n")
cat("  2. Can baseline MPA levels predict future rejection risk?\n")
cat("  3. Are patient-level differences independent of rejection state?\n\n")

cat("Comparison with other analyses:\n")
cat("  08_MPA_POD.R:  Within-patient comparison (0R vs 2R+ in same patients)\n")
cat("                 → Tests rejection STATE effects\n\n")
cat("  08a_MPA_POD.R: 0R/1R-only patients vs 2R+ samples from 2R patients\n")
cat("                 → Confounds patient phenotype AND rejection state\n\n")
cat("  08b_MPA_POD.R: 0R-only patients vs ALL samples from 2R patients\n")
cat("                 → Patient phenotype with mixed rejection states\n\n")
cat("  08c_MMF_POD.R: 0R from 0R-only patients vs 0R from 2R patients (THIS SCRIPT)\n")
cat("                 → Pure patient phenotype, rejection state controlled ✓\n\n")

cat("Clinical Implications:\n")
cat("  - If significant: MPA levels during quiescence differ by rejection phenotype\n")
cat("    → Suggests inherent patient-level pharmacokinetic/pharmacodynamic differences\n")
cat("    → Potential for early risk stratification\n")
cat("  - If not significant: Rejection appears to be STATE-dependent, not phenotype\n")
cat("    → Focus on monitoring during at-risk periods\n\n")

# ============================================================================
#- PART 5: ADDITIONAL COMPARISON - 0R/1R samples
# ============================================================================

cat("\n\n" , rep("=", 76), "\n", sep = "")
cat("=== ADDITIONAL ANALYSIS: 0R/1R SAMPLE COMPARISON ===\n")
cat(rep("=", 76), "\n\n", sep = "")

cat("Research Question: Do 0R/1R samples differ between patients who only have 0R/1R\n")
cat("                   vs patients who will develop 2R?\n\n")

cat("Patient Groups:\n")
cat("  - 0R/1R-only patients: Never develop 2R (n=", 
    length(unique(patients_with_0R_1R_only$H)), " patients, ",
    nrow(patients_with_0R_1R_only), " samples)\n", sep = "")
cat("  - 2R patients: Develop 2R at some point (n=", 
    length(unique(patients_with_2R$H)), " patients)\n", sep = "")
cat("    - 0R/1R samples from 2R patients: ", 
    sum(patients_with_2R$ACR %in% c("0R", "1R")), " samples\n\n", sep = "")

# ============================================================================
# PART 5A: Overall 0R/1R Sample Comparison
# ============================================================================

cat("=== OVERALL 0R/1R SAMPLE COMPARISON ===\n")
cat("Comparing: All 0R/1R samples from 0R/1R-only patients vs 0R/1R samples from 2R patients\n\n")

# Get all samples from 0R/1R-only patients (they only have 0R/1R samples)
mpa_0r1r_only <- patients_with_0R_1R_only %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "0R/1R-only patients")

# Get ONLY 0R/1R samples from patients who develop 2R
mpa_0r1r_from_2r <- patients_with_2R %>%
    filter(ACR %in% c("0R", "1R")) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "2R patients (0R/1R samples)")

# Combine for comparison
mmf_data_0r1r <- bind_rows(mpa_0r1r_only, mpa_0r1r_from_2r)

cat("Sample sizes:\n")
cat("  0R/1R samples from 0R/1R-only patients:", nrow(mpa_0r1r_only), "samples from", 
    length(unique(mpa_0r1r_only$H)), "patients\n")
cat("    - 0R:", sum(mpa_0r1r_only$ACR == "0R"), ", 1R:", sum(mpa_0r1r_only$ACR == "1R"), "\n")
cat("  0R/1R samples from 2R patients:", nrow(mpa_0r1r_from_2r), "samples from", 
    length(unique(mpa_0r1r_from_2r$H)), "patients\n")
cat("    - 0R:", sum(mpa_0r1r_from_2r$ACR == "0R"), ", 1R:", sum(mpa_0r1r_from_2r$ACR == "1R"), "\n\n")

# Run Wilcoxon tests
cat("MPA (C18):\n")
wilcox_c18_0r1r <- wilcox.test(
    `Mycophenolate..C18.` ~ Patient_Phenotype, 
    data = mmf_data_0r1r, 
    exact = FALSE
)
print(wilcox_c18_0r1r)

cat("\nMPA (HILIC):\n")
wilcox_hilic_0r1r <- wilcox.test(
    `Mycophenolate..HILIC.` ~ Patient_Phenotype, 
    data = mmf_data_0r1r, 
    exact = FALSE
)
print(wilcox_hilic_0r1r)

# Get medians for interpretation
medians_0r1r <- mmf_data_0r1r %>%
    group_by(Patient_Phenotype) %>%
    summarise(
        n = n(),
        median_C18 = median(`Mycophenolate..C18.`, na.rm = TRUE),
        median_HILIC = median(`Mycophenolate..HILIC.`, na.rm = TRUE)
    )
print(medians_0r1r)

# Plot overall comparison
mmf_long_0r1r <- mmf_data_0r1r %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))

# Get sample sizes for x-axis labels
n_0r1r_only <- nrow(mpa_0r1r_only)
n_0r1r_from_2r <- nrow(mpa_0r1r_from_2r)

p_overall_0r1r <- ggplot(mmf_long_0r1r, aes(x = Patient_Phenotype, y = Level, fill = Patient_Phenotype)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    scale_fill_manual(values = c("0R/1R-only patients" = "lightgreen", 
                                  "2R patients (0R/1R samples)" = "lightcoral")) +
    scale_x_discrete(labels = c("0R/1R-only patients" = paste0("0R/1R-only\npatients\n(n=", n_0r1r_only, ")"),
                                "2R patients (0R/1R samples)" = paste0("2R patients\n(0R/1R samples)\n(n=", n_0r1r_from_2r, ")"))) +
    labs(title = "MPA Levels in 0R/1R Samples: Baseline State Comparison",
         subtitle = paste0("0R/1R-only patients vs 0R/1R samples from 2R patients | C18: p=", 
                          round(wilcox_c18_0r1r$p.value, 3), "; HILIC: p=", round(wilcox_hilic_0r1r$p.value, 3)),
         x = "Patient Phenotype", y = "Level") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 9))

print(p_overall_0r1r)

ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R1R_Samples_Overall.png", p_overall_0r1r, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 5B: POD-Stratified 0R/1R Sample Comparison
# ============================================================================

cat("\n\n=== POD-STRATIFIED 0R/1R SAMPLE COMPARISON ===\n\n")

# Add POD stratum
mmf_data_0r1r <- mmf_data_0r1r %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (≤30 days)", "Late (>30 days)"))

# Function to run tests for 0R/1R stratum
analyze_stratum_0r1r <- function(data, stratum) {
    cat("--- ", stratum, " ---\n")
    
    stratum_data <- data %>% filter(POD_Stratum == stratum)
    n_0r1r_only <- sum(stratum_data$Patient_Phenotype == "0R/1R-only patients")
    n_0r1r_from_2r <- sum(stratum_data$Patient_Phenotype == "2R patients (0R/1R samples)")
    
    cat("Sample size:", nrow(stratum_data), 
        "(0R/1R-only:", n_0r1r_only, ", 0R/1R from 2R patients:", n_0r1r_from_2r, ")\n")
    
    # Check if we have enough samples
    if (n_0r1r_only < 2 || n_0r1r_from_2r < 2) {
        cat("Insufficient samples for comparison\n\n")
        return(NULL)
    }
    
    # Wilcoxon tests
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ Patient_Phenotype, 
                         data = stratum_data, exact = FALSE)
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ Patient_Phenotype, 
                           data = stratum_data, exact = FALSE)
    
    cat("C18:   p =", round(w_c18$p.value, 3), "\n")
    cat("HILIC: p =", round(w_hilic$p.value, 3), "\n\n")
    
    # Get medians
    meds <- stratum_data %>%
        group_by(Patient_Phenotype) %>%
        summarise(
            median_C18 = median(`Mycophenolate..C18.`, na.rm = TRUE),
            median_HILIC = median(`Mycophenolate..HILIC.`, na.rm = TRUE)
        )
    print(meds)
    cat("\n")
    
    # Create plot
    stratum_long <- stratum_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite", values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(stratum_long, aes(x = Patient_Phenotype, y = Level, fill = Patient_Phenotype)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R/1R-only patients" = "lightgreen", 
                                      "2R patients (0R/1R samples)" = "lightcoral")) +
        scale_x_discrete(labels = c("0R/1R-only patients" = paste0("0R/1R-only\npatients\n(n=", n_0r1r_only, ")"),
                                    "2R patients (0R/1R samples)" = paste0("2R patients\n(0R/1R samples)\n(n=", n_0r1r_from_2r, ")"))) +
        labs(title = paste("MPA Levels in 0R/1R Samples:", stratum),
             subtitle = paste0("Patient phenotype comparison | C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "Patient Phenotype", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 9))
    
    return(list(plot = p, p_c18 = w_c18$p.value, p_hilic = w_hilic$p.value))
}

# Run for each stratum
early_results_0r1r <- analyze_stratum_0r1r(mmf_data_0r1r, "Early (≤30 days)")
if (!is.null(early_results_0r1r)) {
    print(early_results_0r1r$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R1R_Samples_Early.png", 
           early_results_0r1r$plot, width = 10, height = 6, dpi = 300)
}

late_results_0r1r <- analyze_stratum_0r1r(mmf_data_0r1r, "Late (>30 days)")
if (!is.null(late_results_0r1r)) {
    print(late_results_0r1r$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R1R_Samples_Late.png", 
           late_results_0r1r$plot, width = 10, height = 6, dpi = 300)
}

# ============================================================================
# PART 5C: POD-Stratified 0R/1R Analysis EXCLUDING POD<10
# ============================================================================

cat("\n\n=== POD-STRATIFIED 0R/1R SAMPLE COMPARISON (EXCLUDING POD<10) ===\n")
cat("NOTE: Removing very early post-transplant samples (POD<10)\n\n")

# Filter out POD<10 samples
mmf_0r1r_only_pod10plus <- patients_with_0R_1R_only %>%
    filter(POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "0R/1R-only patients")

mmf_0r1r_from_2r_pod10plus <- patients_with_2R %>%
    filter(ACR %in% c("0R", "1R"), POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Phenotype = "2R patients (0R/1R samples)")

mmf_data_0r1r_pod10plus <- bind_rows(mmf_0r1r_only_pod10plus, mmf_0r1r_from_2r_pod10plus)

cat("Sample sizes with POD >= 10:\n")
cat("  0R/1R samples from 0R/1R-only patients:", nrow(mmf_0r1r_only_pod10plus), "\n")
cat("  0R/1R samples from 2R patients:", nrow(mmf_0r1r_from_2r_pod10plus), "\n\n")

# Add POD stratum
mmf_data_0r1r_pod10plus <- mmf_data_0r1r_pod10plus %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (>30 days)"))

# Run for each stratum
early_results_0r1r_pod10 <- analyze_stratum_0r1r(mmf_data_0r1r_pod10plus, "Early (10-30 days)")
if (!is.null(early_results_0r1r_pod10)) {
    print(early_results_0r1r_pod10$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R1R_Samples_Early_POD10plus.png", 
           early_results_0r1r_pod10$plot, width = 10, height = 6, dpi = 300)
}

late_results_0r1r_pod10 <- analyze_stratum_0r1r(mmf_data_0r1r_pod10plus, "Late (>30 days)")
if (!is.null(late_results_0r1r_pod10)) {
    print(late_results_0r1r_pod10$plot)
    ggsave("Results2/Feedback_Analysis/Patient_Phenotype/MMF_0R1R_Samples_Late_POD10plus.png", 
           late_results_0r1r_pod10$plot, width = 10, height = 6, dpi = 300)
}

# ============================================================================
# FINAL SUMMARY INCLUDING 0R/1R COMPARISON
# ============================================================================

cat("\n\n", rep("=", 76), "\n", sep = "")
cat("=== FINAL SUMMARY (ALL COMPARISONS) ===\n")
cat(rep("=", 76), "\n\n", sep = "")

cat("--- 0R SAMPLE COMPARISON (0R-only vs 2R patients) ---\n")
cat("Overall:      C18 p=", round(wilcox_c18$p.value, 3), 
    ", HILIC p=", round(wilcox_hilic$p.value, 3), "\n")
if (!is.null(early_results)) {
    cat("Early ≤30d:   C18 p=", round(early_results$p_c18, 3), 
        ", HILIC p=", round(early_results$p_hilic, 3), "\n")
}
if (!is.null(late_results)) {
    cat("Late >30d:    C18 p=", round(late_results$p_c18, 3), 
        ", HILIC p=", round(late_results$p_hilic, 3), "\n")
}

cat("\n--- 0R/1R SAMPLE COMPARISON (0R/1R-only vs 2R patients) ---\n")
cat("Overall:      C18 p=", round(wilcox_c18_0r1r$p.value, 3), 
    ", HILIC p=", round(wilcox_hilic_0r1r$p.value, 3), "\n")
if (!is.null(early_results_0r1r)) {
    cat("Early ≤30d:   C18 p=", round(early_results_0r1r$p_c18, 3), 
        ", HILIC p=", round(early_results_0r1r$p_hilic, 3), "\n")
}
if (!is.null(late_results_0r1r)) {
    cat("Late >30d:    C18 p=", round(late_results_0r1r$p_c18, 3), 
        ", HILIC p=", round(late_results_0r1r$p_hilic, 3), "\n")
}

cat("\nAll plots saved to: Results2/Feedback_Analysis/Patient_Phenotype/\n\n")

cat("=== KEY INSIGHTS ===\n\n")
cat("Comparing 0R samples only:\n")
cat("  - Strictest comparison (0R-only patients vs 0R from 2R patients)\n")
cat("  - Controls for both rejection state AND mild rejection (1R)\n")
cat("  - Tests pure patient phenotype during complete quiescence\n\n")

cat("Comparing 0R/1R samples:\n")
cat("  - Broader comparison (includes mild rejection)\n")
cat("  - More statistical power (larger sample size)\n")
cat("  - Tests patient phenotype during non-severe rejection states\n")
cat("  - Clinically relevant: 1R is often not treated aggressively\n\n")

cat("If both analyses show similar patterns:\n")
cat("  → Strong evidence for inherent patient-level differences\n\n")
cat("If 0R/1R analysis is significant but 0R-only is not:\n")
cat("  → 1R may contribute to or reflect underlying patient differences\n\n")
cat("If 0R-only is significant but 0R/1R is not:\n")
cat("  → Including 1R may dilute true quiescent-state differences\n\n")
