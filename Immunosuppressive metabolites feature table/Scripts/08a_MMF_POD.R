# MMF Analysis with POD Stratification - Between-Patient Comparison
# Comparing samples from patients who NEVER reject (0R/1R only) vs patients who DO reject (2R+)
# This is a BETWEEN-PATIENT comparison (different patient populations)

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# NOTE: patients_with_0R_1R_only and patients_with_2R are already created in 02_setup

# ============================================================================
# VERIFY PATIENT GROUPS
# ============================================================================

cat("\n=== PATIENT GROUPS (from 02_setup) ===\n\n")

cat("Patients who NEVER develop 2R (0R/1R only):", 
    length(unique(patients_with_0R_1R_only$H)), "patients,", 
    nrow(patients_with_0R_1R_only), "samples\n")
cat("Patients who develop 2R at any point:", 
    length(unique(patients_with_2R$H)), "patients,", 
    nrow(patients_with_2R), "samples\n\n")

# ============================================================================
# PART 1: Overall MMF comparison (0R/1R-only patients vs 2R+ samples)
# ============================================================================

cat("\n=== OVERALL MMF ANALYSIS (Between-Patient Comparison) ===\n")
cat("Comparing: All samples from 0R/1R-only patients vs 2R+ samples from patients who reject\n\n")

# Get data from 0R/1R-only patients (all their samples)
mmf_0r1r_only <- patients_with_0R_1R_only %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "0R/1R-only patients")

# Get only 2R+ samples from patients who develop rejection
mmf_2r_samples <- patients_with_2R %>%
    filter(grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "2R+ samples")

# Combine for comparison
mmf_data <- bind_rows(mmf_0r1r_only, mmf_2r_samples)

cat("Total samples in comparison:\n")
cat("  0R/1R-only patients:", nrow(mmf_0r1r_only), "samples from", 
    length(unique(mmf_0r1r_only$H)), "patients\n")
cat("  2R+ samples:", nrow(mmf_2r_samples), "samples from", 
    length(unique(mmf_2r_samples$H)), "patients\n\n")

# Run Wilcoxon tests
cat("Mycophenolate C18:\n")
wilcox_c18 <- wilcox.test(
    `Mycophenolate..C18.` ~ Patient_Group, 
    data = mmf_data, 
    exact = FALSE
)
print(wilcox_c18)

cat("\nMycophenolate HILIC:\n")
wilcox_hilic <- wilcox.test(
    `Mycophenolate..HILIC.` ~ Patient_Group, 
    data = mmf_data, 
    exact = FALSE
)
print(wilcox_hilic)

# Plot overall comparison
mmf_long <- mmf_data %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))

# Get sample sizes for x-axis labels
n_0r1r <- nrow(mmf_0r1r_only)
n_2r <- nrow(mmf_2r_samples)

p_overall <- ggplot(mmf_long, aes(x = Patient_Group, y = Level, fill = Patient_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    scale_fill_manual(values = c("0R/1R-only patients" = "lightblue", "2R+ samples" = "lightcoral")) +
    scale_x_discrete(labels = c("0R/1R-only patients" = paste0("0R/1R-only\n(n=", n_0r1r, ")"),
                                "2R+ samples" = paste0("2R+\n(n=", n_2r, ")"))) +
    labs(title = "MMF Levels: 0R/1R-only Patients vs 2R+ Samples",
         subtitle = paste0("Between-patient comparison | C18: p=", round(wilcox_c18$p.value, 3), 
                          "; HILIC: p=", round(wilcox_hilic$p.value, 3)),
         x = "Patient Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "none")

print(p_overall)

# Create output directory
dir.create("Results2/Feedback_Analysis/Between_Patient_POD_Stratified", showWarnings = FALSE, recursive = TRUE)
ggsave("Results2/Feedback_Analysis/Between_Patient_POD_Stratified/MMF_Overall.png", p_overall, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 2: POD-Stratified Analysis (Early ≤30 days vs Late >30 days)
# ============================================================================

cat("\n\n=== POD-STRATIFIED MMF ANALYSIS ===\n\n")

# Add POD stratum
mmf_data <- mmf_data %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (≤30 days)", "Late (>30 days)"))

# Function to run tests for one stratum
analyze_stratum <- function(data, stratum) {
    cat("--- ", stratum, " ---\n")
    
    stratum_data <- data %>% filter(POD_Stratum == stratum)
    n_0r1r <- sum(stratum_data$Patient_Group == "0R/1R-only patients")
    n_2r <- sum(stratum_data$Patient_Group == "2R+ samples")
    
    cat("Sample size:", nrow(stratum_data), 
        "(0R/1R-only:", n_0r1r, ", 2R+:", n_2r, ")\n\n")
    
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
        scale_fill_manual(values = c("0R/1R-only patients" = "lightblue", "2R+ samples" = "lightcoral")) +
        scale_x_discrete(labels = c("0R/1R-only patients" = paste0("0R/1R-only\n(n=", n_0r1r, ")"),
                                    "2R+ samples" = paste0("2R+\n(n=", n_2r, ")"))) +
        labs(title = paste("MMF Levels:", stratum, "(Between-Patient)"),
             subtitle = paste0("C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "Patient Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(plot = p, p_c18 = w_c18$p.value, p_hilic = w_hilic$p.value))
}

# Run for each stratum
early_results <- analyze_stratum(mmf_data, "Early (≤30 days)")
print(early_results$plot)
ggsave("Results2/Feedback_Analysis/Between_Patient_POD_Stratified/MMF_Early.png", early_results$plot, 
       width = 10, height = 6, dpi = 300)

late_results <- analyze_stratum(mmf_data, "Late (>30 days)")
print(late_results$plot)
ggsave("Results2/Feedback_Analysis/Between_Patient_POD_Stratified/MMF_Late.png", late_results$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 3: POD-Stratified Analysis EXCLUDING POD<10
# ============================================================================

cat("\n\n=== POD-STRATIFIED MMF ANALYSIS (EXCLUDING POD<10) ===\n")
cat("NOTE: Removing very early post-transplant samples (POD<10)\n\n")

# Filter out POD<10 samples
mmf_0r1r_pod10plus <- patients_with_0R_1R_only %>%
    filter(POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "0R/1R-only patients")

mmf_2r_pod10plus <- patients_with_2R %>%
    filter(grepl("^2R", ACR, ignore.case = TRUE)) %>%
    filter(POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "2R+ samples")

mmf_data_pod10plus <- bind_rows(mmf_0r1r_pod10plus, mmf_2r_pod10plus)

cat("Total samples with POD >= 10:", nrow(mmf_data_pod10plus), "\n")
cat("0R/1R-only patients:", nrow(mmf_0r1r_pod10plus), "\n")
cat("2R+ samples:", nrow(mmf_2r_pod10plus), "\n\n")

# Add POD stratum
mmf_data_pod10plus <- mmf_data_pod10plus %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (>30 days)"))

# Run for each stratum
early_results_pod10 <- analyze_stratum(mmf_data_pod10plus, "Early (10-30 days)")
print(early_results_pod10$plot)
ggsave("Results2/Feedback_Analysis/Between_Patient_POD_Stratified/MMF_Early_POD10plus.png", early_results_pod10$plot, 
       width = 10, height = 6, dpi = 300)

late_results_pod10 <- analyze_stratum(mmf_data_pod10plus, "Late (>30 days)")
print(late_results_pod10$plot)
ggsave("Results2/Feedback_Analysis/Between_Patient_POD_Stratified/MMF_Late_POD10plus.png", late_results_pod10$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 4: POD-Stratified Analysis EXCLUDING POD<10 AND POD>45
# ============================================================================

cat("\n\n=== POD-STRATIFIED MMF ANALYSIS (POD 10-45 ONLY) ===\n")
cat("NOTE: Removing very early (POD<10) and late (POD>45) samples\n")
cat("      to focus on the acute post-transplant window\n\n")

# Filter to POD 10-45 only
mmf_0r1r_pod10_45 <- patients_with_0R_1R_only %>%
    filter(POD >= 10 & POD <= 45) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "0R/1R-only patients")

mmf_2r_pod10_45 <- patients_with_2R %>%
    filter(grepl("^2R", ACR, ignore.case = TRUE)) %>%
    filter(POD >= 10 & POD <= 45) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(Patient_Group = "2R+ samples")

mmf_data_pod10_45 <- bind_rows(mmf_0r1r_pod10_45, mmf_2r_pod10_45)

cat("Total samples with POD 10-45:", nrow(mmf_data_pod10_45), "\n")
cat("0R/1R-only patients:", nrow(mmf_0r1r_pod10_45), "\n")
cat("2R+ samples:", nrow(mmf_2r_pod10_45), "\n\n")

# Add POD stratum
mmf_data_pod10_45 <- mmf_data_pod10_45 %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (31-45 days)"))

# Run for each stratum
early_results_pod10_45 <- analyze_stratum(mmf_data_pod10_45, "Early (10-30 days)")
print(early_results_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/Between_Patient_POD_Stratified/MMF_Early_POD10_45.png", early_results_pod10_45$plot, 
       width = 10, height = 6, dpi = 300)

late_results_pod10_45 <- analyze_stratum(mmf_data_pod10_45, "Late (31-45 days)")
print(late_results_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/Between_Patient_POD_Stratified/MMF_Late_POD10_45.png", late_results_pod10_45$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== SUMMARY (BETWEEN-PATIENT COMPARISON) ===\n")
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
cat("Early 10-45d: C18 p=", round(early_results_pod10_45$p_c18, 3), 
    ", HILIC p=", round(early_results_pod10_45$p_hilic, 3), "\n")
cat("Late 31-45d:  C18 p=", round(late_results_pod10_45$p_c18, 3), 
    ", HILIC p=", round(late_results_pod10_45$p_hilic, 3), "\n")
cat("\nPlots saved to: Results2/Feedback_Analysis/Between_Patient_POD_Stratified/\n\n")

cat("\n=== INTERPRETATION ===\n")
cat("This analysis compares TWO DIFFERENT PATIENT POPULATIONS:\n")
cat("  - 0R/1R-only patients: Never develop clinically significant rejection\n")
cat("  - 2R+ samples: From patients who DO develop clinically significant rejection\n")
cat("\nThis differs from 08_MMF_POD.R which compares samples WITHIN the same\n")
cat("patient population (patients who eventually develop 2R).\n\n")



