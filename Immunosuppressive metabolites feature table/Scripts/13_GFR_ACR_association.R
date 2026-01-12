# GFR and ACR Association Analysis
# Testing whether kidney function (GFR) is associated with rejection status (ACR)

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# ============================================================================
# LOAD AND PREPARE DATA
# ============================================================================

cat("\n=== LOADING CLINICAL DATA ===\n")

# Load clinical data
clinical_data <- read_csv("OHT_Clinical.csv")

cat("Clinical data loaded:", nrow(clinical_data), "patients\n\n")

# Define GFR categories (same as script 12)
categorize_gfr <- function(gfr) {
    case_when(
        is.na(gfr) ~ NA_character_,
        gfr >= 90 ~ "Normal (≥90)",
        gfr >= 60 ~ "Mild Reduction (60-89)",
        gfr >= 30 ~ "Moderate Reduction (30-59)",
        TRUE ~ "Severe (<30)"
    )
}

clinical_data <- clinical_data %>%
    mutate(GFR_Category = categorize_gfr(GFR))

# Merge with metabolite data
mmf_gfr_data <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, S, Sample_ID) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R")) %>%
    left_join(clinical_data %>% select(H, GFR, GFR_Category, Age, Sex, Race, BMI), by = "H")

cat("Total samples:", nrow(mmf_gfr_data), "\n")
cat("Samples with GFR data:", sum(!is.na(mmf_gfr_data$GFR)), "\n\n")

# ============================================================================
# PART 1: DESCRIPTIVE STATISTICS
# ============================================================================

cat("========================================================================\n")
cat("PART 1: GFR DESCRIPTIVE STATISTICS BY ACR GROUP\n")
cat("========================================================================\n\n")

# Overall GFR summary
cat("--- Overall GFR Summary ---\n")
cat("Mean GFR:", round(mean(mmf_gfr_data$GFR, na.rm = TRUE), 2), "mL/min/1.73m²\n")
cat("Median GFR:", round(median(mmf_gfr_data$GFR, na.rm = TRUE), 2), "mL/min/1.73m²\n")
cat("SD:", round(sd(mmf_gfr_data$GFR, na.rm = TRUE), 2), "\n")
cat("Range:", round(min(mmf_gfr_data$GFR, na.rm = TRUE), 1), "-", 
    round(max(mmf_gfr_data$GFR, na.rm = TRUE), 1), "\n\n")

# GFR by ACR group
cat("--- GFR by ACR Group ---\n")
gfr_by_acr <- mmf_gfr_data %>%
    filter(!is.na(GFR)) %>%
    group_by(ACR_Group) %>%
    summarise(
        n = n(),
        Mean = round(mean(GFR), 2),
        Median = round(median(GFR), 2),
        SD = round(sd(GFR), 2),
        Min = round(min(GFR), 1),
        Max = round(max(GFR), 1),
        Q1 = round(quantile(GFR, 0.25), 1),
        Q3 = round(quantile(GFR, 0.75), 1)
    )

print(gfr_by_acr)
cat("\n")

# ============================================================================
# PART 2: STATISTICAL TEST - GFR COMPARISON BETWEEN ACR GROUPS
# ============================================================================

cat("========================================================================\n")
cat("PART 2: GFR COMPARISON BETWEEN ACR GROUPS\n")
cat("========================================================================\n\n")

cat("Research Question: Do patients with 2R+ have different GFR than those with 0R?\n\n")

# Wilcoxon rank-sum test
cat("--- Wilcoxon Rank-Sum Test ---\n")
wilcox_result <- wilcox.test(GFR ~ ACR_Group, 
                              data = mmf_gfr_data %>% filter(!is.na(GFR)), 
                              exact = FALSE)
print(wilcox_result)

cat("\nInterpretation:\n")
if(wilcox_result$p.value < 0.05) {
    cat("✓ Significant difference in GFR between 0R and 2R+ groups (p =", 
        round(wilcox_result$p.value, 4), ")\n")
    if(median(mmf_gfr_data$GFR[mmf_gfr_data$ACR_Group == "0R"], na.rm = TRUE) > 
       median(mmf_gfr_data$GFR[mmf_gfr_data$ACR_Group == "2R+"], na.rm = TRUE)) {
        cat("  → 0R samples have higher GFR than 2R+ samples\n")
    } else {
        cat("  → 2R+ samples have higher GFR than 0R samples\n")
    }
} else {
    cat("✗ No significant difference in GFR between groups (p =", 
        round(wilcox_result$p.value, 4), ")\n")
}
cat("\n")

# T-test for comparison
cat("--- T-test (parametric comparison) ---\n")
t_result <- t.test(GFR ~ ACR_Group, data = mmf_gfr_data %>% filter(!is.na(GFR)))
print(t_result)
cat("\n")

# ============================================================================
# PART 3: CONTINGENCY TABLE - GFR CATEGORY BY ACR GROUP
# ============================================================================

cat("========================================================================\n")
cat("PART 3: GFR CATEGORY DISTRIBUTION BY ACR GROUP\n")
cat("========================================================================\n\n")

cat("Research Question: Is the distribution of GFR categories different between 0R and 2R+?\n\n")

# Create contingency table
cat("--- Contingency Table: GFR Category × ACR Group ---\n")
gfr_order <- c("Normal (≥90)", "Mild Reduction (60-89)", 
               "Moderate Reduction (30-59)", "Severe (<30)")

contingency_data <- mmf_gfr_data %>%
    filter(!is.na(GFR_Category)) %>%
    mutate(GFR_Category = factor(GFR_Category, levels = gfr_order))

cont_table <- table(contingency_data$GFR_Category, contingency_data$ACR_Group)
print(cont_table)

cat("\n--- Proportions by ACR Group ---\n")
prop_table <- prop.table(cont_table, margin = 2) * 100  # Proportions by column (ACR group)
print(round(prop_table, 1))
cat("\n")

# Chi-square test
cat("--- Chi-Square Test ---\n")
if(all(cont_table >= 5)) {
    chisq_result <- chisq.test(cont_table)
    print(chisq_result)
    cat("\nInterpretation:\n")
    if(chisq_result$p.value < 0.05) {
        cat("✓ Significant association between GFR category and ACR status (p =", 
            round(chisq_result$p.value, 4), ")\n")
        cat("  → The distribution of GFR categories differs between 0R and 2R+ groups\n")
    } else {
        cat("✗ No significant association (p =", round(chisq_result$p.value, 4), ")\n")
        cat("  → GFR category distribution is similar between 0R and 2R+ groups\n")
    }
} else {
    cat("WARNING: Some expected frequencies < 5, using Fisher's Exact Test instead\n\n")
    fisher_result <- fisher.test(cont_table)
    print(fisher_result)
    cat("\nInterpretation:\n")
    if(fisher_result$p.value < 0.05) {
        cat("✓ Significant association between GFR category and ACR status (p =", 
            round(fisher_result$p.value, 4), ")\n")
    } else {
        cat("✗ No significant association (p =", round(fisher_result$p.value, 4), ")\n")
    }
}
cat("\n")

# ============================================================================
# PART 4: VISUALIZATION
# ============================================================================

cat("========================================================================\n")
cat("PART 4: VISUALIZATIONS\n")
cat("========================================================================\n\n")

# Create output directory
dir.create("Results2/Feedback_Analysis/GFR_ACR_Association", showWarnings = FALSE, recursive = TRUE)

# 1. Boxplot: GFR by ACR group
p_boxplot <- ggplot(mmf_gfr_data %>% filter(!is.na(GFR)), 
                    aes(x = ACR_Group, y = GFR, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    labs(title = "GFR Distribution by ACR Status",
         subtitle = paste0("Wilcoxon p = ", round(wilcox_result$p.value, 4)),
         x = "ACR Group",
         y = "GFR (mL/min/1.73m²)") +
    theme_minimal() +
    theme(legend.position = "none")

print(p_boxplot)
ggsave("Results2/Feedback_Analysis/GFR_ACR_Association/GFR_by_ACR_Boxplot.png", 
       p_boxplot, width = 8, height = 6, dpi = 300)
cat("Saved: GFR_by_ACR_Boxplot.png\n")

# 2. Violin plot with sample sizes
n_0r <- sum(mmf_gfr_data$ACR_Group == "0R" & !is.na(mmf_gfr_data$GFR))
n_2r <- sum(mmf_gfr_data$ACR_Group == "2R+" & !is.na(mmf_gfr_data$GFR))

p_violin <- ggplot(mmf_gfr_data %>% filter(!is.na(GFR)), 
                   aes(x = ACR_Group, y = GFR, fill = ACR_Group)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.1, alpha = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r, ")"),
                                "2R+" = paste0("2R+\n(n=", n_2r, ")"))) +
    labs(title = "GFR Distribution by ACR Status (Violin Plot)",
         subtitle = paste0("Median 0R: ", round(median(mmf_gfr_data$GFR[mmf_gfr_data$ACR_Group == "0R"], na.rm = TRUE), 1),
                          " | Median 2R+: ", round(median(mmf_gfr_data$GFR[mmf_gfr_data$ACR_Group == "2R+"], na.rm = TRUE), 1),
                          " | p = ", round(wilcox_result$p.value, 4)),
         x = "ACR Group",
         y = "GFR (mL/min/1.73m²)") +
    theme_minimal() +
    theme(legend.position = "none")

print(p_violin)
ggsave("Results2/Feedback_Analysis/GFR_ACR_Association/GFR_by_ACR_Violin.png", 
       p_violin, width = 8, height = 6, dpi = 300)
cat("Saved: GFR_by_ACR_Violin.png\n")

# 3. Stacked bar chart: GFR category by ACR group
gfr_prop_data <- contingency_data %>%
    group_by(ACR_Group, GFR_Category) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(ACR_Group) %>%
    mutate(prop = n / sum(n) * 100)

p_stacked <- ggplot(gfr_prop_data, aes(x = ACR_Group, y = prop, fill = GFR_Category)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
    geom_text(aes(label = paste0(round(prop, 1), "%\n(n=", n, ")")), 
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    labs(title = "GFR Category Distribution by ACR Status",
         x = "ACR Group",
         y = "Percentage (%)",
         fill = "GFR Category") +
    theme_minimal() +
    theme(legend.position = "right")

print(p_stacked)
ggsave("Results2/Feedback_Analysis/GFR_ACR_Association/GFR_Category_by_ACR_Stacked.png", 
       p_stacked, width = 10, height = 6, dpi = 300)
cat("Saved: GFR_Category_by_ACR_Stacked.png\n")

# 4. Side-by-side bar chart
p_sidebyside <- ggplot(gfr_prop_data, aes(x = GFR_Category, y = prop, fill = ACR_Group)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(aes(label = paste0(round(prop, 1), "%")), 
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    labs(title = "GFR Category Distribution: 0R vs 2R+",
         x = "GFR Category",
         y = "Percentage (%)",
         fill = "ACR Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_sidebyside)
ggsave("Results2/Feedback_Analysis/GFR_ACR_Association/GFR_Category_by_ACR_SideBySide.png", 
       p_sidebyside, width = 10, height = 6, dpi = 300)
cat("Saved: GFR_Category_by_ACR_SideBySide.png\n\n")

# ============================================================================
# PART 5: PATIENT-LEVEL ANALYSIS
# ============================================================================

cat("========================================================================\n")
cat("PART 5: PATIENT-LEVEL GFR COMPARISON\n")
cat("========================================================================\n\n")

cat("NOTE: Above analyses use sample-level data (repeated measures per patient)\n")
cat("      This section uses patient-level data (one GFR per patient)\n\n")

# Create patient-level dataset
patient_level <- mmf_gfr_data %>%
    group_by(H) %>%
    summarise(
        GFR = first(GFR),  # GFR is patient-level (same for all samples)
        GFR_Category = first(GFR_Category),
        Has_2R = any(ACR_Group == "2R+"),
        .groups = "drop"
    ) %>%
    mutate(Patient_Group = ifelse(Has_2R, "Patients with 2R+", "Patients with only 0R"))

cat("Total patients:", nrow(patient_level), "\n")
cat("Patients with 2R+:", sum(patient_level$Has_2R), "\n")
cat("Patients with only 0R:", sum(!patient_level$Has_2R), "\n\n")

# Patient-level GFR comparison
cat("--- Patient-Level GFR by Rejection History ---\n")
patient_gfr_summary <- patient_level %>%
    filter(!is.na(GFR)) %>%
    group_by(Patient_Group) %>%
    summarise(
        n = n(),
        Mean = round(mean(GFR), 2),
        Median = round(median(GFR), 2),
        SD = round(sd(GFR), 2),
        Min = round(min(GFR), 1),
        Max = round(max(GFR), 1)
    )

print(patient_gfr_summary)
cat("\n")

# Patient-level statistical test
cat("--- Wilcoxon Test (Patient-Level) ---\n")
wilcox_patient <- wilcox.test(GFR ~ Patient_Group, 
                               data = patient_level %>% filter(!is.na(GFR)), 
                               exact = FALSE)
print(wilcox_patient)

cat("\nInterpretation:\n")
if(wilcox_patient$p.value < 0.05) {
    cat("✓ Patients who develop 2R+ have significantly different GFR (p =", 
        round(wilcox_patient$p.value, 4), ")\n")
} else {
    cat("✗ No significant difference in patient-level GFR by rejection history (p =", 
        round(wilcox_patient$p.value, 4), ")\n")
}
cat("\n")

# Patient-level visualization
p_patient <- ggplot(patient_level %>% filter(!is.na(GFR)), 
                    aes(x = Patient_Group, y = GFR, fill = Patient_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    scale_fill_manual(values = c("Patients with only 0R" = "lightblue", 
                                  "Patients with 2R+" = "lightcoral")) +
    labs(title = "Patient-Level GFR by Rejection History",
         subtitle = paste0("p = ", round(wilcox_patient$p.value, 4)),
         x = "Patient Group",
         y = "GFR (mL/min/1.73m²)") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 9))

print(p_patient)
ggsave("Results2/Feedback_Analysis/GFR_ACR_Association/GFR_by_Patient_Rejection_History.png", 
       p_patient, width = 8, height = 6, dpi = 300)
cat("Saved: GFR_by_Patient_Rejection_History.png\n\n")

# ============================================================================
# SUMMARY AND INTERPRETATION
# ============================================================================

cat("========================================================================\n")
cat("SUMMARY AND CLINICAL INTERPRETATION\n")
cat("========================================================================\n\n")

cat("=== Key Findings ===\n\n")

cat("1. SAMPLE-LEVEL COMPARISON (0R vs 2R+ samples):\n")
cat("   - Wilcoxon p-value:", round(wilcox_result$p.value, 4), "\n")
cat("   - Median GFR (0R):", round(median(mmf_gfr_data$GFR[mmf_gfr_data$ACR_Group == "0R"], na.rm = TRUE), 1), 
    "mL/min/1.73m²\n")
cat("   - Median GFR (2R+):", round(median(mmf_gfr_data$GFR[mmf_gfr_data$ACR_Group == "2R+"], na.rm = TRUE), 1), 
    "mL/min/1.73m²\n\n")

cat("2. PATIENT-LEVEL COMPARISON (by rejection history):\n")
cat("   - Wilcoxon p-value:", round(wilcox_patient$p.value, 4), "\n")
cat("   - Median GFR (0R-only patients):", 
    round(median(patient_level$GFR[patient_level$Patient_Group == "Patients with only 0R"], na.rm = TRUE), 1), 
    "mL/min/1.73m²\n")
cat("   - Median GFR (2R+ patients):", 
    round(median(patient_level$GFR[patient_level$Patient_Group == "Patients with 2R+"], na.rm = TRUE), 1), 
    "mL/min/1.73m²\n\n")

cat("=== Clinical Implications ===\n\n")

if(wilcox_result$p.value < 0.05 | wilcox_patient$p.value < 0.05) {
    cat("⚠️  CONFOUNDING ALERT:\n")
    cat("   - GFR differs between ACR groups\n")
    cat("   - GFR-stratified analyses in script 12 are IMPORTANT to control for this\n")
    cat("   - Any observed MPA differences by ACR could be partly due to GFR differences\n\n")
    
    cat("Recommendations:\n")
    cat("   1. Report GFR differences in baseline tables\n")
    cat("   2. Use GFR-stratified analyses (script 12) as PRIMARY analysis\n")
    cat("   3. Consider multivariable models adjusting for GFR\n")
    cat("   4. Discuss GFR as potential confounder in interpretation\n\n")
} else {
    cat("✓ GOOD NEWS:\n")
    cat("   - GFR does not differ significantly between ACR groups\n")
    cat("   - GFR is unlikely to be a major confounder\n")
    cat("   - Overall MPA comparisons are less likely confounded by kidney function\n\n")
    
    cat("Note:\n")
    cat("   - Still valuable to report GFR-stratified analyses for completeness\n")
    cat("   - Shows results are consistent across kidney function levels\n\n")
}

cat("=== Files Generated ===\n")
cat("Results2/Feedback_Analysis/GFR_ACR_Association/\n")
cat("  - GFR_by_ACR_Boxplot.png\n")
cat("  - GFR_by_ACR_Violin.png\n")
cat("  - GFR_Category_by_ACR_Stacked.png\n")
cat("  - GFR_Category_by_ACR_SideBySide.png\n")
cat("  - GFR_by_Patient_Rejection_History.png\n\n")

cat("========================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================================================\n\n")
