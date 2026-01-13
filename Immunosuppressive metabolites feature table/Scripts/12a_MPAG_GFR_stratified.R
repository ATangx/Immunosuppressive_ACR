# GFR-Stratified Analysis: 0R vs 2R+ Mycophenolic Acid O-Acyl Glucuronide Levels
# NOTE: MPAG is the major metabolite of MPA and is RENALLY CLEARED
#       Expected to show stronger GFR effects than parent drug (MPA)

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# NOTE: Comparing 0R vs 2R+ only (excluding 1R for clearer comparison)

# ============================================================================
# LOAD CLINICAL DATA WITH GFR INFORMATION
# ============================================================================

cat("\n=== LOADING CLINICAL DATA ===\n")

# Load clinical data
clinical_data <- read_csv("OHT_Clinical.csv")

cat("Clinical data loaded:", nrow(clinical_data), "patients\n")
cat("\nGFR summary statistics:\n")
print(summary(clinical_data$GFR))

# ============================================================================
# DEFINE CLINICALLY RELEVANT GFR CATEGORIES
# ============================================================================

# Clinical GFR categories based on CKD staging (4 groups):
# Normal (≥90): CKD Stage 1 (normal kidney function)
# Mild (60-89): CKD Stage 2 (mildly decreased)
# Moderate (30-59): CKD Stage 3 (moderate reduction, combines 3a and 3b)
# Severe (<30): CKD Stage 4-5 (severe reduction/kidney failure)

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

cat("\nGFR Category distribution in clinical data:\n")
print(table(clinical_data$GFR_Category, useNA = "ifany"))

# ============================================================================
# MERGE METABOLITE DATA WITH GFR INFORMATION
# ============================================================================

cat("\n=== MERGING METABOLITE AND CLINICAL DATA ===\n\n")

# Get MPAG data and merge with GFR
mpag_gfr_data <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, `Mycophenolic.acid.O.acyl.glucuronide..C18.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R")) %>%
    left_join(clinical_data %>% select(H, GFR, GFR_Category, Age, Sex, Race, BMI), by = "H")

cat("Total samples after merge:", nrow(mpag_gfr_data), "\n")
cat("Samples with GFR information:", sum(!is.na(mpag_gfr_data$GFR)), "\n")
cat("Samples with MPAG data:", sum(!is.na(mpag_gfr_data$`Mycophenolic.acid.O.acyl.glucuronide..C18.`)), "\n\n")

# Check GFR category distribution in merged data
cat("GFR Category distribution in merged data:\n")
print(table(mpag_gfr_data$GFR_Category, useNA = "ifany"))

cat("\nGFR Category by ACR group:\n")
print(table(mpag_gfr_data$GFR_Category, mpag_gfr_data$ACR_Group, useNA = "ifany"))

# ============================================================================
# GFR-STRATIFIED ANALYSIS
# ============================================================================

cat("\n\n=== GFR-STRATIFIED MPAG ANALYSIS (0R vs 2R+) ===\n\n")
cat("NOTE: MPAG is renally cleared - expect stronger GFR effects!\n\n")

# Function to analyze MPAG levels by GFR category
analyze_mpag_by_gfr <- function(data, gfr_category) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("GFR CATEGORY:", gfr_category, "\n")
    cat(rep("=", 70), "\n", sep = "")
    
    # Filter for this GFR category
    gfr_data <- data %>% filter(GFR_Category == gfr_category)
    
    n_total <- nrow(gfr_data)
    n_0r <- sum(gfr_data$ACR_Group == "0R")
    n_2r <- sum(gfr_data$ACR_Group == "2R+")
    
    cat("Total samples:", n_total, "(0R:", n_0r, ", 2R+:", n_2r, ")\n")
    
    if(n_total > 0) {
        cat("GFR range in this category:", 
            round(min(gfr_data$GFR, na.rm = TRUE), 1), "-", 
            round(max(gfr_data$GFR, na.rm = TRUE), 1), "\n\n")
    }
    
    # Check if we have enough samples
    if(n_0r < 3 | n_2r < 3) {
        cat("WARNING: Insufficient sample size for analysis (need ≥3 per group)\n")
        return(NULL)
    }
    
    # Wilcoxon test for MPAG
    cat("--- Mycophenolic Acid O-Acyl Glucuronide C18 ---\n")
    w_mpag <- wilcox.test(`Mycophenolic.acid.O.acyl.glucuronide..C18.` ~ ACR_Group, 
                          data = gfr_data, exact = FALSE)
    cat("p-value:", round(w_mpag$p.value, 4), "\n")
    cat("Median 0R:", round(median(gfr_data$`Mycophenolic.acid.O.acyl.glucuronide..C18.`[gfr_data$ACR_Group == "0R"], na.rm = TRUE), 3), "\n")
    cat("Median 2R+:", round(median(gfr_data$`Mycophenolic.acid.O.acyl.glucuronide..C18.`[gfr_data$ACR_Group == "2R+"], na.rm = TRUE), 3), "\n\n")
    
    # Create plot
    p <- ggplot(gfr_data, aes(x = ACR_Group, y = `Mycophenolic.acid.O.acyl.glucuronide..C18.`, fill = ACR_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
        scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r, ")"),
                                    "2R+" = paste0("2R+\n(n=", n_2r, ")"))) +
        labs(title = paste("MPAG Levels:", gfr_category, "(0R vs 2R+)"),
             subtitle = paste0("p=", round(w_mpag$p.value, 3)),
             x = "ACR Group", y = "MPAG C18 Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(
        plot = p, 
        gfr_category = gfr_category,
        n_0r = n_0r,
        n_2r = n_2r,
        p_value = w_mpag$p.value,
        median_0r = median(gfr_data$`Mycophenolic.acid.O.acyl.glucuronide..C18.`[gfr_data$ACR_Group == "0R"], na.rm = TRUE),
        median_2r = median(gfr_data$`Mycophenolic.acid.O.acyl.glucuronide..C18.`[gfr_data$ACR_Group == "2R+"], na.rm = TRUE),
        gfr_min = min(gfr_data$GFR, na.rm = TRUE),
        gfr_max = max(gfr_data$GFR, na.rm = TRUE)
    ))
}

# Get unique GFR categories with data (in logical order)
gfr_order <- c("Normal (≥90)", "Mild Reduction (60-89)", 
               "Moderate Reduction (30-59)", "Severe (<30)")
gfr_categories <- gfr_order[gfr_order %in% unique(mpag_gfr_data$GFR_Category[!is.na(mpag_gfr_data$GFR_Category)])]

cat("GFR categories to analyze:", paste(gfr_categories, collapse = ", "), "\n")

# Create output directory
dir.create("Results2/Feedback_Analysis/MPAG_GFR_Stratified", showWarnings = FALSE, recursive = TRUE)

# Run analysis for each GFR category
mpag_gfr_results <- list()
for(gfr_cat in gfr_categories) {
    result <- analyze_mpag_by_gfr(mpag_gfr_data, gfr_cat)
    if(!is.null(result)) {
        mpag_gfr_results[[gfr_cat]] <- result
    }
}

# Save all plots
cat("\n=== SAVING PLOTS ===\n")
for(gfr_cat in names(mpag_gfr_results)) {
    filename <- paste0("Results2/Feedback_Analysis/MPAG_GFR_Stratified/MPAG_", 
                      gsub("[^A-Za-z0-9]", "_", gfr_cat), "_0R_vs_2R.png")
    ggsave(filename, mpag_gfr_results[[gfr_cat]]$plot, width = 8, height = 6, dpi = 300)
    cat("Saved:", filename, "\n")
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n\n=== GFR-STRATIFIED MPAG SUMMARY ===\n\n")

if(length(mpag_gfr_results) > 0) {
    summary_df <- data.frame(
        GFR_Category = sapply(mpag_gfr_results, function(x) x$gfr_category),
        GFR_Range = paste0(sapply(mpag_gfr_results, function(x) round(x$gfr_min, 1)), 
                          " - ", 
                          sapply(mpag_gfr_results, function(x) round(x$gfr_max, 1))),
        N_0R = sapply(mpag_gfr_results, function(x) x$n_0r),
        N_2R_plus = sapply(mpag_gfr_results, function(x) x$n_2r),
        Median_0R = round(sapply(mpag_gfr_results, function(x) x$median_0r), 0),
        Median_2R = round(sapply(mpag_gfr_results, function(x) x$median_2r), 0),
        P_value = round(sapply(mpag_gfr_results, function(x) x$p_value), 4)
    )
    
    # Reorder by GFR category
    summary_df$GFR_Category <- factor(summary_df$GFR_Category, levels = gfr_order)
    summary_df <- summary_df %>% arrange(GFR_Category)
    summary_df$GFR_Category <- as.character(summary_df$GFR_Category)
    
    print(summary_df)
    
    # Save summary table
    write_csv(summary_df, "Results2/Feedback_Analysis/MPAG_GFR_Stratified/MPAG_GFR_Stratified_Summary.csv")
    cat("\nSummary table saved to: Results2/Feedback_Analysis/MPAG_GFR_Stratified/MPAG_GFR_Stratified_Summary.csv\n")
}

# ============================================================================
# COMBINED VISUALIZATION
# ============================================================================

cat("\n=== CREATING COMBINED GFR COMPARISON PLOT ===\n\n")

# Create a combined plot showing all GFR categories together
mpag_gfr_long <- mpag_gfr_data %>%
    filter(!is.na(GFR_Category)) %>%
    mutate(GFR_Category = factor(GFR_Category, levels = gfr_order))

p_combined <- ggplot(mpag_gfr_long, aes(x = ACR_Group, y = `Mycophenolic.acid.O.acyl.glucuronide..C18.`, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    facet_wrap(~ GFR_Category, scales = "free_y") +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    labs(title = "MPAG Levels by GFR Category: 0R vs 2R+",
        subtitle = "Renally cleared metabolite",
         x = "ACR Group", y = "MPAG C18 Level") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

print(p_combined)
ggsave("Results2/Feedback_Analysis/MPAG_GFR_Stratified/MPAG_All_GFR_Categories_Combined.png", p_combined, 
       width = 12, height = 6, dpi = 300)

cat("\n=== DISPLAYING ALL GFR-STRATIFIED PLOTS ===\n")
cat("Total plots generated:", length(mpag_gfr_results), "\n\n")

# Display all individual GFR plots
for(gfr_cat in names(mpag_gfr_results)) {
    cat("Displaying:", gfr_cat, "\n")
    print(mpag_gfr_results[[gfr_cat]]$plot)
}

# ============================================================================
# CORRELATION ANALYSIS: GFR vs MPAG LEVELS
# ============================================================================

cat("\n\n=== CORRELATION: GFR vs MPA-GLUCURONIDE LEVELS ===\n\n")
cat("NOTE: Expecting STRONGER negative correlation than MPA (renal clearance!)\n\n")

# Overall correlation (all samples)
cat("--- Overall (All Samples) ---\n")
mpag_complete <- mpag_gfr_data %>% filter(!is.na(GFR) & !is.na(`Mycophenolic.acid.O.acyl.glucuronide..C18.`))
cor_overall <- cor.test(mpag_complete$GFR, mpag_complete$`Mycophenolic.acid.O.acyl.glucuronide..C18.`)
cat("Correlation: r =", round(cor_overall$estimate, 3), ", p =", round(cor_overall$p.value, 4), "\n\n")

# Correlation within 0R group
cat("--- 0R Group ---\n")
mpag_0r <- mpag_gfr_data %>% filter(ACR_Group == "0R" & !is.na(GFR) & !is.na(`Mycophenolic.acid.O.acyl.glucuronide..C18.`))
cor_0r <- cor.test(mpag_0r$GFR, mpag_0r$`Mycophenolic.acid.O.acyl.glucuronide..C18.`)
cat("Correlation: r =", round(cor_0r$estimate, 3), ", p =", round(cor_0r$p.value, 4), "\n\n")

# Correlation within 2R+ group
cat("--- 2R+ Group ---\n")
mpag_2r <- mpag_gfr_data %>% filter(ACR_Group == "2R+" & !is.na(GFR) & !is.na(`Mycophenolic.acid.O.acyl.glucuronide..C18.`))
cor_2r <- cor.test(mpag_2r$GFR, mpag_2r$`Mycophenolic.acid.O.acyl.glucuronide..C18.`)
cat("Correlation: r =", round(cor_2r$estimate, 3), ", p =", round(cor_2r$p.value, 4), "\n\n")

# Scatter plot with correlation annotations
p_scatter <- ggplot(mpag_gfr_data %>% filter(!is.na(GFR)), 
                    aes(x = GFR, y = `Mycophenolic.acid.O.acyl.glucuronide..C18.`, color = ACR_Group)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = c("0R" = "blue", "2R+" = "red")) +
    labs(title = "GFR vs MPAG Levels",
         subtitle = sprintf("Overall: r=%.3f, p=%.3f  |  0R: r=%.3f, p=%.3f  |  2R+: r=%.3f, p=%.3f",
                           cor_overall$estimate, cor_overall$p.value,
                           cor_0r$estimate, cor_0r$p.value,
                           cor_2r$estimate, cor_2r$p.value),
         x = "GFR (mL/min/1.73m²)", 
         y = "MPAG C18 Level") +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 9))

print(p_scatter)

ggsave("Results2/Feedback_Analysis/MPAG_GFR_Stratified/GFR_vs_MPAG_Scatter.png", p_scatter, 
       width = 10, height = 6, dpi = 300)

# Create separate scatter plots for better visibility
p_scatter_0r <- ggplot(mpag_0r, aes(x = GFR, y = `Mycophenolic.acid.O.acyl.glucuronide..C18.`)) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(title = "GFR vs MPAG: 0R Samples Only",
         subtitle = sprintf("r = %.3f, p = %.3f", cor_0r$estimate, cor_0r$p.value),
         x = "GFR (mL/min/1.73m²)", 
         y = "MPAG C18 Level") +
    theme_minimal()

p_scatter_2r <- ggplot(mpag_2r, aes(x = GFR, y = `Mycophenolic.acid.O.acyl.glucuronide..C18.`)) +
    geom_point(alpha = 0.6, color = "red") +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(title = "GFR vs MPAG: 2R+ Samples Only",
         subtitle = sprintf("r = %.3f, p = %.3f", cor_2r$estimate, cor_2r$p.value),
         x = "GFR (mL/min/1.73m²)", 
         y = "MPAG C18 Level") +
    theme_minimal()

print(p_scatter_0r)
print(p_scatter_2r)

ggsave("Results2/Feedback_Analysis/MPAG_GFR_Stratified/GFR_vs_MPAG_0R_Scatter.png", p_scatter_0r, 
       width = 8, height = 6, dpi = 300)
ggsave("Results2/Feedback_Analysis/MPAG_GFR_Stratified/GFR_vs_MPAG_2R_Scatter.png", p_scatter_2r, 
       width = 8, height = 6, dpi = 300)

# ============================================================================
# STATISTICAL INTERACTION TEST
# ============================================================================

cat("\n\n=== TESTING FOR GFR x ACR INTERACTION ===\n\n")

# Test if the effect of ACR differs by GFR (as continuous variable)
mpag_gfr_complete <- mpag_gfr_data %>% 
    filter(!is.na(GFR) & !is.na(ACR_Group) & !is.na(`Mycophenolic.acid.O.acyl.glucuronide..C18.`))

# Linear model with interaction
cat("--- MPAG C18 ---\n")
model_mpag <- lm(`Mycophenolic.acid.O.acyl.glucuronide..C18.` ~ GFR * ACR_Group, data = mpag_gfr_complete)
cat("Model coefficients:\n")
print(summary(model_mpag)$coefficients)

cat("\nInterpretation:\n")
interaction_p <- summary(model_mpag)$coefficients["GFR:ACR_Group2R+", "Pr(>|t|)"]
if(interaction_p < 0.05) {
    cat("✓ Significant interaction (p =", round(interaction_p, 4), ")\n")
    cat("  → The relationship between GFR and MPAG differs by ACR status\n")
    cat("  → Stratified analyses are ESSENTIAL\n\n")
} else {
    cat("✗ No significant interaction (p =", round(interaction_p, 4), ")\n")
    cat("  → GFR affects MPAG similarly in both 0R and 2R+ groups\n")
    cat("  → Main effects of GFR and ACR are independent\n\n")
}

# ============================================================================
# COMPARISON WITH MPA (PARENT DRUG)
# ============================================================================

cat("\n\n=== COMPARISON: MPAG vs MPA (Parent Drug) ===\n\n")

cat("Expected differences:\n")
cat("1. MPAG (this script):\n")
cat("   - Renally cleared → stronger negative correlation with GFR\n")
cat("   - Accumulates in kidney dysfunction\n")
cat("   - Higher levels in low GFR groups\n\n")

cat("2. MPA (script 12):\n")
cat("   - Less dependent on renal function\n")
cat("   - Weaker GFR correlation expected\n\n")

cat("To compare, run script 12 and compare:\n")
cat("   - Correlation coefficients (expect MPAG > MPA in magnitude)\n")
cat("   - GFR-stratified differences (expect larger MPAG effects)\n")
cat("   - Interaction effects\n\n")

# ============================================================================
# SUMMARY AND INTERPRETATION
# ============================================================================

cat("========================================================================\n")
cat("SUMMARY AND CLINICAL INTERPRETATION\n")
cat("========================================================================\n\n")

cat("=== Key Findings ===\n\n")

cat("1. GFR-STRATIFIED COMPARISONS (0R vs 2R+):\n")
if(length(mpag_gfr_results) > 0) {
    for(gfr_cat in names(mpag_gfr_results)) {
        cat("   ", gfr_cat, ": p =", round(mpag_gfr_results[[gfr_cat]]$p_value, 4), "\n")
    }
} else {
    cat("   No results generated (insufficient sample sizes)\n")
}
cat("\n")

cat("2. GFR-MPAG CORRELATION:\n")
cat("   Overall: r =", round(cor_overall$estimate, 3), ", p =", round(cor_overall$p.value, 4), "\n")
cat("   0R group: r =", round(cor_0r$estimate, 3), ", p =", round(cor_0r$p.value, 4), "\n")
cat("   2R+ group: r =", round(cor_2r$estimate, 3), ", p =", round(cor_2r$p.value, 4), "\n\n")

cat("3. GFR × ACR INTERACTION:\n")
cat("   p =", round(interaction_p, 4), "\n\n")

cat("=== Clinical Implications ===\n\n")

if(abs(cor_overall$estimate) > 0.3 & cor_overall$p.value < 0.05) {
    cat("⚠️  STRONG GFR DEPENDENCE DETECTED:\n")
    cat("   - MPAG levels are significantly affected by kidney function\n")
    cat("   - Renal clearance is the primary elimination route\n")
    cat("   - GFR-stratified analyses are CRITICAL for interpretation\n\n")
    
    cat("Clinical Recommendations:\n")
    cat("   1. Consider GFR when interpreting MPAG levels\n")
    cat("   2. Low GFR patients may accumulate metabolite\n")
    cat("   3. May need dose adjustments in kidney dysfunction\n")
    cat("   4. Monitor for toxicity in patients with reduced GFR\n\n")
} else {
    cat("✓ MODERATE GFR DEPENDENCE:\n")
    cat("   - Some relationship with kidney function but not dominant\n")
    cat("   - Other factors also influence MPAG levels\n\n")
}

if(interaction_p < 0.05) {
    cat("⚠️  SIGNIFICANT INTERACTION:\n")
    cat("   - GFR affects MPAG differently in 0R vs 2R+ groups\n")
    cat("   - Rejection status modifies the GFR-MPAG relationship\n")
    cat("   - MUST report stratified results\n\n")
}

cat("=== Files Generated ===\n")
cat("Results2/Feedback_Analysis/MPAG_GFR_Stratified/\n")
cat("  - MPAG_[Category]_0R_vs_2R.png (individual GFR categories)\n")
cat("  - MPAG_All_GFR_Categories_Combined.png\n")
cat("  - GFR_vs_MPAG_Scatter.png\n")
cat("  - GFR_vs_MPAG_0R_Scatter.png\n")
cat("  - GFR_vs_MPAG_2R_Scatter.png\n")
cat("  - MPAG_GFR_Stratified_Summary.csv\n\n")

cat("========================================================================\n")
cat("MPAG GFR-STRATIFIED ANALYSIS COMPLETE\n")
cat("========================================================================\n\n")
