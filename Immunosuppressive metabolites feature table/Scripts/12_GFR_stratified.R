# GFR-Stratified Analysis: 0R vs 2R+ MPA Levels

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

cat("\nGFR distribution:\n")
hist_data <- hist(clinical_data$GFR, breaks = 10, plot = FALSE)
print(hist_data$counts)
print(hist_data$breaks)

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

# Get MPA data and merge with GFR
mmf_gfr_data <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R")) %>%
    left_join(clinical_data %>% select(H, GFR, GFR_Category, Age, Sex, Race, BMI), by = "H")

cat("Total samples after merge:", nrow(mmf_gfr_data), "\n")
cat("Samples with GFR information:", sum(!is.na(mmf_gfr_data$GFR)), "\n\n")

# Check GFR category distribution in merged data
cat("GFR Category distribution in merged data:\n")
print(table(mmf_gfr_data$GFR_Category, useNA = "ifany"))

cat("\nGFR Category by ACR group:\n")
print(table(mmf_gfr_data$GFR_Category, mmf_gfr_data$ACR_Group, useNA = "ifany"))

# ============================================================================
# GFR-STRATIFIED ANALYSIS
# ============================================================================

cat("\n\n=== GFR-STRATIFIED ANALYSIS (0R vs 2R+) ===\n\n")

# Function to analyze MPA levels by GFR category
analyze_by_gfr <- function(data, gfr_category) {
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
    
    # Wilcoxon tests for C18
    cat("--- MPA (C18) ---\n")
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ ACR_Group, 
                         data = gfr_data, exact = FALSE)
    cat("p-value:", round(w_c18$p.value, 4), "\n")
    cat("Median 0R:", round(median(gfr_data$`Mycophenolate..C18.`[gfr_data$ACR_Group == "0R"], na.rm = TRUE), 3), "\n")
    cat("Median 2R+:", round(median(gfr_data$`Mycophenolate..C18.`[gfr_data$ACR_Group == "2R+"], na.rm = TRUE), 3), "\n\n")
    
    # Wilcoxon tests for HILIC
    cat("--- MPA (HILIC) ---\n")
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ ACR_Group, 
                           data = gfr_data, exact = FALSE)
    cat("p-value:", round(w_hilic$p.value, 4), "\n")
    cat("Median 0R:", round(median(gfr_data$`Mycophenolate..HILIC.`[gfr_data$ACR_Group == "0R"], na.rm = TRUE), 3), "\n")
    cat("Median 2R+:", round(median(gfr_data$`Mycophenolate..HILIC.`[gfr_data$ACR_Group == "2R+"], na.rm = TRUE), 3), "\n\n")
    
    # Create plot
    gfr_long <- gfr_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite",
                     values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(gfr_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
        scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r, ")"),
                                    "2R+" = paste0("2R+\n(n=", n_2r, ")"))) +
        labs(title = paste("MPA Levels:", gfr_category, "(0R vs 2R+)"),
             subtitle = paste0("C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "ACR Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(
        plot = p, 
        gfr_category = gfr_category,
        n_0r = n_0r,
        n_2r = n_2r,
        p_c18 = w_c18$p.value,
        p_hilic = w_hilic$p.value,
        median_0r_c18 = median(gfr_data$`Mycophenolate..C18.`[gfr_data$ACR_Group == "0R"], na.rm = TRUE),
        median_2r_c18 = median(gfr_data$`Mycophenolate..C18.`[gfr_data$ACR_Group == "2R+"], na.rm = TRUE),
        median_0r_hilic = median(gfr_data$`Mycophenolate..HILIC.`[gfr_data$ACR_Group == "0R"], na.rm = TRUE),
        median_2r_hilic = median(gfr_data$`Mycophenolate..HILIC.`[gfr_data$ACR_Group == "2R+"], na.rm = TRUE),
        gfr_min = min(gfr_data$GFR, na.rm = TRUE),
        gfr_max = max(gfr_data$GFR, na.rm = TRUE)
    ))
}

# Get unique GFR categories with data (in logical order)
gfr_order <- c("Normal (≥90)", "Mild Reduction (60-89)", 
               "Moderate Reduction (30-59)", "Severe (<30)")
gfr_categories <- gfr_order[gfr_order %in% unique(mmf_gfr_data$GFR_Category[!is.na(mmf_gfr_data$GFR_Category)])]

cat("GFR categories to analyze:", paste(gfr_categories, collapse = ", "), "\n")

# Create output directory
dir.create("Results2/Feedback_Analysis/GFR_Stratified", showWarnings = FALSE, recursive = TRUE)

# Run analysis for each GFR category
gfr_results <- list()
for(gfr_cat in gfr_categories) {
    result <- analyze_by_gfr(mmf_gfr_data, gfr_cat)
    if(!is.null(result)) {
        gfr_results[[gfr_cat]] <- result
    }
}

# Save all plots
cat("\n=== SAVING PLOTS ===\n")
for(gfr_cat in names(gfr_results)) {
    filename <- paste0("Results2/Feedback_Analysis/GFR_Stratified/MMF_", 
                      gsub("[^A-Za-z0-9]", "_", gfr_cat), "_0R_vs_2R.png")
    ggsave(filename, gfr_results[[gfr_cat]]$plot, width = 10, height = 6, dpi = 300)
    cat("Saved:", filename, "\n")
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n\n=== GFR-STRATIFIED SUMMARY ===\n\n")

if(length(gfr_results) > 0) {
    summary_df <- data.frame(
        GFR_Category = sapply(gfr_results, function(x) x$gfr_category),
        GFR_Range = paste0(sapply(gfr_results, function(x) round(x$gfr_min, 1)), 
                          " - ", 
                          sapply(gfr_results, function(x) round(x$gfr_max, 1))),
        N_0R = sapply(gfr_results, function(x) x$n_0r),
        N_2R_plus = sapply(gfr_results, function(x) x$n_2r),
        Median_0R_C18 = round(sapply(gfr_results, function(x) x$median_0r_c18), 0),
        Median_2R_C18 = round(sapply(gfr_results, function(x) x$median_2r_c18), 0),
        P_C18 = round(sapply(gfr_results, function(x) x$p_c18), 4),
        Median_0R_HILIC = round(sapply(gfr_results, function(x) x$median_0r_hilic), 0),
        Median_2R_HILIC = round(sapply(gfr_results, function(x) x$median_2r_hilic), 0),
        P_HILIC = round(sapply(gfr_results, function(x) x$p_hilic), 4)
    )
    
    # Reorder by GFR category
    summary_df$GFR_Category <- factor(summary_df$GFR_Category, levels = gfr_order)
    summary_df <- summary_df %>% arrange(GFR_Category)
    summary_df$GFR_Category <- as.character(summary_df$GFR_Category)
    
    print(summary_df)
    
    # Save summary table
    write_csv(summary_df, "Results2/Feedback_Analysis/GFR_Stratified/GFR_Stratified_Summary.csv")
    cat("\nSummary table saved to: Results2/Feedback_Analysis/GFR_Stratified/GFR_Stratified_Summary.csv\n")
}

# ============================================================================
# COMBINED VISUALIZATION
# ============================================================================

cat("\n=== CREATING COMBINED GFR COMPARISON PLOT ===\n\n")

# Create a combined plot showing all GFR categories together
mmf_gfr_long <- mmf_gfr_data %>%
    filter(!is.na(GFR_Category)) %>%
    mutate(GFR_Category = factor(GFR_Category, levels = gfr_order)) %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))

p_combined <- ggplot(mmf_gfr_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    facet_grid(Metabolite ~ GFR_Category, scales = "free_y") +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    labs(title = "MPA Levels by GFR Category: 0R vs 2R+",
         x = "ACR Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x = element_text(size = 8))

print(p_combined)
ggsave("Results2/Feedback_Analysis/GFR_Stratified/MMF_All_GFR_Categories_Combined.png", p_combined, 
       width = 14, height = 8, dpi = 300)

cat("\n=== DISPLAYING ALL GFR-STRATIFIED PLOTS ===\n")
cat("Total plots generated:", length(gfr_results), "\n\n")

# Display all individual GFR plots
for(gfr_cat in names(gfr_results)) {
    cat("Displaying:", gfr_cat, "\n")
    print(gfr_results[[gfr_cat]]$plot)
}

# ============================================================================
# CORRELATION ANALYSIS: GFR vs MPA LEVELS
# ============================================================================

cat("\n\n=== CORRELATION: GFR vs MPA LEVELS ===\n\n")

# Correlation within 0R group
cat("--- 0R Group ---\n")
gfr_0r <- mmf_gfr_data %>% filter(ACR_Group == "0R" & !is.na(GFR))
cor_0r_c18 <- cor.test(gfr_0r$GFR, gfr_0r$`Mycophenolate..C18.`)
cor_0r_hilic <- cor.test(gfr_0r$GFR, gfr_0r$`Mycophenolate..HILIC.`)
cat("C18 correlation: r =", round(cor_0r_c18$estimate, 3), ", p =", round(cor_0r_c18$p.value, 4), "\n")
cat("HILIC correlation: r =", round(cor_0r_hilic$estimate, 3), ", p =", round(cor_0r_hilic$p.value, 4), "\n\n")

# Correlation within 2R+ group
cat("--- 2R+ Group ---\n")
gfr_2r <- mmf_gfr_data %>% filter(ACR_Group == "2R+" & !is.na(GFR))
cor_2r_c18 <- cor.test(gfr_2r$GFR, gfr_2r$`Mycophenolate..C18.`)
cor_2r_hilic <- cor.test(gfr_2r$GFR, gfr_2r$`Mycophenolate..HILIC.`)
cat("C18 correlation: r =", round(cor_2r_c18$estimate, 3), ", p =", round(cor_2r_c18$p.value, 4), "\n")
cat("HILIC correlation: r =", round(cor_2r_hilic$estimate, 3), ", p =", round(cor_2r_hilic$p.value, 4), "\n\n")

# Scatter plots with correlation annotations
p_scatter_c18 <- ggplot(mmf_gfr_data %>% filter(!is.na(GFR)), 
                        aes(x = GFR, y = `Mycophenolate..C18.`, color = ACR_Group)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = c("0R" = "blue", "2R+" = "red")) +
    labs(title = "GFR vs MPA C18",
         subtitle = sprintf("0R: r=%.3f, p=%.3f  |  2R+: r=%.3f, p=%.3f",
                           cor_0r_c18$estimate, cor_0r_c18$p.value,
                           cor_2r_c18$estimate, cor_2r_c18$p.value),
         x = "GFR (mL/min/1.73m²)", 
         y = "MPA C18 Level") +
    theme_minimal()

p_scatter_hilic <- ggplot(mmf_gfr_data %>% filter(!is.na(GFR)), 
                          aes(x = GFR, y = `Mycophenolate..HILIC.`, color = ACR_Group)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = c("0R" = "blue", "2R+" = "red")) +
    labs(title = "GFR vs MPA HILIC",
         subtitle = sprintf("0R: r=%.3f, p=%.3f  |  2R+: r=%.3f, p=%.3f",
                           cor_0r_hilic$estimate, cor_0r_hilic$p.value,
                           cor_2r_hilic$estimate, cor_2r_hilic$p.value),
         x = "GFR (mL/min/1.73m²)", 
         y = "MPA HILIC Level") +
    theme_minimal()

print(p_scatter_c18)
print(p_scatter_hilic)

ggsave("Results2/Feedback_Analysis/GFR_Stratified/GFR_vs_MMF_C18_Scatter.png", p_scatter_c18, 
       width = 10, height = 6, dpi = 300)
ggsave("Results2/Feedback_Analysis/GFR_Stratified/GFR_vs_MMF_HILIC_Scatter.png", p_scatter_hilic, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# STATISTICAL INTERACTION TEST
# ============================================================================

cat("\n\n=== TESTING FOR GFR x ACR INTERACTION ===\n\n")

# Test if the effect of ACR differs by GFR (as continuous variable)
mmf_gfr_complete <- mmf_gfr_data %>% 
    filter(!is.na(GFR) & !is.na(ACR_Group))

# Linear model with interaction for C18
cat("--- MPA (C18) ---\n")
model_c18 <- lm(`Mycophenolate..C18.` ~ GFR * ACR_Group, data = mmf_gfr_complete)
cat("Interaction effect:\n")
print(summary(model_c18)$coefficients)

# Linear model with interaction for HILIC
cat("\n--- MPA HILIC ---\n")
model_hilic <- lm(`Mycophenolate..HILIC.` ~ GFR * ACR_Group, data = mmf_gfr_complete)
cat("Interaction effect:\n")
print(summary(model_hilic)$coefficients)

cat("\n=== GFR-STRATIFIED ANALYSIS COMPLETE ===\n")
cat("Results saved to: Results2/Feedback_Analysis/GFR_Stratified/\n\n")
