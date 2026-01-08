# Sex-Stratified Analysis: 0R vs 2R+ Mycophenolate Levels

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# NOTE: Comparing 0R vs 2R+ only (excluding 1R for clearer comparison)

# ============================================================================
# LOAD CLINICAL DATA WITH SEX INFORMATION
# ============================================================================

cat("\n=== LOADING CLINICAL DATA ===\n")

# Load clinical data
clinical_data <- read_csv("OHT_Clinical.csv")

cat("Clinical data loaded:", nrow(clinical_data), "patients\n")
cat("Sex distribution in clinical data:\n")
print(table(clinical_data$Sex))

# ============================================================================
# MERGE METABOLITE DATA WITH SEX INFORMATION
# ============================================================================

cat("\n=== MERGING METABOLITE AND CLINICAL DATA ===\n\n")

# Get MMF data and merge with sex
mmf_sex_data <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R")) %>%
    left_join(clinical_data %>% select(H, Sex, Age, Race, BMI), by = "H")

cat("Total samples after merge:", nrow(mmf_sex_data), "\n")
cat("Samples with sex information:", sum(!is.na(mmf_sex_data$Sex)), "\n\n")

# Check sex distribution in merged data
cat("Sex distribution in merged data:\n")
print(table(mmf_sex_data$Sex, useNA = "ifany"))

cat("\nSex by ACR group:\n")
print(table(mmf_sex_data$Sex, mmf_sex_data$ACR_Group, useNA = "ifany"))

# ============================================================================
# SEX-STRATIFIED ANALYSIS
# ============================================================================

cat("\n\n=== SEX-STRATIFIED ANALYSIS (0R vs 2R+) ===\n\n")

# Function to analyze MMF levels by sex
analyze_by_sex <- function(data, sex_group) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("SEX:", sex_group, "\n")
    cat(rep("=", 70), "\n", sep = "")
    
    # Filter for this sex
    sex_data <- data %>% filter(Sex == sex_group)
    
    n_total <- nrow(sex_data)
    n_0r <- sum(sex_data$ACR_Group == "0R")
    n_2r <- sum(sex_data$ACR_Group == "2R+")
    
    cat("Total samples:", n_total, "(0R:", n_0r, ", 2R+:", n_2r, ")\n\n")
    
    # Check if we have enough samples
    if(n_0r < 3 | n_2r < 3) {
        cat("WARNING: Insufficient sample size for analysis (need ≥3 per group)\n")
        return(NULL)
    }
    
    # Wilcoxon tests for C18
    cat("--- Mycophenolate C18 ---\n")
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ ACR_Group, 
                         data = sex_data, exact = FALSE)
    cat("p-value:", round(w_c18$p.value, 4), "\n")
    cat("Median 0R:", round(median(sex_data$`Mycophenolate..C18.`[sex_data$ACR_Group == "0R"], na.rm = TRUE), 3), "\n")
    cat("Median 2R+:", round(median(sex_data$`Mycophenolate..C18.`[sex_data$ACR_Group == "2R+"], na.rm = TRUE), 3), "\n\n")
    
    # Wilcoxon tests for HILIC
    cat("--- Mycophenolate HILIC ---\n")
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ ACR_Group, 
                           data = sex_data, exact = FALSE)
    cat("p-value:", round(w_hilic$p.value, 4), "\n")
    cat("Median 0R:", round(median(sex_data$`Mycophenolate..HILIC.`[sex_data$ACR_Group == "0R"], na.rm = TRUE), 3), "\n")
    cat("Median 2R+:", round(median(sex_data$`Mycophenolate..HILIC.`[sex_data$ACR_Group == "2R+"], na.rm = TRUE), 3), "\n\n")
    
    # Create plot
    sex_long <- sex_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite",
                     values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(sex_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
        scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r, ")"),
                                    "2R+" = paste0("2R+\n(n=", n_2r, ")"))) +
        labs(title = paste("MMF Levels:", sex_group, "patients (0R vs 2R+)"),
             subtitle = paste0("C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "ACR Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(
        plot = p, 
        sex = sex_group,
        n_0r = n_0r,
        n_2r = n_2r,
        p_c18 = w_c18$p.value,
        p_hilic = w_hilic$p.value,
        median_0r_c18 = median(sex_data$`Mycophenolate..C18.`[sex_data$ACR_Group == "0R"], na.rm = TRUE),
        median_2r_c18 = median(sex_data$`Mycophenolate..C18.`[sex_data$ACR_Group == "2R+"], na.rm = TRUE),
        median_0r_hilic = median(sex_data$`Mycophenolate..HILIC.`[sex_data$ACR_Group == "0R"], na.rm = TRUE),
        median_2r_hilic = median(sex_data$`Mycophenolate..HILIC.`[sex_data$ACR_Group == "2R+"], na.rm = TRUE)
    ))
}

# Get unique sex groups with sufficient data
sex_groups <- unique(mmf_sex_data$Sex[!is.na(mmf_sex_data$Sex)])
cat("Sex groups to analyze:", paste(sex_groups, collapse = ", "), "\n")

# Create output directory
dir.create("Results2/Feedback_Analysis/Sex_Stratified", showWarnings = FALSE, recursive = TRUE)

# Run analysis for each sex group
sex_results <- list()
for(sex_group in sex_groups) {
    result <- analyze_by_sex(mmf_sex_data, sex_group)
    if(!is.null(result)) {
        sex_results[[sex_group]] <- result
    }
}

# Save all plots
cat("\n=== SAVING PLOTS ===\n")
for(sex_group in names(sex_results)) {
    filename <- paste0("Results2/Feedback_Analysis/Sex_Stratified/MMF_", 
                      gsub(" ", "_", sex_group), "_0R_vs_2R.png")
    ggsave(filename, sex_results[[sex_group]]$plot, width = 10, height = 6, dpi = 300)
    cat("Saved:", filename, "\n")
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n\n=== SEX-STRATIFIED SUMMARY ===\n\n")

if(length(sex_results) > 0) {
    summary_df <- data.frame(
        Sex = sapply(sex_results, function(x) x$sex),
        N_0R = sapply(sex_results, function(x) x$n_0r),
        N_2R_plus = sapply(sex_results, function(x) x$n_2r),
        Median_0R_C18 = round(sapply(sex_results, function(x) x$median_0r_c18), 3),
        Median_2R_C18 = round(sapply(sex_results, function(x) x$median_2r_c18), 3),
        P_C18 = round(sapply(sex_results, function(x) x$p_c18), 4),
        Median_0R_HILIC = round(sapply(sex_results, function(x) x$median_0r_hilic), 3),
        Median_2R_HILIC = round(sapply(sex_results, function(x) x$median_2r_hilic), 3),
        P_HILIC = round(sapply(sex_results, function(x) x$p_hilic), 4)
    )
    
    print(summary_df)
    
    # Save summary table
    write_csv(summary_df, "Results2/Feedback_Analysis/Sex_Stratified/Sex_Stratified_Summary.csv")
    cat("\nSummary table saved to: Results2/Feedback_Analysis/Sex_Stratified/Sex_Stratified_Summary.csv\n")
}

# ============================================================================
# COMBINED VISUALIZATION
# ============================================================================

cat("\n=== CREATING COMBINED SEX COMPARISON PLOT ===\n\n")

# Create a combined plot showing all sex groups together
mmf_sex_long <- mmf_sex_data %>%
    filter(!is.na(Sex)) %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))

p_combined <- ggplot(mmf_sex_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    facet_grid(Metabolite ~ Sex, scales = "free_y") +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    labs(title = "MMF Levels by Sex: 0R vs 2R+",
         x = "ACR Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

print(p_combined)
ggsave("Results2/Feedback_Analysis/Sex_Stratified/MMF_All_Sex_Groups_Combined.png", p_combined, 
       width = 10, height = 8, dpi = 300)

cat("\n=== DISPLAYING ALL SEX-STRATIFIED PLOTS ===\n")
cat("Total plots generated:", length(sex_results), "\n\n")

# Display all individual sex plots
for(sex_group in names(sex_results)) {
    cat("Displaying:", sex_group, "\n")
    print(sex_results[[sex_group]]$plot)
}

# ============================================================================
# STATISTICAL INTERACTION TEST
# ============================================================================

cat("\n\n=== TESTING FOR SEX x ACR INTERACTION ===\n\n")

# Test if the effect of ACR differs by sex
mmf_sex_complete <- mmf_sex_data %>% 
    filter(!is.na(Sex) & !is.na(ACR_Group))

# Linear model with interaction for C18
cat("--- Mycophenolate C18 ---\n")
model_c18 <- lm(`Mycophenolate..C18.` ~ Sex * ACR_Group, data = mmf_sex_complete)
cat("Interaction effect:\n")
print(summary(model_c18)$coefficients)

# Linear model with interaction for HILIC
cat("\n--- Mycophenolate HILIC ---\n")
model_hilic <- lm(`Mycophenolate..HILIC.` ~ Sex * ACR_Group, data = mmf_sex_complete)
cat("Interaction effect:\n")
print(summary(model_hilic)$coefficients)

cat("\n=== SEX-STRATIFIED ANALYSIS COMPLETE ===\n")
cat("Results saved to: Results2/Feedback_Analysis/Sex_Stratified/\n\n")
