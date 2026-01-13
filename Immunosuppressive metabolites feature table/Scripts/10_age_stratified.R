# Age-Stratified Analysis: 0R vs 2R+ MPA Levels

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# NOTE: Comparing 0R vs 2R+ only (excluding 1R for clearer comparison)

# ============================================================================
# LOAD CLINICAL DATA WITH AGE INFORMATION
# ============================================================================

cat("\n=== LOADING CLINICAL DATA ===\n")

# Load clinical data
clinical_data <- read_csv("OHT_Clinical.csv")

cat("Clinical data loaded:", nrow(clinical_data), "patients\n")

# Check age distribution
cat("\nAge distribution:\n")
cat("Min:", min(clinical_data$Age, na.rm = TRUE), "\n")
cat("Max:", max(clinical_data$Age, na.rm = TRUE), "\n")
cat("Mean:", round(mean(clinical_data$Age, na.rm = TRUE), 1), "\n")
cat("Median:", round(median(clinical_data$Age, na.rm = TRUE), 1), "\n")
cat("Quartiles:\n")
print(quantile(clinical_data$Age, na.rm = TRUE))

# ============================================================================
# MERGE METABOLITE DATA WITH AGE INFORMATION
# ============================================================================

cat("\n=== MERGING METABOLITE AND CLINICAL DATA ===\n\n")

# Get MPA data and merge with age
mmf_age_data <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R")) %>%
    left_join(clinical_data %>% select(H, Age, Race, Sex, BMI), by = "H")

cat("Total samples after merge:", nrow(mmf_age_data), "\n")
cat("Samples with age information:", sum(!is.na(mmf_age_data$Age)), "\n\n")

# Create age groups
# Using clinically meaningful cutoffs: <40, 40-60, >60 years
mmf_age_data <- mmf_age_data %>%
    mutate(Age_Group = case_when(
        Age < 40 ~ "Young (<40 years)",
        Age >= 40 & Age < 60 ~ "Middle-aged (40-60 years)",
        Age >= 60 ~ "Older (≥60 years)",
        TRUE ~ NA_character_
    ))

# Alternative: binary classification at age 50
mmf_age_data <- mmf_age_data %>%
    mutate(Age_Binary = ifelse(Age < 50, "Younger (<50 years)", "Older (≥50 years)"))

cat("Age group distribution (3 groups):\n")
print(table(mmf_age_data$Age_Group, useNA = "ifany"))

cat("\nAge group distribution (binary):\n")
print(table(mmf_age_data$Age_Binary, useNA = "ifany"))

cat("\nAge group by ACR status (3 groups):\n")
print(table(mmf_age_data$Age_Group, mmf_age_data$ACR_Group, useNA = "ifany"))

cat("\nAge group by ACR status (binary):\n")
print(table(mmf_age_data$Age_Binary, mmf_age_data$ACR_Group, useNA = "ifany"))

# ============================================================================
# AGE-STRATIFIED ANALYSIS (3 GROUPS: <40, 40-60, ≥60)
# ============================================================================

cat("\n\n=== AGE-STRATIFIED ANALYSIS (0R vs 2R+) - 3 AGE GROUPS ===\n\n")

# Function to analyze MPA levels by age group
analyze_by_age <- function(data, age_group_name, age_var = "Age_Group") {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("AGE GROUP:", age_group_name, "\n")
    cat(rep("=", 70), "\n", sep = "")
    
    # Filter for this age group
    age_data <- data %>% filter(.data[[age_var]] == age_group_name)
    
    n_total <- nrow(age_data)
    n_0r <- sum(age_data$ACR_Group == "0R")
    n_2r <- sum(age_data$ACR_Group == "2R+")
    
    cat("Total samples:", n_total, "(0R:", n_0r, ", 2R+:", n_2r, ")\n")
    
    if(n_total > 0) {
        cat("Age range:", round(min(age_data$Age, na.rm = TRUE), 1), "-", 
            round(max(age_data$Age, na.rm = TRUE), 1), "years\n")
        cat("Mean age:", round(mean(age_data$Age, na.rm = TRUE), 1), "years\n\n")
    }
    
    # Check if we have enough samples
    if(n_0r < 3 | n_2r < 3) {
        cat("WARNING: Insufficient sample size for analysis (need ≥3 per group)\n")
        return(NULL)
    }
    
    # Wilcoxon tests for C18
    cat("--- MPA (C18) ---\n")
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ ACR_Group, 
                         data = age_data, exact = FALSE)
    cat("p-value:", round(w_c18$p.value, 4), "\n")
    cat("Median 0R:", round(median(age_data$`Mycophenolate..C18.`[age_data$ACR_Group == "0R"], na.rm = TRUE), 3), "\n")
    cat("Median 2R+:", round(median(age_data$`Mycophenolate..C18.`[age_data$ACR_Group == "2R+"], na.rm = TRUE), 3), "\n\n")
    
    # Wilcoxon tests for HILIC
    cat("--- MPA (HILIC) ---\n")
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ ACR_Group, 
                           data = age_data, exact = FALSE)
    cat("p-value:", round(w_hilic$p.value, 4), "\n")
    cat("Median 0R:", round(median(age_data$`Mycophenolate..HILIC.`[age_data$ACR_Group == "0R"], na.rm = TRUE), 3), "\n")
    cat("Median 2R+:", round(median(age_data$`Mycophenolate..HILIC.`[age_data$ACR_Group == "2R+"], na.rm = TRUE), 3), "\n\n")
    
    # Create plot
    age_long <- age_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite",
                     values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(age_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
        scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r, ")"),
                                    "2R+" = paste0("2R+\n(n=", n_2r, ")"))) +
        labs(title = paste("MPA Levels:", age_group_name, "(0R vs 2R+)"),
             subtitle = paste0("C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "ACR Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(
        plot = p, 
        age_group = age_group_name,
        n_0r = n_0r,
        n_2r = n_2r,
        p_c18 = w_c18$p.value,
        p_hilic = w_hilic$p.value
    ))
}

# Create output directory
dir.create("Results2/Feedback_Analysis/Age_Stratified", showWarnings = FALSE, recursive = TRUE)

# Get unique age groups (3-group stratification) in logical order
age_groups_3 <- c("Young (<40 years)", "Middle-aged (40-60 years)", "Older (≥60 years)")
# Keep only groups that exist in the data
age_groups_3 <- age_groups_3[age_groups_3 %in% unique(mmf_age_data$Age_Group)]
cat("Age groups (3-level) to analyze:", paste(age_groups_3, collapse = ", "), "\n")

# Run analysis for each age group (store results, don't print yet)
age_results_3 <- list()
for(age_grp in age_groups_3) {
    result <- analyze_by_age(mmf_age_data, age_grp, age_var = "Age_Group")
    if(!is.null(result)) {
        age_results_3[[age_grp]] <- result
    }
}

# Save all plots (3-group)
cat("\n=== SAVING PLOTS (3 AGE GROUPS) ===\n")
for(age_grp in names(age_results_3)) {
    filename <- paste0("Results2/Feedback_Analysis/Age_Stratified/MMF_", 
                      gsub("[^A-Za-z0-9]", "_", age_grp), "_0R_vs_2R.png")
    ggsave(filename, age_results_3[[age_grp]]$plot, width = 10, height = 6, dpi = 300)
    cat("Saved:", filename, "\n")
}

# ============================================================================
# AGE-STRATIFIED ANALYSIS (BINARY: <50 vs ≥50)
# ============================================================================

cat("\n\n=== AGE-STRATIFIED ANALYSIS (0R vs 2R+) - BINARY AGE GROUPS ===\n\n")

# Get unique binary age groups in logical order
age_groups_binary <- c("Younger (<50 years)", "Older (≥50 years)")
# Keep only groups that exist in the data
age_groups_binary <- age_groups_binary[age_groups_binary %in% unique(mmf_age_data$Age_Binary)]
cat("Age groups (binary) to analyze:", paste(age_groups_binary, collapse = ", "), "\n")

# Run analysis for binary age groups
age_results_binary <- list()
for(age_grp in age_groups_binary) {
    result <- analyze_by_age(mmf_age_data, age_grp, age_var = "Age_Binary")
    if(!is.null(result)) {
        age_results_binary[[age_grp]] <- result
    }
}

# Save all plots (binary)
cat("\n=== SAVING PLOTS (BINARY AGE GROUPS) ===\n")
for(age_grp in names(age_results_binary)) {
    filename <- paste0("Results2/Feedback_Analysis/Age_Stratified/MMF_", 
                      gsub("[^A-Za-z0-9]", "_", age_grp), "_0R_vs_2R.png")
    ggsave(filename, age_results_binary[[age_grp]]$plot, width = 10, height = 6, dpi = 300)
    cat("Saved:", filename, "\n")
}

# ============================================================================
# SUMMARY TABLES
# ============================================================================

cat("\n\n=== AGE-STRATIFIED SUMMARY (3 GROUPS) ===\n\n")

if(length(age_results_3) > 0) {
    summary_df_3 <- data.frame(
        Age_Group = sapply(age_results_3, function(x) x$age_group),
        N_0R = sapply(age_results_3, function(x) x$n_0r),
        N_2R_plus = sapply(age_results_3, function(x) x$n_2r),
        P_C18 = round(sapply(age_results_3, function(x) x$p_c18), 4),
        P_HILIC = round(sapply(age_results_3, function(x) x$p_hilic), 4)
    )
    
    print(summary_df_3)
    
    # Save summary table
    write_csv(summary_df_3, "Results2/Feedback_Analysis/Age_Stratified/Age_Stratified_3groups_Summary.csv")
    cat("\nSummary table saved to: Results2/Feedback_Analysis/Age_Stratified/Age_Stratified_3groups_Summary.csv\n")
}

cat("\n=== AGE-STRATIFIED SUMMARY (BINARY) ===\n\n")

if(length(age_results_binary) > 0) {
    summary_df_binary <- data.frame(
        Age_Group = sapply(age_results_binary, function(x) x$age_group),
        N_0R = sapply(age_results_binary, function(x) x$n_0r),
        N_2R_plus = sapply(age_results_binary, function(x) x$n_2r),
        P_C18 = round(sapply(age_results_binary, function(x) x$p_c18), 4),
        P_HILIC = round(sapply(age_results_binary, function(x) x$p_hilic), 4)
    )
    
    print(summary_df_binary)
    
    # Save summary table
    write_csv(summary_df_binary, "Results2/Feedback_Analysis/Age_Stratified/Age_Stratified_Binary_Summary.csv")
    cat("\nSummary table saved to: Results2/Feedback_Analysis/Age_Stratified/Age_Stratified_Binary_Summary.csv\n")
}

# ============================================================================
# COMBINED VISUALIZATIONS
# ============================================================================

cat("\n=== CREATING COMBINED AGE COMPARISON PLOTS ===\n\n")

# Create a combined plot showing all age groups together (3-group)
mmf_age_long_3 <- mmf_age_data %>%
    filter(!is.na(Age_Group)) %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite),
           Age_Group = factor(Age_Group, 
                             levels = c("Young (<40 years)", 
                                       "Middle-aged (40-60 years)", 
                                       "Older (≥60 years)")))

p_combined_3 <- ggplot(mmf_age_long_3, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    facet_grid(Metabolite ~ Age_Group, scales = "free_y") +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    labs(title = "MPA Levels by Age Group (3 groups): 0R vs 2R+",
         x = "ACR Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Results2/Feedback_Analysis/Age_Stratified/MMF_All_Ages_3groups_Combined.png", p_combined_3, 
       width = 12, height = 8, dpi = 300)

# Create a combined plot for binary age groups
mmf_age_long_binary <- mmf_age_data %>%
    filter(!is.na(Age_Binary)) %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))

p_combined_binary <- ggplot(mmf_age_long_binary, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    facet_grid(Metabolite ~ Age_Binary, scales = "free_y") +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    labs(title = "MPA Levels by Age Group (Binary): 0R vs 2R+",
         x = "ACR Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Results2/Feedback_Analysis/Age_Stratified/MMF_All_Ages_Binary_Combined.png", p_combined_binary, 
       width = 10, height = 8, dpi = 300)

# ============================================================================
# DISPLAY ALL PLOTS
# ============================================================================

cat("\n=== DISPLAYING ALL AGE-STRATIFIED PLOTS ===\n")

# Display 3-group plots
cat("\n3-GROUP STRATIFICATION:\n")
cat("Total plots generated:", length(age_results_3), "\n\n")
for(age_grp in names(age_results_3)) {
    cat("Displaying:", age_grp, "\n")
    print(age_results_3[[age_grp]]$plot)
}

cat("\n✓ All", length(age_results_3), "age-stratified plots (3 groups) displayed.\n\n")

# Display binary plots
cat("BINARY STRATIFICATION:\n")
cat("Total plots generated:", length(age_results_binary), "\n\n")
for(age_grp in names(age_results_binary)) {
    cat("Displaying:", age_grp, "\n")
    print(age_results_binary[[age_grp]]$plot)
}

cat("\n✓ All", length(age_results_binary), "age-stratified plots (binary) displayed.\n\n")

# Display combined plots
cat("COMBINED PLOTS:\n")
print(p_combined_3)
print(p_combined_binary)

cat("\n✓ Use the ← → arrow buttons in the Plots pane to navigate.\n")
cat("✓ Files saved to: Results2/Feedback_Analysis/Age_Stratified/\n\n")

cat("Analysis complete!\n")
