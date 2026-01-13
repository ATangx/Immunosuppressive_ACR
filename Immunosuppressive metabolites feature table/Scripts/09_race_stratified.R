# Race-Stratified Analysis: 0R vs 2R+ MPA Levels (Feedback Q3)

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# NOTE: Comparing 0R vs 2R+ only (excluding 1R for clearer comparison)

# ============================================================================
# LOAD CLINICAL DATA WITH RACE INFORMATION
# ============================================================================

cat("\n=== LOADING CLINICAL DATA ===\n")

# Load clinical data
clinical_data <- read_csv("OHT_Clinical.csv")

cat("Clinical data loaded:", nrow(clinical_data), "patients\n")
cat("Race distribution in clinical data:\n")
print(table(clinical_data$Race))

# ============================================================================
# MERGE METABOLITE DATA WITH RACE INFORMATION
# ============================================================================

cat("\n=== MERGING METABOLITE AND CLINICAL DATA ===\n\n")

# Get MPA data and merge with race
mmf_race_data <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R")) %>%
    left_join(clinical_data %>% select(H, Race, Age, Sex, BMI), by = "H")

cat("Total samples after merge:", nrow(mmf_race_data), "\n")
cat("Samples with race information:", sum(!is.na(mmf_race_data$Race)), "\n\n")

# Check race distribution in merged data
cat("Race distribution in merged data:\n")
print(table(mmf_race_data$Race, useNA = "ifany"))

cat("\nRace by ACR group:\n")
print(table(mmf_race_data$Race, mmf_race_data$ACR_Group, useNA = "ifany"))

# ============================================================================
# RACE-STRATIFIED ANALYSIS
# ============================================================================

cat("\n\n=== RACE-STRATIFIED ANALYSIS (0R vs 2R+) ===\n\n")

# Function to analyze MPA levels by race
analyze_by_race <- function(data, race_group) {
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("RACE:", race_group, "\n")
    cat(rep("=", 70), "\n", sep = "")
    
    # Filter for this race
    race_data <- data %>% filter(Race == race_group)
    
    n_total <- nrow(race_data)
    n_0r <- sum(race_data$ACR_Group == "0R")
    n_2r <- sum(race_data$ACR_Group == "2R+")
    
    cat("Total samples:", n_total, "(0R:", n_0r, ", 2R+:", n_2r, ")\n\n")
    
    # Check if we have enough samples
    if(n_0r < 3 | n_2r < 3) {
        cat("WARNING: Insufficient sample size for analysis (need ≥3 per group)\n")
        return(NULL)
    }
    
    # Wilcoxon tests for C18
    cat("--- MPA (C18) ---\n")
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ ACR_Group, 
                         data = race_data, exact = FALSE)
    cat("p-value:", round(w_c18$p.value, 4), "\n")
    cat("Median 0R:", round(median(race_data$`Mycophenolate..C18.`[race_data$ACR_Group == "0R"], na.rm = TRUE), 3), "\n")
    cat("Median 2R+:", round(median(race_data$`Mycophenolate..C18.`[race_data$ACR_Group == "2R+"], na.rm = TRUE), 3), "\n\n")
    
    # Wilcoxon tests for HILIC
    cat("--- MPA (HILIC) ---\n")
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ ACR_Group, 
                           data = race_data, exact = FALSE)
    cat("p-value:", round(w_hilic$p.value, 4), "\n")
    cat("Median 0R:", round(median(race_data$`Mycophenolate..HILIC.`[race_data$ACR_Group == "0R"], na.rm = TRUE), 3), "\n")
    cat("Median 2R+:", round(median(race_data$`Mycophenolate..HILIC.`[race_data$ACR_Group == "2R+"], na.rm = TRUE), 3), "\n\n")
    
    # Create plot
    race_long <- race_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite",
                     values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(race_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
        scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r, ")"),
                                    "2R+" = paste0("2R+\n(n=", n_2r, ")"))) +
        labs(title = paste("MPA Levels:", race_group, "patients (0R vs 2R+)"),
             subtitle = paste0("C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "ACR Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(
        plot = p, 
        race = race_group,
        n_0r = n_0r,
        n_2r = n_2r,
        p_c18 = w_c18$p.value,
        p_hilic = w_hilic$p.value
    ))
}

# Get unique races with sufficient data
races <- unique(mmf_race_data$Race[!is.na(mmf_race_data$Race)])
cat("Races to analyze:", paste(races, collapse = ", "), "\n")

# Create output directory
dir.create("Results2/Feedback_Analysis/Race_Stratified", showWarnings = FALSE, recursive = TRUE)

# Run analysis for each race (store results, don't print yet)
race_results <- list()
for(race in races) {
    result <- analyze_by_race(mmf_race_data, race)
    if(!is.null(result)) {
        race_results[[race]] <- result
    }
}

# Save all plots
cat("\n=== SAVING PLOTS ===\n")
for(race in names(race_results)) {
    filename <- paste0("Results2/Feedback_Analysis/Race_Stratified/MMF_", 
                      gsub(" ", "_", race), "_0R_vs_2R.png")
    ggsave(filename, race_results[[race]]$plot, width = 10, height = 6, dpi = 300)
    cat("Saved:", filename, "\n")
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n\n=== RACE-STRATIFIED SUMMARY ===\n\n")

if(length(race_results) > 0) {
    summary_df <- data.frame(
        Race = sapply(race_results, function(x) x$race),
        N_0R = sapply(race_results, function(x) x$n_0r),
        N_2R_plus = sapply(race_results, function(x) x$n_2r),
        P_C18 = round(sapply(race_results, function(x) x$p_c18), 4),
        P_HILIC = round(sapply(race_results, function(x) x$p_hilic), 4)
    )
    
    print(summary_df)
    
    # Save summary table
    write_csv(summary_df, "Results2/Feedback_Analysis/Race_Stratified/Race_Stratified_Summary.csv")
    cat("\nSummary table saved to: Results2/Feedback_Analysis/Race_Stratified/Race_Stratified_Summary.csv\n")
}

# ============================================================================
# COMBINED VISUALIZATION
# ============================================================================

cat("\n=== CREATING COMBINED RACE COMPARISON PLOT ===\n\n")

# Create a combined plot showing all races together
mmf_race_long <- mmf_race_data %>%
    filter(!is.na(Race)) %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))

p_combined <- ggplot(mmf_race_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    facet_grid(Metabolite ~ Race, scales = "free_y") +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    labs(title = "MPA Levels by Race: 0R vs 2R+",
         x = "ACR Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

print(p_combined)
ggsave("Results2/Feedback_Analysis/Race_Stratified/MMF_All_Races_Combined.png", p_combined, 
       width = 12, height = 8, dpi = 300)

cat("\n=== DISPLAYING ALL RACE-STRATIFIED PLOTS ===\n")
cat("Total plots generated:", length(race_results), "\n\n")

# Display all individual race plots
for(race in names(race_results)) {
    cat("Displaying:", race, "\n")
    print(race_results[[race]]$plot)
}

cat("\n✓ All", length(race_results), "race-stratified plots displayed.\n")
cat("✓ Use the ← → arrow buttons in the Plots pane to navigate.\n")
cat("✓ Files saved to: Results2/Feedback_Analysis/Race_Stratified/\n\n")

cat("All plots saved to: Results2/Feedback_Analysis/Race_Stratified/\n\n")
cat("Analysis complete!\n")
