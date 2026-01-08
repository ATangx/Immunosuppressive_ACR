# Simplified MMF Analysis with POD Stratification

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# ============================================================================
# PART 1: Overall MMF comparison (0R/1R vs 2R+)
# ============================================================================

cat("\n=== OVERALL MMF ANALYSIS (0R/1R vs 2R+) ===\n\n")

# Get MMF data
mmf_data <- patients_with_2R %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R/1R"))

# Run Wilcoxon tests
cat("Mycophenolate C18:\n")
wilcox_c18 <- wilcox.test(
    `Mycophenolate..C18.` ~ ACR_Group, 
    data = mmf_data, 
    exact = FALSE
)
print(wilcox_c18)

cat("\nMycophenolate HILIC:\n")
wilcox_hilic <- wilcox.test(
    `Mycophenolate..HILIC.` ~ ACR_Group, 
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

p_overall <- ggplot(mmf_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    scale_fill_manual(values = c("0R/1R" = "lightblue", "2R+" = "lightcoral")) +
    labs(title = "MMF Levels: 0R/1R vs 2R+",
         subtitle = paste0("C18: p=", round(wilcox_c18$p.value, 3), 
                          "; HILIC: p=", round(wilcox_hilic$p.value, 3)),
         x = "ACR Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "none")

print(p_overall)
ggsave("Results2/Feedback_Analysis/MMF_Overall.png", p_overall, 
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
    cat("Sample size:", nrow(stratum_data), 
        "(0R/1R:", sum(stratum_data$ACR_Group == "0R/1R"), 
        ", 2R+:", sum(stratum_data$ACR_Group == "2R+"), ")\n\n")
    
    # Wilcoxon tests
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ ACR_Group, 
                         data = stratum_data, exact = FALSE)
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ ACR_Group, 
                           data = stratum_data, exact = FALSE)
    
    cat("C18:   p =", round(w_c18$p.value, 3), "\n")
    cat("HILIC: p =", round(w_hilic$p.value, 3), "\n\n")
    
    # Create plot
    stratum_long <- stratum_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite", values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(stratum_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R/1R" = "lightblue", "2R+" = "lightcoral")) +
        labs(title = paste("MMF Levels:", stratum),
             subtitle = paste0("C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "ACR Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(plot = p, p_c18 = w_c18$p.value, p_hilic = w_hilic$p.value))
}

# Run for each stratum
dir.create("Results2/Feedback_Analysis", showWarnings = FALSE, recursive = TRUE)

early_results <- analyze_stratum(mmf_data, "Early (≤30 days)")
print(early_results$plot)
ggsave("Results2/Feedback_Analysis/MMF_Early.png", early_results$plot, 
       width = 10, height = 6, dpi = 300)

late_results <- analyze_stratum(mmf_data, "Late (>30 days)")
print(late_results$plot)
ggsave("Results2/Feedback_Analysis/MMF_Late.png", late_results$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n=== SUMMARY ===\n")
cat("Overall:     C18 p=", round(wilcox_c18$p.value, 3), 
    ", HILIC p=", round(wilcox_hilic$p.value, 3), "\n")
cat("Early ≤30d:  C18 p=", round(early_results$p_c18, 3), 
    ", HILIC p=", round(early_results$p_hilic, 3), "\n")
cat("Late >30d:   C18 p=", round(late_results$p_c18, 3), 
    ", HILIC p=", round(late_results$p_hilic, 3), "\n")
cat("\nPlots saved to: Results2/Feedback_Analysis/\n")
cat("- MMF_Overall.png\n")
cat("- MMF_Early.png\n")
cat("- MMF_Late.png\n\n")
