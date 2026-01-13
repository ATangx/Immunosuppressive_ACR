# MPA Analysis with POD Stratification (Feedback Q1 and Q2)
#! This is only comparing samples within patients who have at least one 2R episode
    # Question in mind: Do MPA levels differ between rejection states, and does this vary by POD?

# Load data (run 00_source first)
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# NOTE: Comparing 0R vs 2R+ only (excluding 1R for clearer comparison)
# This focuses on normal biopsies vs clinically significant rejection

# ============================================================================
# PART 1: Overall MPA comparison (0R vs 2R+)
# ============================================================================

cat("\n=== OVERALL MPA ANALYSIS (0R vs 2R+) ===\n")
cat("NOTE: Excluding 1R samples for clearer comparison\n\n")

# Get MPA data - filter to keep only 0R and 2R+
mpa_data <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R"))

# Run Wilcoxon tests
cat("MPA (C18):\n")
wilcox_c18 <- wilcox.test(
    `Mycophenolate..C18.` ~ ACR_Group, 
    data = mpa_data, 
    exact = FALSE
)
print(wilcox_c18)

cat("\nMPA (HILIC):\n")
wilcox_hilic <- wilcox.test(
    `Mycophenolate..HILIC.` ~ ACR_Group, 
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
n_0r_overall <- sum(mpa_data$ACR_Group == "0R")
n_2r_overall <- sum(mpa_data$ACR_Group == "2R+")

p_overall <- ggplot(mpa_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r_overall, ")"),
                                "2R+" = paste0("2R+\n(n=", n_2r_overall, ")"))) +
    labs(title = "MPA Levels: 0R vs 2R+",
         subtitle = paste0("C18: p=", round(wilcox_c18$p.value, 3), 
                          "; HILIC: p=", round(wilcox_hilic$p.value, 3)),
         x = "ACR Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "none")

print(p_overall)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Overall.png", p_overall, 
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
    n_0r <- sum(stratum_data$ACR_Group == "0R")
    n_2r <- sum(stratum_data$ACR_Group == "2R+")
    
    cat("Sample size:", nrow(stratum_data), 
        "(0R:", n_0r, ", 2R+:", n_2r, ")\n\n")
    
    # Wilcoxon tests
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ ACR_Group, 
                         data = stratum_data, exact = FALSE)
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ ACR_Group, 
                           data = stratum_data, exact = FALSE)
    
    cat("C18:   p =", round(w_c18$p.value, 3), "\n")
    cat("HILIC: p =", round(w_hilic$p.value, 3), "\n\n")
    
    # Create plot with sample sizes in x-axis labels
    stratum_long <- stratum_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite", values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(stratum_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
        scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r, ")"),
                                    "2R+" = paste0("2R+\n(n=", n_2r, ")"))) +
        labs(title = paste("MPA Levels:", stratum),
             subtitle = paste0("C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "ACR Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(plot = p, p_c18 = w_c18$p.value, p_hilic = w_hilic$p.value))
}

# Run for each stratum
dir.create("Results2/Feedback_Analysis/POD_Stratified_new", showWarnings = FALSE, recursive = TRUE)

early_results <- analyze_stratum(mpa_data, "Early (≤30 days)")
print(early_results$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Early.png", early_results$plot, 
       width = 10, height = 6, dpi = 300)

late_results <- analyze_stratum(mpa_data, "Late (>30 days)")
print(late_results$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Late.png", late_results$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 3: POD-Stratified Analysis EXCLUDING POD<10 (Early ≤30 days vs Late >30 days)
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA ANALYSIS (EXCLUDING POD<10) ===\n")
cat("NOTE: Removing very early post-transplant samples (POD<10)\n\n")

# Filter out POD<10 samples
mmf_data_pod10plus <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    filter(POD >= 10) %>%  # Exclude POD < 10
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R"))

cat("Total samples with POD >= 10:", nrow(mmf_data_pod10plus), "\n")
cat("0R samples:", sum(mmf_data_pod10plus$ACR_Group == "0R"), "\n")
cat("2R+ samples:", sum(mmf_data_pod10plus$ACR_Group == "2R+"), ")\n\n")

# Add POD stratum
mmf_data_pod10plus <- mmf_data_pod10plus %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (>30 days)"))

# Run for each stratum (using the same analyze_stratum function)
early_results_pod10 <- analyze_stratum(mmf_data_pod10plus, "Early (10-30 days)")
print(early_results_pod10$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Early_POD10plus.png", early_results_pod10$plot, 
       width = 10, height = 6, dpi = 300)

late_results_pod10 <- analyze_stratum(mmf_data_pod10plus, "Late (>30 days)")
print(late_results_pod10$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Late_POD10plus.png", late_results_pod10$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 4: POD-Stratified Analysis EXCLUDING POD<10 AND POD>45
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA ANALYSIS (POD 10-45 ONLY) ===\n")
cat("NOTE: Removing very early (POD<10) and late (POD>45) samples\n")
cat("      to focus on the acute post-transplant window\n\n")

# Filter to POD 10-45 only
mmf_data_pod10_45 <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    filter(POD >= 10 & POD <= 45) %>%  # Keep only POD 10-45
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R"))

cat("Total samples with POD 10-45:", nrow(mmf_data_pod10_45), "\n")
cat("0R samples:", sum(mmf_data_pod10_45$ACR_Group == "0R"), "\n")
cat("2R+ samples:", sum(mmf_data_pod10_45$ACR_Group == "2R+"), "\n\n")

# Add POD stratum
mmf_data_pod10_45 <- mmf_data_pod10_45 %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (31-45 days)"))

# Run for each stratum
early_results_pod10_45 <- analyze_stratum(mmf_data_pod10_45, "Early (10-30 days)")
print(early_results_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Early_POD10_45.png", early_results_pod10_45$plot, 
       width = 10, height = 6, dpi = 300)

late_results_pod10_45 <- analyze_stratum(mmf_data_pod10_45, "Late (31-45 days)")
print(late_results_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Late_POD10_45.png", late_results_pod10_45$plot, 
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
cat("Early 10-30d: C18 p=", round(early_results_pod10$p_c18, 3), 
    ", HILIC p=", round(early_results_pod10$p_hilic, 3), "\n")
cat("Late >30d:    C18 p=", round(late_results_pod10$p_c18, 3), 
    ", HILIC p=", round(late_results_pod10$p_hilic, 3), "\n")
cat("Early 10-45d: C18 p=", round(early_results_pod10_45$p_c18, 3), 
    ", HILIC p=", round(early_results_pod10_45$p_hilic, 3), "\n")
cat("Late 31-45d:  C18 p=", round(late_results_pod10_45$p_c18, 3), 
    ", HILIC p=", round(late_results_pod10_45$p_hilic, 3), "\n")
cat("\nPlots saved to: Results2/Feedback_Analysis/POD_Stratified_new/\n")
cat("- MPA_Overall.png\n")
cat("- MPA_Early.png\n")
cat("- MPA_Late.png\n")
cat("- MPA_Early_POD10plus.png\n")
cat("- MPA_Late_POD10plus.png\n")
cat("- MPA_Early_POD10_45.png\n")
cat("- MPA_Late_POD10_45.png\n\n")

# ============================================================================
# SUPPLEMENTARY ANALYSIS: 0R/1R vs 2R+ (Including 1R samples)
# ============================================================================

cat("\n\n")
cat("================================================================================\n")
cat("SUPPLEMENTARY ANALYSIS: MPA LEVELS 0R/1R vs 2R+\n")
cat("NOTE: This analysis INCLUDES 1R samples (grouped with 0R)\n")
cat("================================================================================\n\n")

# ============================================================================
# PART 1: Overall MPA comparison (0R/1R vs 2R+)
# ============================================================================

cat("\n=== OVERALL MPA ANALYSIS (0R/1R vs 2R+) ===\n")
cat("NOTE: Grouping 0R and 1R together vs 2R+\n\n")

# Get MPA data - include 0R, 1R, and 2R+
mmf_data_01r <- patients_with_2R %>%
    filter(ACR %in% c("0R", "1R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R/1R"))

# Run Wilcoxon tests
cat("MPA (C18):\n")
wilcox_c18_01r <- wilcox.test(
    `Mycophenolate..C18.` ~ ACR_Group, 
    data = mmf_data_01r, 
    exact = FALSE
)
print(wilcox_c18_01r)

cat("\nMPA (HILIC):\n")
wilcox_hilic_01r <- wilcox.test(
    `Mycophenolate..HILIC.` ~ ACR_Group, 
    data = mmf_data_01r, 
    exact = FALSE
)
print(wilcox_hilic_01r)

# Plot overall comparison
mmf_long_01r <- mmf_data_01r %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                 names_to = "Metabolite",
                 values_to = "Level") %>%
    mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))

# Get sample sizes for x-axis labels
n_01r_overall <- sum(mmf_data_01r$ACR_Group == "0R/1R")
n_2r_overall_01r <- sum(mmf_data_01r$ACR_Group == "2R+")

p_overall_01r <- ggplot(mmf_long_01r, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    scale_fill_manual(values = c("0R/1R" = "lightgreen", "2R+" = "lightcoral")) +
    scale_x_discrete(labels = c("0R/1R" = paste0("0R/1R\n(n=", n_01r_overall, ")"),
                                "2R+" = paste0("2R+\n(n=", n_2r_overall_01r, ")"))) +
    labs(title = "MPA Levels: 0R/1R vs 2R+",
         subtitle = paste0("C18: p=", round(wilcox_c18_01r$p.value, 3), 
                          "; HILIC: p=", round(wilcox_hilic_01r$p.value, 3)),
         x = "ACR Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "none")

print(p_overall_01r)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Overall_0R1R_vs_2R.png", p_overall_01r, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 2: POD-Stratified Analysis (Early ≤30 days vs Late >30 days)
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA ANALYSIS (0R/1R vs 2R+) ===\n\n")

# Add POD stratum
mmf_data_01r <- mmf_data_01r %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (≤30 days)", "Late (>30 days)"))

# Function to run tests for one stratum (0R/1R version)
analyze_stratum_01r <- function(data, stratum) {
    cat("--- ", stratum, " ---\n")
    
    stratum_data <- data %>% filter(POD_Stratum == stratum)
    n_01r <- sum(stratum_data$ACR_Group == "0R/1R")
    n_2r <- sum(stratum_data$ACR_Group == "2R+")
    
    cat("Sample size:", nrow(stratum_data), 
        "(0R/1R:", n_01r, ", 2R+:", n_2r, ")\n\n")
    
    # Wilcoxon tests
    w_c18 <- wilcox.test(`Mycophenolate..C18.` ~ ACR_Group, 
                         data = stratum_data, exact = FALSE)
    w_hilic <- wilcox.test(`Mycophenolate..HILIC.` ~ ACR_Group, 
                           data = stratum_data, exact = FALSE)
    
    cat("C18:   p =", round(w_c18$p.value, 3), "\n")
    cat("HILIC: p =", round(w_hilic$p.value, 3), "\n\n")
    
    # Create plot with sample sizes in x-axis labels
    stratum_long <- stratum_data %>%
        pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`),
                     names_to = "Metabolite", values_to = "Level") %>%
        mutate(Metabolite = gsub("\\.\\.", " ", Metabolite))
    
    p <- ggplot(stratum_long, aes(x = ACR_Group, y = Level, fill = ACR_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        facet_wrap(~ Metabolite, scales = "free_y") +
        scale_fill_manual(values = c("0R/1R" = "lightgreen", "2R+" = "lightcoral")) +
        scale_x_discrete(labels = c("0R/1R" = paste0("0R/1R\n(n=", n_01r, ")"),
                                    "2R+" = paste0("2R+\n(n=", n_2r, ")"))) +
        labs(title = paste("MPA Levels:", stratum, "(0R/1R vs 2R+)"),
             subtitle = paste0("C18: p=", round(w_c18$p.value, 3), 
                              "; HILIC: p=", round(w_hilic$p.value, 3)),
             x = "ACR Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(plot = p, p_c18 = w_c18$p.value, p_hilic = w_hilic$p.value))
}

# Run for each stratum
early_results_01r <- analyze_stratum_01r(mmf_data_01r, "Early (≤30 days)")
print(early_results_01r$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Early_0R1R_vs_2R.png", early_results_01r$plot, 
       width = 10, height = 6, dpi = 300)

late_results_01r <- analyze_stratum_01r(mmf_data_01r, "Late (>30 days)")
print(late_results_01r$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Late_0R1R_vs_2R.png", late_results_01r$plot, 
       width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 3: POD-Stratified Analysis EXCLUDING POD<10
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA ANALYSIS (0R/1R vs 2R+, EXCLUDING POD<10) ===\n")
cat("NOTE: Removing very early post-transplant samples (POD<10)\n\n")

# Filter out POD<10 samples
mmf_data_01r_pod10plus <- patients_with_2R %>%
    filter(ACR %in% c("0R", "1R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    filter(POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R/1R"))

cat("Total samples with POD >= 10:", nrow(mmf_data_01r_pod10plus), "\n")
cat("0R/1R samples:", sum(mmf_data_01r_pod10plus$ACR_Group == "0R/1R"), "\n")
cat("2R+ samples:", sum(mmf_data_01r_pod10plus$ACR_Group == "2R+"), "\n\n")

# Add POD stratum
mmf_data_01r_pod10plus <- mmf_data_01r_pod10plus %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (>30 days)"))

# Run for each stratum
early_results_01r_pod10 <- analyze_stratum_01r(mmf_data_01r_pod10plus, "Early (10-30 days)")
print(early_results_01r_pod10$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Early_0R1R_vs_2R_POD10plus.png", 
       early_results_01r_pod10$plot, width = 10, height = 6, dpi = 300)

late_results_01r_pod10 <- analyze_stratum_01r(mmf_data_01r_pod10plus, "Late (>30 days)")
print(late_results_01r_pod10$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Late_0R1R_vs_2R_POD10plus.png", 
       late_results_01r_pod10$plot, width = 10, height = 6, dpi = 300)

# ============================================================================
# PART 4: POD-Stratified Analysis POD 10-45 ONLY
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA ANALYSIS (0R/1R vs 2R+, POD 10-45 ONLY) ===\n")
cat("NOTE: Removing very early (POD<10) and late (POD>45) samples\n")
cat("      to focus on the acute post-transplant window\n\n")

# Filter to POD 10-45 only
mmf_data_01r_pod10_45 <- patients_with_2R %>%
    filter(ACR %in% c("0R", "1R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    filter(POD >= 10 & POD <= 45) %>%
    select(H, POD, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R/1R"))

cat("Total samples with POD 10-45:", nrow(mmf_data_01r_pod10_45), "\n")
cat("0R/1R samples:", sum(mmf_data_01r_pod10_45$ACR_Group == "0R/1R"), "\n")
cat("2R+ samples:", sum(mmf_data_01r_pod10_45$ACR_Group == "2R+"), "\n\n")

# Add POD stratum
mmf_data_01r_pod10_45 <- mmf_data_01r_pod10_45 %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (31-45 days)"))

# Run for each stratum
early_results_01r_pod10_45 <- analyze_stratum_01r(mmf_data_01r_pod10_45, "Early (10-30 days)")
print(early_results_01r_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Early_0R1R_vs_2R_POD10_45.png", 
       early_results_01r_pod10_45$plot, width = 10, height = 6, dpi = 300)

late_results_01r_pod10_45 <- analyze_stratum_01r(mmf_data_01r_pod10_45, "Late (31-45 days)")
print(late_results_01r_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPA_Late_0R1R_vs_2R_POD10_45.png", 
       late_results_01r_pod10_45$plot, width = 10, height = 6, dpi = 300)

# ============================================================================
# SUMMARY: 0R/1R vs 2R+ Analysis
# ============================================================================

cat("\n=== SUMMARY: 0R/1R vs 2R+ ANALYSIS ===\n\n")

cat("1. All samples (including POD<10):\n")
cat("   Overall:     C18 p=", round(wilcox_c18_01r$p.value, 3), 
    ", HILIC p=", round(wilcox_hilic_01r$p.value, 3), "\n")
cat("   Early ≤30d:  C18 p=", round(early_results_01r$p_c18, 3), 
    ", HILIC p=", round(early_results_01r$p_hilic, 3), "\n")
cat("   Late >30d:   C18 p=", round(late_results_01r$p_c18, 3), 
    ", HILIC p=", round(late_results_01r$p_hilic, 3), "\n\n")

cat("2. Excluding POD<10:\n")
cat("   Early 10-30d: C18 p=", round(early_results_01r_pod10$p_c18, 3), 
    ", HILIC p=", round(early_results_01r_pod10$p_hilic, 3), "\n")
cat("   Late >30d:    C18 p=", round(late_results_01r_pod10$p_c18, 3), 
    ", HILIC p=", round(late_results_01r_pod10$p_hilic, 3), "\n\n")

cat("3. POD 10-45 only:\n")
cat("   Early 10-30d:  C18 p=", round(early_results_01r_pod10_45$p_c18, 3), 
    ", HILIC p=", round(early_results_01r_pod10_45$p_hilic, 3), "\n")
cat("   Late 31-45d:   C18 p=", round(late_results_01r_pod10_45$p_c18, 3), 
    ", HILIC p=", round(late_results_01r_pod10_45$p_hilic, 3), "\n\n")

cat("Plots saved to: Results2/Feedback_Analysis/POD_Stratified_new/\n")
cat("Overall comparison:\n")
cat("  - MPA_Overall_0R1R_vs_2R.png\n\n")
cat("Original stratification (all samples):\n")
cat("  - MPA_Early_0R1R_vs_2R.png\n")
cat("  - MPA_Late_0R1R_vs_2R.png\n\n")
cat("Excluding POD<10:\n")
cat("  - MPA_Early_0R1R_vs_2R_POD10plus.png\n")
cat("  - MPA_Late_0R1R_vs_2R_POD10plus.png\n\n")
cat("POD 10-45 only:\n")
cat("  - MPA_Early_0R1R_vs_2R_POD10_45.png\n")
cat("  - MPA_Late_0R1R_vs_2R_POD10_45.png\n\n")

# ============================================================================
# MYCOPHENOLIC ACID O-ACYL GLUCURONIDE ANALYSIS
# ============================================================================

cat("\n\n")
cat("================================================================================\n")
cat("MYCOPHENOLIC ACID O-ACYL GLUCURONIDE (MPA METABOLITE) ANALYSIS\n")
cat("================================================================================\n\n")

# ============================================================================
# PART 1: Overall MPAG comparison (0R vs 2R+)
# ============================================================================

cat("\n=== OVERALL MPA-GLUCURONIDE ANALYSIS (0R vs 2R+) ===\n")
cat("NOTE: Excluding 1R samples for clearer comparison\n\n")

# Get MPAG data - filter to keep only 0R and 2R+
mpa_data <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    select(H, POD, ACR, `Mycophenolic.acid.O.acyl.glucuronide..C18.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R"))

# Run Wilcoxon test
cat("Mycophenolic acid O-acyl glucuronide C18:\n")
wilcox_mpa <- wilcox.test(
    `Mycophenolic.acid.O.acyl.glucuronide..C18.` ~ ACR_Group, 
    data = mpa_data, 
    exact = FALSE
)
print(wilcox_mpa)

# Plot overall comparison
n_0r_mpa <- sum(mpa_data$ACR_Group == "0R")
n_2r_mpa <- sum(mpa_data$ACR_Group == "2R+")

p_mpa_overall <- ggplot(mpa_data, aes(x = ACR_Group, y = `Mycophenolic.acid.O.acyl.glucuronide..C18.`, fill = ACR_Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
    scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
    scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r_mpa, ")"),
                                "2R+" = paste0("2R+\n(n=", n_2r_mpa, ")"))) +
    labs(title = "MPAG Levels: 0R vs 2R+",
         subtitle = paste0("p=", round(wilcox_mpa$p.value, 3)),
         x = "ACR Group", y = "Level") +
    theme_minimal() +
    theme(legend.position = "none")

print(p_mpa_overall)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPAG_Overall.png", p_mpa_overall, 
       width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 2: POD-Stratified Analysis (Early ≤30 days vs Late >30 days)
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA-GLUCURONIDE ANALYSIS ===\n\n")

# Add POD stratum
mpa_data <- mpa_data %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (≤30 days)", "Late (>30 days)"))

# Function to run tests for one stratum (MPA version)
analyze_stratum_mpa <- function(data, stratum) {
    cat("--- ", stratum, " ---\n")
    
    stratum_data <- data %>% filter(POD_Stratum == stratum)
    n_0r <- sum(stratum_data$ACR_Group == "0R")
    n_2r <- sum(stratum_data$ACR_Group == "2R+")
    
    cat("Sample size:", nrow(stratum_data), 
        "(0R:", n_0r, ", 2R+:", n_2r, ")\n\n")
    
    # Wilcoxon test
    w_mpa <- wilcox.test(`Mycophenolic.acid.O.acyl.glucuronide..C18.` ~ ACR_Group, 
                         data = stratum_data, exact = FALSE)
    
    cat("MPAG: p =", round(w_mpa$p.value, 3), "\n\n")
    
    # Create plot
    p <- ggplot(stratum_data, aes(x = ACR_Group, y = `Mycophenolic.acid.O.acyl.glucuronide..C18.`, fill = ACR_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", color = "red", size = 3, shape = 18) +
        scale_fill_manual(values = c("0R" = "lightblue", "2R+" = "lightcoral")) +
        scale_x_discrete(labels = c("0R" = paste0("0R\n(n=", n_0r, ")"),
                                    "2R+" = paste0("2R+\n(n=", n_2r, ")"))) +
        labs(title = paste("MPAG Levels:", stratum),
             subtitle = paste0("p=", round(w_mpa$p.value, 3)),
             x = "ACR Group", y = "Level") +
        theme_minimal() +
        theme(legend.position = "none")
    
    return(list(plot = p, p_value = w_mpa$p.value))
}

# Run for each stratum
early_mpa <- analyze_stratum_mpa(mpa_data, "Early (≤30 days)")
print(early_mpa$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPAG_Early.png", early_mpa$plot, 
       width = 8, height = 6, dpi = 300)

late_mpa <- analyze_stratum_mpa(mpa_data, "Late (>30 days)")
print(late_mpa$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPAG_Late.png", late_mpa$plot, 
       width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 3: POD-Stratified Analysis EXCLUDING POD<10 (Early 10-30 vs Late >30)
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA-GLUCURONIDE ANALYSIS (EXCLUDING POD<10) ===\n")
cat("NOTE: Removing very early post-transplant samples (POD<10)\n\n")

# Filter out POD<10 samples
mpa_data_pod10plus <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    filter(POD >= 10) %>%
    select(H, POD, ACR, `Mycophenolic.acid.O.acyl.glucuronide..C18.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R"))

cat("Total samples with POD >= 10:", nrow(mpa_data_pod10plus), "\n")
cat("0R samples:", sum(mpa_data_pod10plus$ACR_Group == "0R"), "\n")
cat("2R+ samples:", sum(mpa_data_pod10plus$ACR_Group == "2R+"), "\n\n")

# Add POD stratum
mpa_data_pod10plus <- mpa_data_pod10plus %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (>30 days)"))

# Run for each stratum
early_mpa_pod10 <- analyze_stratum_mpa(mpa_data_pod10plus, "Early (10-30 days)")
print(early_mpa_pod10$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPAG_Early_POD10plus.png", early_mpa_pod10$plot, 
       width = 8, height = 6, dpi = 300)

late_mpa_pod10 <- analyze_stratum_mpa(mpa_data_pod10plus, "Late (>30 days)")
print(late_mpa_pod10$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPAG_Late_POD10plus.png", late_mpa_pod10$plot, 
       width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 4: POD-Stratified Analysis POD 10-45 (Early 10-30 vs Late 31-45)
# ============================================================================

cat("\n\n=== POD-STRATIFIED MPA-GLUCURONIDE ANALYSIS (POD 10-45 ONLY) ===\n")
cat("NOTE: Removing very early (POD<10) and late (POD>45) samples\n")
cat("      to focus on the acute post-transplant window\n\n")

# Filter to POD 10-45 only
mpa_data_pod10_45 <- patients_with_2R %>%
    filter(ACR %in% c("0R") | grepl("^2R", ACR, ignore.case = TRUE)) %>%
    filter(POD >= 10 & POD <= 45) %>%
    select(H, POD, ACR, `Mycophenolic.acid.O.acyl.glucuronide..C18.`) %>%
    mutate(ACR_Group = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R"))

cat("Total samples with POD 10-45:", nrow(mpa_data_pod10_45), "\n")
cat("0R samples:", sum(mpa_data_pod10_45$ACR_Group == "0R"), "\n")
cat("2R+ samples:", sum(mpa_data_pod10_45$ACR_Group == "2R+"), "\n\n")

# Add POD stratum
mpa_data_pod10_45 <- mpa_data_pod10_45 %>%
    mutate(POD_Stratum = ifelse(POD <= 30, "Early (10-30 days)", "Late (31-45 days)"))

# Run for each stratum
early_mpa_pod10_45 <- analyze_stratum_mpa(mpa_data_pod10_45, "Early (10-30 days)")
print(early_mpa_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPAG_Early_POD10_45.png", early_mpa_pod10_45$plot, 
       width = 8, height = 6, dpi = 300)

late_mpa_pod10_45 <- analyze_stratum_mpa(mpa_data_pod10_45, "Late (31-45 days)")
print(late_mpa_pod10_45$plot)
ggsave("Results2/Feedback_Analysis/POD_Stratified_new/MPAG_Late_POD10_45.png", late_mpa_pod10_45$plot, 
       width = 8, height = 6, dpi = 300)

# ============================================================================
# MPA-GLUCURONIDE SUMMARY
# ============================================================================

cat("\n=== MPA-GLUCURONIDE COMPREHENSIVE SUMMARY ===\n\n")

cat("1. All samples (including POD<10):\n")
cat("   Overall:     p=", round(wilcox_mpa$p.value, 3), "\n")
cat("   Early ≤30d:  p=", round(early_mpa$p_value, 3), "\n")
cat("   Late >30d:   p=", round(late_mpa$p_value, 3), "\n\n")

cat("2. Excluding POD<10:\n")
cat("   Early 10-30d: p=", round(early_mpa_pod10$p_value, 3), "\n")
cat("   Late >30d:    p=", round(late_mpa_pod10$p_value, 3), "\n\n")

cat("3. POD 10-45 only (excluding POD<10 and POD>45):\n")
cat("   Early 10-30d:  p=", round(early_mpa_pod10_45$p_value, 3), "\n")
cat("   Late 31-45d:   p=", round(late_mpa_pod10_45$p_value, 3), "\n\n")

cat("MPAG plots saved to: Results2/Feedback_Analysis/POD_Stratified_new/\n")
cat("Original stratification:\n")
cat("  - MPAG_Overall.png\n")
cat("  - MPAG_Early.png\n")
cat("  - MPAG_Late.png\n\n")
cat("Excluding POD<10:\n")
cat("  - MPAG_Early_POD10plus.png\n")
cat("  - MPAG_Late_POD10plus.png\n\n")
cat("POD 10-45 only:\n")
cat("  - MPAG_Early_POD10_45.png\n")
cat("  - MPAG_Late_POD10_45.png\n\n")



