# Focusing on MMF and Prograf metabolites and relevant metabolite pathways 

# run 00_source first
# source("Immunosuppressive metabolites feature table/Scripts/00_source")

# Ensure output is visible
options(width = 120)  # Wider console output

# Q1: in patients that develop 2R or greater rejection, are the MMF/prograf levels and prograf levels different compared to those with 0R?

## 2R+ vs 0R
#! need to ask about Prograf data as we don't have this metabolite in the dataset

#+ Create MMF and Prograf tables for 0R and 2R+ groups
acr_0r_log2_MMF_C18 <- acr_0r_log2 %>%
    select(`Mycophenolate..C18.`)

acr_0r_log2_MMF_HILIC <- acr_0r_log2 %>%
    select(`Mycophenolate..HILIC.`)

acr_2r_log2_MMF_C18 <- acr_2r_log2 %>%
    select(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`)

acr_2r_log2_MMF_HILIC <- acr_2r_log2 %>%
    select(`Mycophenolate..HILIC.`)

# Perform Wilcoxon tests for MMF_C18 and MMF_HILIC levels between 0R and 2R+ groups
wilcox_MMF_C18 <- wilcox.test(acr_0r_log2_MMF_C18$`Mycophenolate..C18.`, 
                              acr_2r_log2_MMF_C18$`Mycophenolate..C18.`)

wilcox_MMF_HILIC <- wilcox.test(acr_0r_log2_MMF_HILIC$`Mycophenolate..HILIC.`, 
                                 acr_2r_log2_MMF_HILIC$`Mycophenolate..HILIC.`)
# Output results
cat("Wilcoxon test for MMF C18 levels between 0R and 2R+ groups:\n")
cat("W =", wilcox_MMF_C18$statistic, ", p-value =", wilcox_MMF_C18$p.value, "\n\n")

cat("Wilcoxon test for MMF HILIC levels between 0R and 2R+ groups:\n")
cat("W =", wilcox_MMF_HILIC$statistic, ", p-value =", wilcox_MMF_HILIC$p.value, "\n\n")

# Visualize MMF levels in 0R vs 2R+ groups
mmf_data <- nonlog2_data %>%
    filter(H %in% patients_with_2R$H) %>%
    select(H, ACR, `Mycophenolate..C18.`, `Mycophenolate..HILIC.`) %>%
    pivot_longer(cols = c(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`), 
                 names_to = "Metabolite", 
                 values_to = "Level")
p_mmf <- mmf_data %>%
    ggplot(aes(x = ACR, y = Level, fill = ACR)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    facet_wrap(~ Metabolite, scales = "free_y") +
    labs(title = "MMF Levels by ACR Status",
         y = "MMF Level (Non-Log2)",
         x = "ACR Status") +
    theme_minimal() +
    theme(legend.position = "none")

ggsave("MMF_Levels_by_ACR.png", p_mmf, width = 10, height = 6, dpi = 300)


#=== POD-STRATIFIED MMF ANALYSIS ===

cat("\n=== POD-STRATIFIED MMF ANALYSIS ===\n")

# Define POD strata: POD <= 30 and POD > 30
pod_strata <- list(
    "Early (≤30 days)" = c(0, 30),
    "Late (>30 days)" = c(31, Inf)
)

# Function to perform POD-stratified analysis and create plots
perform_pod_mmf_analysis <- function(data, strata_list) {
    
    all_results <- list()
    all_plots <- list()
    
    for(stratum_name in names(strata_list)) {
        range_vals <- strata_list[[stratum_name]]
        
        cat("\n--- Analyzing", stratum_name, "(POD", range_vals[1], "to", 
            ifelse(is.infinite(range_vals[2]), "max", range_vals[2]), ") ---\n")
        
        # Filter data for this POD stratum
        stratum_data <- data %>%
            filter(POD >= range_vals[1] & POD <= range_vals[2])
        
        cat("Total samples in stratum:", nrow(stratum_data), "\n")
        
        if(nrow(stratum_data) < 10) {
            cat("Skipping stratum - insufficient samples\n")
            next
        }
        
        # Create ACR groups for MMF metabolites only
        acr_0r_mmf <- stratum_data %>%
            filter(ACR %in% c("0R", "1R")) %>%
            select(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`)
        
        acr_2r_mmf <- stratum_data %>%
            filter(grepl("^2R", ACR, ignore.case = TRUE)) %>%
            select(`Mycophenolate..C18.`, `Mycophenolate..HILIC.`)
        
        cat("Sample sizes: 0R/1R =", nrow(acr_0r_mmf), ", 2R+ =", nrow(acr_2r_mmf), "\n")
        
        if(nrow(acr_0r_mmf) < 3 | nrow(acr_2r_mmf) < 3) {
            cat("Skipping stratum - insufficient samples per group\n")
            next
        }
        
        # Perform Wilcoxon tests for MMF metabolites
        mmf_results <- lapply(names(acr_0r_mmf), function(metabolite) {
            wilcox_test <- wilcox.test(acr_0r_mmf[[metabolite]], 
                                       acr_2r_mmf[[metabolite]], 
                                       exact = FALSE)
            tibble(
                Stratum = stratum_name,
                Metabolite = metabolite,
                Median_0_1R = median(acr_0r_mmf[[metabolite]], na.rm = TRUE),
                Median_2R_plus = median(acr_2r_mmf[[metabolite]], na.rm = TRUE),
                median_diff = median(acr_2r_mmf[[metabolite]], na.rm = TRUE) - 
                             median(acr_0r_mmf[[metabolite]], na.rm = TRUE),
                p_value = wilcox_test$p.value
            )
        })
        mmf_results <- bind_rows(mmf_results)
        
        # Print results
        cat("\n=== MMF Wilcoxon test results ===\n")
        print(mmf_results)
        cat("\n")
        cat("Column names in mmf_results:", paste(mmf_results$Metabolite, collapse=", "), "\n")
        cat("Column names in acr_0r_mmf:", paste(names(acr_0r_mmf), collapse=", "), "\n")
        cat("\n")
        
        # Create plots using plot_ttest_bars function from 03a_plotting_functions
        if(exists("plot_ttest_bars")) {
            cat("Creating plots for", nrow(mmf_results), "metabolites...\n")
            cat("Metabolites to plot:", paste(mmf_results$Metabolite, collapse = ", "), "\n")
            cat("Data columns available:", paste(names(acr_0r_mmf), collapse = ", "), "\n")
            
            stratum_plots <- plot_ttest_bars(
                sig_results = mmf_results,
                group1_data = acr_0r_mmf,
                group2_data = acr_2r_mmf,
                use_median = TRUE  # Use median for non-parametric approach
            )
            
            cat("Generated", length(stratum_plots), "plots\n")
            cat("Plot names:", paste(names(stratum_plots), collapse = ", "), "\n")
            
            # Save plots for this stratum
            dir.create("Results2/Feedback_Analysis/POD_Stratified", 
                      showWarnings = FALSE, recursive = TRUE)
            
            if(length(stratum_plots) > 0) {
                for(metabolite in names(stratum_plots)) {
                    filename <- paste0("Results2/Feedback_Analysis/POD_Stratified/",
                                     gsub("[^A-Za-z0-9]", "_", stratum_name), "_",
                                     gsub("[^A-Za-z0-9]", "_", metabolite), "_plot.png")
                    ggsave(filename, stratum_plots[[metabolite]], 
                          width = 8, height = 6, dpi = 300)
                    cat("Saved plot:", filename, "\n")
                }
                
                all_plots[[stratum_name]] <- stratum_plots
            } else {
                cat("WARNING: No plots were generated. Check metabolite names.\n")
            }
        } else {
            cat("WARNING: plot_ttest_bars() function not found. Run 00_source first.\n")
        }
        
        all_results[[stratum_name]] <- mmf_results
    }
    
    return(list(results = all_results, plots = all_plots))
}

# Run POD-stratified analysis
cat("\nRunning POD-stratified MMF analysis...\n")
pod_mmf_analysis <- perform_pod_mmf_analysis(
    data = patients_with_2R,  # Use the subset with 2R+ patients
    strata_list = pod_strata
)



# Export combined results
if(length(pod_mmf_analysis$results) > 0) {
    combined_pod_results <- bind_rows(pod_mmf_analysis$results)
    write_csv(combined_pod_results, 
             "Results2/Feedback_Analysis/MMF_POD_Stratified_Results.csv")
    cat("\nPOD-stratified results saved to: Results2/Feedback_Analysis/MMF_POD_Stratified_Results.csv\n")
} else {
    cat("\nWARNING: No POD-stratified results to export.\n")
}

cat("\n=== POD-STRATIFIED MMF ANALYSIS COMPLETE ===\n")

# Flush output to console
flush.console()

# Print summary
if(exists("pod_mmf_analysis") && length(pod_mmf_analysis$results) > 0) {
    cat("\n=== SUMMARY ===\n")
    cat("Number of POD strata analyzed:", length(pod_mmf_analysis$results), "\n")
    cat("Total plots generated:", sum(sapply(pod_mmf_analysis$plots, length)), "\n")
    cat("\nResults files saved to: Results2/Feedback_Analysis/\n")
} else {
    cat("\nWARNING: No results were generated. Check your data and POD ranges.\n")
}

# Display all plots in viewer
if(!is.null(pod_mmf_analysis$plots) && length(pod_mmf_analysis$plots) > 0) {
    cat("\n=== DISPLAYING ALL POD-STRATIFIED PLOTS ===\n")
    cat("Total plots generated:", sum(sapply(pod_mmf_analysis$plots, length)), "\n")
    cat("All plots saved to: Results2/Feedback_Analysis/POD_Stratified/\n\n")
    
    # Display all plots (they go into plot history)
    for(stratum_name in names(pod_mmf_analysis$plots)) {
        for(metabolite in names(pod_mmf_analysis$plots[[stratum_name]])) {
            cat("Displaying:", stratum_name, "-", metabolite, "\n")
            print(pod_mmf_analysis$plots[[stratum_name]][[metabolite]])
        }
    }
    
    cat("\n✓ All", sum(sapply(pod_mmf_analysis$plots, length)), "plots displayed.\n")
    cat("✓ Use the ← → arrow buttons in the Plots pane to navigate.\n")
    cat("✓ Files saved to: Results2/Feedback_Analysis/POD_Stratified/\n\n")
}



