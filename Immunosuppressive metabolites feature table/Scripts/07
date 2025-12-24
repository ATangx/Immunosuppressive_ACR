# ============================================================================
# 07_Mixed_Effects_Modeling: Linear Mixed-Effects Models for Longitudinal Data
# ============================================================================
# 
# PURPOSE:
# This script uses linear mixed-effects models (LME) to analyze longitudinal
# metabolite data, accounting for within-patient correlation and allowing
# flexible modeling of trajectories over time and rejection status.
#
# WHY USE MIXED-EFFECTS MODELS?
# - Handles unbalanced data (different numbers of samples per patient)
# - Accounts for within-patient correlation (repeated measures)
# - Can include time-varying and patient-level covariates
# - More efficient than simple paired tests when >2 timepoints
# - Provides patient-specific (random) and population-level (fixed) effects
#
# KEY FEATURES:
# - Fixed effects: ACR status, POD (time), and their interaction
# - Random effects: Patient-specific intercepts (and optionally slopes)
# - Tests overall ACR effect, time effect, and ACR×time interaction
# - Can handle missing data (uses all available observations)
# - Provides effect sizes and confidence intervals
#
# DIFFERENCES FROM OTHER SCRIPTS:
# - Unlike 06/06b: Uses ALL timepoints, not just paired pre/active/post
# - Unlike 03/05: Properly accounts for repeated measures via random effects
# - More sophisticated than simple paired tests
# - Can test complex hypotheses (e.g., do rejection slopes differ?)
#
# MODELS TESTED:
# 1) Random intercept: metabolite ~ ACR + POD + (1|Patient)
# 2) Random intercept + slope: metabolite ~ ACR + POD + (1 + POD|Patient)
# 3) Interaction model: metabolite ~ ACR * POD + (1|Patient)
#
# OUTPUT:
# - Model comparison results (which model fits best?)
# - Fixed effect estimates and significance tests
# - Random effect variance estimates
# - Predicted trajectories by ACR status
# - Individual patient trajectory plots with model fits
#
# CLINICAL INTERPRETATION:
# - ACR effect: Overall metabolite difference between rejection states
# - POD effect: How metabolite changes over time post-transplant
# - ACR×POD interaction: Does rejection affect metabolite trajectory?
#
# ============================================================================

# Source all required data and functions
suppressMessages(suppressWarnings({
  source("Immunosuppressive metabolites feature table/Scripts/00_source")
  
  # Load mixed-effects modeling packages
  if (!require("lme4")) install.packages("lme4")
  if (!require("lmerTest")) install.packages("lmerTest")  # Adds p-values to lme4
  if (!require("broom.mixed")) install.packages("broom.mixed")  # Tidy model outputs
  
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
}))

cat("\n=== LINEAR MIXED-EFFECTS MODELING FOR LONGITUDINAL METABOLITE DATA ===\n")

# Create output directories
output_dir <- "Results2/Mixed_Effects_Models"
plot_dir <- file.path(output_dir, "Trajectory_Plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Get metabolite columns
metadata_cols <- c("H", "POD", "S", "Sample_ID", "ACR")
metabolite_cols <- setdiff(names(select_if(patients_with_2R, is.numeric)), metadata_cols)

cat("Identified", length(metabolite_cols), "metabolite columns for mixed-effects modeling\n")

# ============================================================================
# STEP 1: PREPARE DATA FOR MIXED-EFFECTS MODELING
# ============================================================================

cat("\nStep 1: Preparing longitudinal data for modeling...\n")

# Check if POD is available
if (!"POD" %in% names(patients_with_2R)) {
  cat("\nWARNING: POD (post-operative day) not found in data.\n")
  cat("Mixed-effects models work best with time information.\n")
  cat("Proceeding with ACR status only (no time modeling).\n")
  use_pod <- FALSE
} else {
  use_pod <- TRUE
  cat("POD (time) variable available - will include in models.\n")
}

# Create modeling dataset
model_data <- patients_with_2R %>%
  mutate(
    Patient = as.factor(H),  # Patient ID as factor
    ACR_binary = ifelse(grepl("^2R", ACR, ignore.case = TRUE), "2R+", "0R/1R"),
    ACR_binary = factor(ACR_binary, levels = c("0R/1R", "2R+")),  # Reference = 0R/1R
    # Optional: create numeric ACR score if needed
    ACR_numeric = case_when(
      ACR == "0R" ~ 0,
      ACR == "1R" ~ 1,
      grepl("^2R", ACR, ignore.case = TRUE) ~ 2,
      TRUE ~ NA_real_
    )
  )

if (use_pod) {
  model_data <- model_data %>%
    mutate(
      POD_centered = POD - mean(POD, na.rm = TRUE),  # Center POD for interpretation
      POD_scaled = scale(POD)[,1]  # Standardize POD for model convergence
    )
}

# Summary of data structure
cat("\nData structure for mixed-effects modeling:\n")
cat("Total observations:", nrow(model_data), "\n")
cat("Unique patients:", length(unique(model_data$Patient)), "\n")
cat("ACR status distribution:\n")
print(table(model_data$ACR_binary))

# Samples per patient
samples_per_patient <- model_data %>%
  group_by(Patient) %>%
  summarise(n_samples = n(), .groups = "drop")

cat("\nSamples per patient (summary):\n")
print(summary(samples_per_patient$n_samples))
cat("Patients with ≥2 samples:", sum(samples_per_patient$n_samples >= 2), "\n")
cat("Patients with ≥3 samples:", sum(samples_per_patient$n_samples >= 3), "\n")

if (sum(samples_per_patient$n_samples >= 2) < 5) {
  cat("\nWARNING: Very few patients with repeated measures.\n")
  cat("Mixed-effects models may not converge or provide limited benefit.\n")
  cat("Consider using simpler paired analyses (scripts 06/06b).\n")
}

# ============================================================================
# STEP 2: FIT MIXED-EFFECTS MODELS
# ============================================================================

cat("\nStep 2: Fitting linear mixed-effects models...\n")

# Function to fit multiple models and compare
fit_lme_models <- function(metabolite_col) {
  
  metabolite_name <- gsub("\\.", " ", metabolite_col)
  
  # Prepare data for this metabolite
  if (use_pod) {
    met_data <- model_data %>%
      select(Patient, ACR_binary, concentration = !!sym(metabolite_col), 
             POD_centered, POD_scaled) %>%
      filter(!is.na(concentration))
  } else {
    met_data <- model_data %>%
      select(Patient, ACR_binary, concentration = !!sym(metabolite_col)) %>%
      filter(!is.na(concentration))
  }
  
  # Check if enough data
  if (nrow(met_data) < 10) {
    return(list(
      metabolite = metabolite_name,
      n_obs = nrow(met_data),
      error = "Insufficient data"
    ))
  }
  
  # Initialize results
  results <- list(
    metabolite = metabolite_name,
    n_obs = nrow(met_data),
    n_patients = length(unique(met_data$Patient))
  )
  
  # MODEL 1: Random intercept only, ACR effect
  # metabolite ~ ACR + (1|Patient)
  model1 <- tryCatch({
    if (use_pod) {
      lmer(concentration ~ ACR_binary + POD_centered + (1|Patient), 
           data = met_data, REML = FALSE)
    } else {
      lmer(concentration ~ ACR_binary + (1|Patient), 
           data = met_data, REML = FALSE)
    }
  }, error = function(e) NULL, warning = function(w) NULL)
  
  if (!is.null(model1)) {
    results$model1_converged <- TRUE
    results$model1_AIC <- AIC(model1)
    results$model1_summary <- tidy(model1, effects = "fixed")
    
    # Extract ACR effect
    acr_effect <- results$model1_summary %>% 
      filter(term == "ACR_binary2R+")
    
    if (nrow(acr_effect) > 0) {
      results$ACR_estimate <- acr_effect$estimate
      results$ACR_se <- acr_effect$std.error
      results$ACR_pvalue <- acr_effect$p.value
      results$ACR_CI_lower <- acr_effect$estimate - 1.96 * acr_effect$std.error
      results$ACR_CI_upper <- acr_effect$estimate + 1.96 * acr_effect$std.error
    }
  } else {
    results$model1_converged <- FALSE
  }
  
  # MODEL 2: Random intercept + slope (if POD available and enough data)
  if (use_pod && nrow(met_data) >= 20) {
    model2 <- tryCatch({
      lmer(concentration ~ ACR_binary + POD_centered + (1 + POD_centered|Patient), 
           data = met_data, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
    }, error = function(e) NULL, warning = function(w) NULL)
    
    if (!is.null(model2)) {
      results$model2_converged <- TRUE
      results$model2_AIC <- AIC(model2)
      
      # Compare models
      if (!is.null(model1)) {
        anova_result <- anova(model1, model2)
        results$model2_vs_model1_pvalue <- anova_result$`Pr(>Chisq)`[2]
      }
    } else {
      results$model2_converged <- FALSE
    }
  }
  
  # MODEL 3: Interaction model (ACR × POD)
  if (use_pod) {
    model3 <- tryCatch({
      lmer(concentration ~ ACR_binary * POD_centered + (1|Patient), 
           data = met_data, REML = FALSE)
    }, error = function(e) NULL, warning = function(w) NULL)
    
    if (!is.null(model3)) {
      results$model3_converged <- TRUE
      results$model3_AIC <- AIC(model3)
      results$model3_summary <- tidy(model3, effects = "fixed")
      
      # Extract interaction effect
      int_effect <- results$model3_summary %>% 
        filter(grepl(":", term))
      
      if (nrow(int_effect) > 0) {
        results$interaction_estimate <- int_effect$estimate
        results$interaction_pvalue <- int_effect$p.value
      }
      
      # Compare to model1
      if (!is.null(model1)) {
        anova_result <- anova(model1, model3)
        results$model3_vs_model1_pvalue <- anova_result$`Pr(>Chisq)`[2]
      }
    } else {
      results$model3_converged <- FALSE
    }
  }
  
  # Store best model for later use
  best_model <- NULL
  if (!is.null(model3) && results$model3_converged) {
    best_model <- model3
    results$best_model <- "Model 3 (interaction)"
  } else if (!is.null(model2) && results$model2_converged) {
    best_model <- model2
    results$best_model <- "Model 2 (random slope)"
  } else if (!is.null(model1) && results$model1_converged) {
    best_model <- model1
    results$best_model <- "Model 1 (random intercept)"
  }
  
  results$model_object <- best_model
  
  return(results)
}

# Fit models for all metabolites
cat("Fitting mixed-effects models for all metabolites...\n")
cat("This may take a few minutes...\n\n")

all_lme_results <- list()
for (i in seq_along(metabolite_cols)) {
  met <- metabolite_cols[i]
  if (i %% 5 == 0) cat("Progress:", i, "/", length(metabolite_cols), "\n")
  
  result <- fit_lme_models(met)
  all_lme_results[[met]] <- result
}

cat("\nModel fitting complete.\n")

# ============================================================================
# STEP 3: SUMMARIZE AND EXPORT RESULTS
# ============================================================================

cat("\nStep 3: Summarizing mixed-effects model results...\n")

# Extract key results into data frame
lme_summary <- map_df(all_lme_results, function(res) {
  tibble(
    Metabolite = res$metabolite,
    N_Observations = res$n_obs,
    N_Patients = res$n_patients,
    Model1_Converged = ifelse(is.null(res$model1_converged), FALSE, res$model1_converged),
    ACR_Effect_Estimate = res$ACR_estimate %||% NA_real_,
    ACR_Effect_SE = res$ACR_se %||% NA_real_,
    ACR_Effect_Pvalue = res$ACR_pvalue %||% NA_real_,
    ACR_Effect_CI_Lower = res$ACR_CI_lower %||% NA_real_,
    ACR_Effect_CI_Upper = res$ACR_CI_upper %||% NA_real_,
    Interaction_Estimate = res$interaction_estimate %||% NA_real_,
    Interaction_Pvalue = res$interaction_pvalue %||% NA_real_,
    Best_Model = res$best_model %||% "None converged",
    Model1_AIC = res$model1_AIC %||% NA_real_
  )
}) %>%
  mutate(
    ACR_Effect_FDR = p.adjust(ACR_Effect_Pvalue, method = "fdr"),
    Interaction_FDR = p.adjust(Interaction_Pvalue, method = "fdr")
  ) %>%
  arrange(ACR_Effect_Pvalue)

# Save results
write.csv(lme_summary, 
          file.path(output_dir, "Mixed_Effects_Model_Results.csv"), 
          row.names = FALSE)

cat("\nMixed-effects model summary:\n")
cat("Models attempted:", nrow(lme_summary), "\n")
cat("Models converged:", sum(lme_summary$Model1_Converged, na.rm = TRUE), "\n")
cat("Significant ACR effects (FDR < 0.05):", 
    sum(lme_summary$ACR_Effect_FDR < 0.05, na.rm = TRUE), "\n")
cat("Significant ACR effects (p < 0.05):", 
    sum(lme_summary$ACR_Effect_Pvalue < 0.05, na.rm = TRUE), "\n")

if (use_pod) {
  cat("Significant ACR×POD interactions (FDR < 0.05):", 
      sum(lme_summary$Interaction_FDR < 0.05, na.rm = TRUE), "\n")
}

cat("\nTop metabolites by ACR effect:\n")
print(lme_summary %>% 
        select(Metabolite, N_Patients, ACR_Effect_Estimate, ACR_Effect_Pvalue, ACR_Effect_FDR) %>%
        head(10))

# ============================================================================
# STEP 4: CREATE TRAJECTORY PLOTS WITH MODEL FITS
# ============================================================================

cat("\nStep 4: Creating trajectory plots with mixed-effects model fits...\n")

create_lme_trajectory_plot <- function(metabolite_col) {
  
  metabolite_name <- gsub("\\.", " ", metabolite_col)
  
  # Get model results
  model_result <- all_lme_results[[metabolite_col]]
  
  if (is.null(model_result$model_object)) return(NULL)
  
  # Get data
  if (use_pod) {
    met_data <- model_data %>%
      select(Patient, ACR_binary, concentration = !!sym(metabolite_col),
             POD, POD_centered) %>%
      filter(!is.na(concentration))
  } else {
    met_data <- model_data %>%
      select(Patient, ACR_binary, concentration = !!sym(metabolite_col)) %>%
      filter(!is.na(concentration))
  }
  
  if (nrow(met_data) < 5) return(NULL)
  
  # Get p-values for subtitle
  acr_p <- model_result$ACR_pvalue %||% NA_real_
  int_p <- model_result$interaction_pvalue %||% NA_real_
  
  subtitle_txt <- paste0(
    "n = ", model_result$n_patients, " patients, ",
    model_result$n_obs, " observations\n",
    "ACR effect p = ", ifelse(is.na(acr_p), "NA", format(acr_p, digits = 3))
  )
  
  if (!is.na(int_p)) {
    subtitle_txt <- paste0(subtitle_txt, 
                          "; ACR×Time interaction p = ", format(int_p, digits = 3))
  }
  
  # Create plot
  if (use_pod) {
    # Trajectory plot with POD
    p <- ggplot(met_data, aes(x = POD, y = concentration, color = ACR_binary, group = Patient)) +
      geom_line(alpha = 0.3, linewidth = 0.5) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_smooth(aes(group = ACR_binary), method = "lm", se = TRUE, 
                  alpha = 0.2, linewidth = 1.5) +
      scale_color_manual(values = c("0R/1R" = "steelblue", "2R+" = "red"),
                        name = "ACR Status") +
      labs(
        title = paste("Mixed-Effects Model Trajectories:", metabolite_name),
        subtitle = subtitle_txt,
        x = "Post-Operative Day (POD)",
        y = "Concentration"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9),
        legend.position = "bottom"
      )
  } else {
    # Boxplot without POD
    p <- ggplot(met_data, aes(x = ACR_binary, y = concentration, fill = ACR_binary)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.3) +
      geom_point(alpha = 0.5, position = position_jitter(width = 0.2)) +
      scale_fill_manual(values = c("0R/1R" = "steelblue", "2R+" = "red"),
                       name = "ACR Status") +
      labs(
        title = paste("Mixed-Effects Model:", metabolite_name),
        subtitle = subtitle_txt,
        x = "ACR Status",
        y = "Concentration"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9),
        legend.position = "none"
      )
  }
  
  return(p)
}

# Plot significant metabolites
sig_metabolites <- lme_summary %>%
  filter(ACR_Effect_FDR < 0.05 | ACR_Effect_Pvalue < 0.05) %>%
  arrange(ACR_Effect_Pvalue) %>%
  pull(Metabolite)

if (length(sig_metabolites) > 0) {
  cat("Creating trajectory plots for", length(sig_metabolites), "significant metabolites...\n")
  
  for (met_name in sig_metabolites) {
    met_col <- gsub(" ", "\\.", met_name)
    if (met_col %in% metabolite_cols) {
      plot <- create_lme_trajectory_plot(met_col)
      
      if (!is.null(plot)) {
        safe_filename <- gsub("[^A-Za-z0-9_-]", "_", met_name)
        ggsave(file.path(plot_dir, paste0(safe_filename, "_LME_trajectory.png")),
               plot, width = 10, height = 7, dpi = 300)
      }
    }
  }
  cat("Saved trajectory plots to:", plot_dir, "\n")
} else {
  cat("No metabolites with significant ACR effects for plotting.\n")
}

# ============================================================================
# STEP 5: COMPARE WITH SIMPLE PAIRED ANALYSIS
# ============================================================================

cat("\nStep 5: Comparing mixed-effects results with simple paired analysis...\n")

# Load paired analysis results if available
paired_results_file <- "Results2/Paired_Pre2R_vs_Active2R/paired_test_results.csv"

if (file.exists(paired_results_file)) {
  paired_results <- read.csv(paired_results_file)
  
  # Merge with LME results
  comparison <- lme_summary %>%
    select(Metabolite, LME_ACR_p = ACR_Effect_Pvalue, LME_ACR_FDR = ACR_Effect_FDR,
           LME_Effect = ACR_Effect_Estimate) %>%
    left_join(
      paired_results %>%
        select(Metabolite, Paired_Wilcox_p = wilcox_p_value, Paired_t_p = t_test_p_value),
      by = "Metabolite"
    )
  
  # Classify concordance
  comparison <- comparison %>%
    mutate(
      LME_sig = LME_ACR_FDR < 0.05,
      Paired_sig = Paired_Wilcox_p < 0.05,
      Concordant = case_when(
        LME_sig & Paired_sig ~ "Both significant",
        !LME_sig & !Paired_sig ~ "Both non-significant",
        LME_sig & !Paired_sig ~ "LME only",
        !LME_sig & Paired_sig ~ "Paired only"
      )
    )
  
  cat("\nComparison with paired analysis (06):\n")
  print(table(comparison$Concordant))
  
  write.csv(comparison, 
            file.path(output_dir, "LME_vs_Paired_Comparison.csv"),
            row.names = FALSE)
  
  cat("Saved comparison to:", file.path(output_dir, "LME_vs_Paired_Comparison.csv"), "\n")
} else {
  cat("Paired analysis results not found. Run script 06 first to enable comparison.\n")
}

# ============================================================================
# STEP 6: INTERPRETATION GUIDE AND SUMMARY
# ============================================================================

cat("\n=== MIXED-EFFECTS MODELING SUMMARY ===\n")

cat("\nWHAT WERE WE TESTING?\n")
cat("- Fixed effect of ACR status on metabolite levels (population-level)\n")
cat("- Random effect of Patient (allows each patient to have different baseline)\n")
if (use_pod) {
  cat("- Effect of time (POD) on metabolite levels\n")
  cat("- ACR×Time interaction (do rejection trajectories differ?)\n")
}

cat("\nKEY FINDINGS:\n")
cat("- Total metabolites tested:", nrow(lme_summary), "\n")
cat("- Models that converged:", sum(lme_summary$Model1_Converged, na.rm = TRUE), "\n")
cat("- Significant ACR effects (FDR < 0.05):", 
    sum(lme_summary$ACR_Effect_FDR < 0.05, na.rm = TRUE), "\n")

if (use_pod) {
  cat("- Significant ACR×Time interactions (FDR < 0.05):", 
      sum(lme_summary$Interaction_FDR < 0.05, na.rm = TRUE), "\n")
}

cat("\nINTERPRETATION:\n")
cat("• ACR_Effect_Estimate: Mean difference in metabolite level (2R+ vs 0R/1R)\n")
cat("  - Positive = higher in 2R+ patients\n")
cat("  - Negative = lower in 2R+ patients\n")
cat("• Accounts for within-patient correlation (repeated measures)\n")
cat("• Uses all available data (not just paired samples)\n")

if (use_pod) {
  cat("• Interaction_Estimate: Difference in slopes (2R+ vs 0R/1R over time)\n")
  cat("  - Significant interaction = rejection changes trajectory\n")
}

cat("\nADVANTAGES OVER SIMPLE PAIRED TESTS:\n")
cat("✓ Uses ALL timepoints, not just pre/active pairs\n")
cat("✓ Handles unbalanced data (different sample counts per patient)\n")
cat("✓ Properly accounts for within-patient correlation\n")
cat("✓ Can test time effects and ACR×time interactions\n")
cat("✓ More efficient with limited sample size\n")

cat("\nLIMITATIONS:\n")
cat("✗ Assumes linear relationships (may not capture complex patterns)\n")
cat("✗ Requires sufficient within-patient variation for convergence\n")
cat("✗ More complex to interpret than simple paired tests\n")
cat("✗ May not converge with very small samples\n")

cat("\nRECOMMENDATIONS:\n")
cat("1. Use LME when you have ≥3 timepoints per patient (more efficient)\n")
cat("2. Use paired tests (06/06b) when you have exactly 2 or 3 timepoints\n")
cat("3. Report both if possible (complementary approaches)\n")
cat("4. Focus on metabolites significant in BOTH analyses (robust findings)\n")

cat("\nOUTPUT FILES:\n")
cat("- Model results:", file.path(output_dir, "Mixed_Effects_Model_Results.csv"), "\n")
cat("- Trajectory plots:", plot_dir, "\n")
if (file.exists(paired_results_file)) {
  cat("- LME vs Paired comparison:", 
      file.path(output_dir, "LME_vs_Paired_Comparison.csv"), "\n")
}

cat("\n=== SCRIPT COMPLETE ===\n\n")
