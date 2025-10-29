#* Import data preprocessed in SM4, manipulate
  #+ Targeted
    rename_vector <- read_csv("MSMICA_feature_key.csv") %>%
      select(short_name, Name) %>%
      distinct() %>% # just in case of duplicate entries
      deframe() # makes a named vector: names are File_Name, values are true_sample_ID
    steroids <- read_excel("/Users/jdp2019/Library/CloudStorage/OneDrive-Emory/Research/Manuscripts and Projects/Active Projects/Chan Lab/TPMO/Metabolomics - ACR_Clay/study_metadata.xlsx") %>%
      select(Steroids,original_sample_ID) %>%
      rename(Sample_ID = original_sample_ID) %>%
      rename(steroids = Steroids)
    HILIC_C18_TFT_steroids <- read_csv("full_TFT_with_meta.csv") %>%
      left_join(steroids,by = "Sample_ID") %>%
      mutate(ACR = if_else(is.na(ACR), "Preop", ACR)) %>%
      select(-c(sex,Sample_Type,ACR_relation,ANOVA_pre_active_post)) %>%
      select(Sample_ID,steroids, ACR, Sample_no, everything()) %>%
      rename_with(.fn = ~ rename_vector[.x], .cols = all_of(names(rename_vector))) %>%
      mutate(ACR = as.factor(ACR)) %>%
      mutate(ACR = recode(
        ACR,
        `Preop` = "Preop",
        `0R` = "0R",
        `1R` = "1R",
        `2R` = "2R+",
        `3R` = "2R+"
      )) %>%
      mutate(Patient_no = paste0("H",Patient_no,sep = ""))
  #+ Untargeted
    HILIC_C18_UFT_steroids <- read_csv("full_UFT_with_meta.csv") %>%
      left_join(steroids,by = "Sample_ID") %>%
      mutate(ACR = if_else(is.na(ACR), "Preop", ACR)) %>%
      select(-c(sex,Sample_Type,ACR_relation,ANOVA_pre_active_post)) %>%
      mutate(ACR = as.factor(ACR)) %>%
      select(Sample_ID,steroids, ACR, Sample_no, everything()) %>%
      mutate(ACR = recode(
        ACR,
        `Preop` = "Preop",
        `0R` = "0R",
        `1R` = "1R",
        `2R` = "2R+",
        `3R` = "2R+"
      ))  %>%
      mutate(Patient_no = paste0("H", Patient_no, sep = ""))
      # %>% {
      #   meta_cols <- c("Sample_ID", "steroids", "ACR", "Sample_no", "Patient_no", "paired_POD_EMB_day")
      #   metabolite_cols <- setdiff(names(.), meta_cols)

      #   good_metabolites <- metabolite_cols[!sapply(metabolite_cols, function(col) {
      #     vals <- .[[col]]
      #     min_val <- min(vals, na.rm = TRUE)
      #     mean(vals == min_val) > 0.25
      #   })]

      #   select(., all_of(meta_cols), all_of(good_metabolites))
      # }
#* 0R vs 2R analysis
  #+ PLS-DA
    #- PLS-DA Function
      make_plsda_plot <- function(X, Y, label, filename) {
        plsda_model <- plsda(X, Y, ncomp = 2)
        scores <- plsda_model$variates$X
        explained_variance <- round(plsda_model$prop_expl_var$X[1:2] * 100)
        scores_df <- data.frame(
          Comp1 = scores[, 1],
          Comp2 = scores[, 2],
          Group = factor(Y) # Let ggplot handle colors by default
        )
        p <- ggplot(scores_df, aes(x = Comp1, y = Comp2, color = Group, fill = Group)) +
          geom_point(size = 3, shape = 21, stroke = 0.8) +
          stat_ellipse(geom = "polygon", alpha = 0.3, color = NA) +
          theme_minimal(base_family = "Arial") +
          labs(
            x = paste0("Component 1 (", explained_variance[1], "%)"),
            y = paste0("Component 2 (", explained_variance[2], "%)"),
            color = label,
            fill = label
          ) +
          theme(
            axis.title = element_text(size = 20, face = "bold"),
            axis.text = element_text(size = 18, color = "black"),
            legend.title = element_text(size = 18, face = "bold"),
            legend.text = element_text(size = 16),
            legend.position = c(0.95, 0.95),  # (x, y) from bottom-left corner
            legend.justification = c("right", "top"),
            panel.grid.major = element_line(color = "gray80", linewidth = 0.8),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
          )
        print(p)
        ggsave(
          filename = filename,
          plot = p,
          device = "svg",
          width = 8,
          height = 8,
          dpi = 600
        )
      }
    #- Filter feature tables to 0 v 2R with only no steroids, no washout, preop, or active steroids
      TFT_steroid_filtered_0v2 <- HILIC_C18_TFT_steroids %>%
        filter(!steroids %in% c("Preop","Washout","Yes")) %>%
        filter(ACR != "Preop") %>%
        filter(ACR != "1R") %>%
        select(ACR, everything(),-c(Sample_ID,steroids,Sample_no:paired_POD_EMB_day)) %>%
        mutate(ACR = droplevels(ACR))
      UFT_steroid_filtered_0v2 <- HILIC_C18_UFT_steroids %>%
        filter(!steroids %in% c("Preop", "Washout", "Yes")) %>%
        filter(ACR != "Preop") %>%
        filter(ACR != "1R") %>%
        select(ACR, everything(), -c(steroids, Sample_ID, Sample_no:paired_POD_EMB_day)) %>%
        mutate(ACR = droplevels(ACR))
    #- Run PLS-DA on UFT 0 vs 2R
      X <- UFT_steroid_filtered_0v2[, -1] # all features (remove the ACR column)
      Y <- UFT_steroid_filtered_0v2$ACR
      make_plsda_plot(X,Y, "ACR", "0v2_PLSDA_untargeted.svg")
    #- Run PLS-DA on TFT 0 vs 2R
      X2 <- TFT_steroid_filtered_0v2[, -1] # all features (remove the ACR column)
      Y2 <- TFT_steroid_filtered_0v2$ACR
      make_plsda_plot(X2,Y2, "ACR", "0v2_Targeted_PLSDA.svg")
    #- Export data to generate heatmap in MA
      UFT_steroid_filtered_0v2_MA <- HILIC_C18_UFT_steroids %>%
        filter(!steroids %in% c("Preop", "Washout", "Yes")) %>%
        filter(ACR != "Preop") %>%
        filter(ACR != "1R") %>%
        select(Sample_ID,ACR, everything(), -c(steroids, Sample_no:paired_POD_EMB_day)) %>%
        mutate(ACR = droplevels(ACR))
      write.csv(UFT_steroid_filtered_0v2_MA, "UFT_steroid_filtered_0v2_MA.csv", row.names = FALSE)
      TFT_steroid_filtered_0v2_MA <- HILIC_C18_TFT_steroids %>%
        filter(!steroids %in% c("Preop", "Washout", "Yes")) %>%
        filter(ACR != "Preop") %>%
        filter(ACR != "1R") %>%
        select(Sample_ID, ACR, everything(), -c(steroids, Sample_no:paired_POD_EMB_day)) %>%
        mutate(ACR = droplevels(ACR))
      # Remove duplicated names from targeted/isomers
      is_duplicate <- base::duplicated(as.list(TFT_steroid_filtered_0v2_MA), fromLast = FALSE)
      UFT_unique <- TFT_steroid_filtered_0v2_MA[, !is_duplicate]
      write.csv(UFT_unique, "TFT_steroid_filtered_0v2_MA.csv", row.names = FALSE)
  #+ Pathway enrichment
    #- Ttest for mummichog
      UFT_mummichog <- UFT_steroid_filtered_0v2 %>%
        filter(ACR %in% c("0R", "2R+")) %>%
        select(ACR, where(is.numeric)) %>%
        pivot_longer(-ACR, names_to = "Feature", values_to = "Value") %>%
        group_by(Feature) %>%
        filter(n_distinct(Value[ACR == "0R"]) > 1, n_distinct(Value[ACR == "2R+"]) > 1) %>%
        summarise(p.value = tryCatch(t.test(Value ~ ACR)$p.value, error = function(e) NA_real_), .groups = "drop") %>%
        arrange(p.value) %>%
        separate(Feature, into = c("mode", "m.z", "r.t"), sep = "_") %>%
        mutate(
          mode = case_when(
            mode == "HILIC" ~ "positive",
            mode == "C18" ~ "negative",
            TRUE ~ mode
          ),
          m.z = as.numeric(m.z),
          r.t = as.numeric(r.t)
        ) %>%
        select(m.z,mode,p.value,r.t)
      write.csv(UFT_mummichog, "UFT_mummichog.csv", row.names = FALSE)
    #- Create balloon plot on enrichment data
      pathway_enrich_results <- read_csv("Mummichog/primary_kegg.csv") %>%
        mutate(enrichment_factor = Hits.sig/Expected) %>%
        select(pathway_name, p_gamma, enrichment_factor) %>%
        filter(p_gamma < 0.05) %>%
        arrange(desc(enrichment_factor), p_gamma)
      ggplot(pathway_enrich_results, aes(
        x = 1, y = reorder(pathway_name, enrichment_factor),
        size = enrichment_factor, color = p_gamma
      )) +
        geom_point(alpha = 0.8) + # Bubbles with some transparency
        scale_size_continuous(
          range = c(3, 15), name = "Enrichment Factor"
        ) +
        guides(size = guide_legend(reverse = TRUE)) + # Reverse the legend
        scale_color_gradient(
          low = "#800017", high = "#EFD8DC", name = "P-Value"
        ) +
        theme_minimal(base_family = "Arial") +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 16, face = "bold", color = "black"), # Ensures pure black Y-axis text
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 16),
          legend.key.height = unit(1.5, "cm"), # Adjust space between keys to center title vertically
          legend.text.align = 0.5 # Center-align legend text horizontally
        )
#* 0R vs 1R vs 2R
  #- Filter and structure UFT
    UFT_steroid_filtered_0v1v2 <- HILIC_C18_UFT_steroids %>%
      filter(!steroids %in% c("Preop", "Washout", "Yes")) %>%
      filter(ACR != "Preop") %>%
      select(ACR, everything(), -c(steroids, Sample_ID, Sample_no:paired_POD_EMB_day)) %>%
      mutate(ACR = droplevels(ACR))
    mod_ACR <- read_excel("/Users/jdp2019/Library/CloudStorage/OneDrive-Emory/Research/Manuscripts and Projects/Active Projects/Chan Lab/TPMO/Metabolomics - ACR_Clay/study_metadata.xlsx") %>%
      select(original_sample_ID, mod_ACR) %>%
      rename(Sample_ID = original_sample_ID)
    UFT_steroid_filtered_0vprev2 <- HILIC_C18_UFT_steroids %>%
      left_join(mod_ACR, by = "Sample_ID") %>%
      filter(!steroids %in% c("Preop", "Washout", "Yes")) %>%
      filter(ACR != "Preop") %>%
      select(mod_ACR, everything(), -c(ACR,steroids, Sample_ID, Sample_no:paired_POD_EMB_day))
  #- Run PLS-DA on UFT 0 vs 1 vs 2R
    X1 <- UFT_steroid_filtered_0v1v2[, -1] # all features (remove the ACR column)
    Y1 <- UFT_steroid_filtered_0v1v2$ACR
    make_plsda_plot(X1, Y1, "ACR", "0v1v2_PLSDA.svg")
    X2 <- UFT_steroid_filtered_0vprev2[, -1] # all features (remove the ACR column)
    Y2 <- UFT_steroid_filtered_0vprev2$mod_ACR
    make_plsda_plot(X2, Y2, "mod_ACR", "0vprev2_PLSDA.svg")
    X3 <- UFT_steroid_filtered_0vprevcombo_2[, -1] # all features (remove the ACR column)
    Y3 <- UFT_steroid_filtered_0vprevcombo_2$mod_ACR
    make_plsda_plot(X3, Y3, "mod_ACR", "0vprev2_PLSDA.svg")
    X4 <- UFT_steroid_filtered_mod_ACR_pp[, -1] # all features (remove the ACR column)
    Y4 <- UFT_steroid_filtered_mod_ACR_pp$mod_ACR_pp
    make_plsda_plot(X4, Y4, "mod_ACR_pp", "0vprev2_PLSDA.svg")
#* within subject replicates
#* LASSO
#+ Clean and prep data
  tft_lasso <- TFT_steroid_filtered_0v2 %>%
    filter(ACR %in% c("0R", "2R+")) %>%
    mutate(ACR_bin = ifelse(ACR == "2R+", 1, 0)) %>%
    {
      # Step 1: Keep metadata separate
      meta_cols <- c("ACR_bin", "ACR")
      data_only <- select(., -all_of(meta_cols))

      # Step 2: Remove duplicated columns (by content)
      data_matrix <- as.data.frame(data_only)
      data_unique <- data_matrix[, !base::duplicated(as.data.frame(t(data_matrix)))]

      # Step 3: Remove low-information columns
      good_metabolites <- names(data_unique)[!sapply(data_unique, function(col) {
        min_val <- min(col, na.rm = TRUE)
        mean(col == min_val, na.rm = TRUE) > 0.25
      })]

      # Step 4: Recombine
      bind_cols(select(., all_of(meta_cols)), data_unique[good_metabolites])
    }
#+ LASSO and graph selected features
  X <- tft_lasso %>% select(-ACR, -ACR_bin) %>% as.matrix()
  y <- tft_lasso$ACR_bin
  set.seed(123)
  cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1, standardize = TRUE)
  best_lambda <- cv_fit$lambda.min
  lasso_model <- glmnet(X, y, family = "binomial", alpha = 1, lambda = best_lambda)
  coef_df <- as.matrix(coef(lasso_model))
  selected_features <- coef_df[coef_df[, 1] != 0, , drop = FALSE]
  selected_names <- rownames(selected_features)[-1]  # drop intercept
  pred_probs <- predict(lasso_model, newx = X, type = "response")[, 1]
  roc_obj <- pROC::roc(y, pred_probs)
  auc_val <- pROC::auc(roc_obj)
  coords_best <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  threshold <- as.numeric(coords_best["threshold"])
  pred_labels <- ifelse(as.vector(pred_probs) >= threshold, 1, 0)
  conf_mat <- table(Predicted = pred_labels, Actual = y)
  feature_labels <- c(
    "Phenethylamine_C05332_M+H_HILIC" = "Phenethylamine",
    "Cholest-5-ene-3beta,26-diol_C15610_M+H_HILIC" = "26-hydroxycholesterol",
    "4-O-beta-D-Glucopyranosyl-D-mannose_C02964_M-H_C18" = "Mannobiose"
  )
  top_features <- base::intersect(names(feature_labels), selected_names)
  top_feature_short_names <- unname(feature_labels[top_features])
  plot_data <- tft_lasso %>%
    select(ACR, all_of(top_features)) %>%
    pivot_longer(-ACR, names_to = "Feature", values_to = "Value") %>%
    mutate(Feature = recode(Feature, !!!feature_labels))
  plot_feature <- function(feature_name) {
    ggplot(filter(plot_data, Feature == feature_name), aes(x = ACR, y = Value)) +
      geom_boxplot(aes(fill = ACR), outlier.shape = NA, color = "black", alpha = 1) +
      geom_jitter(aes(fill = ACR), shape = 21, size = 2.5, stroke = 0.2, width = 0.2,
                  color = "black") +
      scale_fill_manual(values = c("0R" = "#b9e7e8", "2R+" = "#fccbc9")) +
      labs(title = feature_name, x = NULL, y = "Intensity") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      )
  }
  plots <- lapply(top_feature_short_names, plot_feature)

  # Arrange: 2 top, 1 centered bottom
  layout <- "
  AB
  .C
  "

  wrap_plots(A = plots[[1]], B = plots[[2]], C = plots[[3]], design = layout)

  #+ ROC curve
  roc_df <- data.frame(
    FPR = rev(roc_obj$specificities),
    TPR = rev(roc_obj$sensitivities)
  )
  auc_val <- pROC::auc(roc_obj)
  coords_best <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  sensitivity <- round(coords_best["sensitivity"] * 100)
  specificity <- round(coords_best["specificity"] * 100)
  auc_label <- paste0("AUC = ", round(auc_val, 2))
  metrics_label <- paste0("Sensitivity = ", sensitivity, "%\nSpecificity = ", specificity, "%")
  # Combined annotation
  annotation_text <- paste(auc_label, metrics_label, sep = "\n")
  # Prism-style ROC plot
  ggplot(roc_df, aes(x = 1 - FPR, y = TPR)) +
    geom_step(color = "#f45e5a", size = 1.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    coord_equal() +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = "ROC Curve: 0R vs 2R+",
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    annotate("text", x = 0.97, y = 0.03, hjust = 1, vjust = 0,
            label = annotation_text, size = 10, fontface = "bold", lineheight = 1.15) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 16),
      panel.grid = element_blank()
    )


    pred_probs <- predict(lasso_model, newx = X, type = "response")[,1]
    roc_obj <- pROC::roc(y, pred_probs)
    auc_val <- pROC::auc(roc_obj)

    # Best threshold for binary classification
    coords_best <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
    threshold <- as.numeric(coords_best["threshold"])
    sensitivity <- round(coords_best["sensitivity"], 2)
    specificity <- round(coords_best["specificity"], 2)

    # Coefficients (excluding intercept)
    coef_vec <- coef(lasso_model)
    coef_df <- as.data.frame(as.matrix(coef_vec))
    coef_df <- coef_df[-1, , drop = FALSE]  # drop intercept
    coef_df <- coef_df[coef_df[,1] != 0, , drop = FALSE]  # only non-zero
    names(coef_df) <- "Coefficient"

    # Rename features for display
    coef_df$Feature <- rownames(coef_df)
    coef_df$Feature <- feature_labels[coef_df$Feature]

    roc_df <- data.frame(
      FPR = rev(roc_obj$specificities),
      TPR = rev(roc_obj$sensitivities)
    )

    # AUC + Metrics
    auc_val <- pROC::auc(roc_obj)
    coords_best <- pROC::coords(roc_obj, "best", ret = c("sensitivity", "specificity"), best.method = "youden")
    sensitivity <- round(coords_best["sensitivity"] * 100)
    specificity <- round(coords_best["specificity"] * 100)

    # Feature names
    feature_names <- feature_labels[top4_features]

    # Expression-based annotations
    feature_exprs <- lapply(feature_names, function(name) bquote(bolditalic(.(name))))
    annotation_lines <- c(
      bquote(bold("AUC = ") * .(round(auc_val, 2))),
      bquote(bold("Sensitivity = ") * .(sensitivity) * "%"),
      bquote(bold("Specificity = ") * .(specificity) * "%"),
      bquote(""),  # line break
      bquote(bold(underline("Selected features"))),
      feature_exprs
    )

    # Create annotate layers
    anno_layers <- Map(function(txt, i) {
      annotate("text", x = 0.97, y = 0.03 + 0.045 * i,
              label = txt, parse = TRUE, hjust = 1, vjust = 0, size = 5.5)
    }, annotation_lines, seq_along(annotation_lines))

    library(purrr)

    # Start with the base plot
    base_plot <- ggplot(roc_df, aes(x = 1 - FPR, y = TPR)) +
      geom_step(color = "#f45e5a", size = 1.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      coord_equal() +
      scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(
        title = "ROC Curve: 0R vs 2R+",
        x = "1 - Specificity",
        y = "Sensitivity"
      ) +
      theme_minimal(base_size = 20) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 16),
        panel.grid = element_blank()
      )

    # Add annotation layers using reduce
    final_plot <- reduce(anno_layers, `+`, .init = base_plot)

    # Display it
    final_plot
#* RM graphs of interesting features
  #+ Untargeted
    HILIC_734 <- HILIC_C18_UFT_steroids %>%
      select(Patient_no,ACR,Sample_no,paired_POD_EMB_day,"HILIC_734.7471_237.7") %>%
      mutate(z_HILIC_734 = as.numeric(scale(HILIC_734.7471_237.7)))
      acr_colors <- c(
        "Preop" = "grey60",
        "0R" = "lightgreen",
        "1R" = "goldenrod",
        "2R+" = "darkred",
        "pre-2R" = "orange"
      )
      ggplot(HILIC_734, aes(x = Sample_no, y = z_HILIC_734)) +
        geom_point(aes(color = ACR), size = 2.5, alpha = 0.8) +
        geom_smooth(aes(group = Patient_no), method = "loess", se = FALSE, color = "gray", linewidth = 0.8) +
        scale_color_manual(values = acr_colors) +
        facet_wrap(~Patient_no, scales = "free_x") +
        labs(
          x = "Days Post-Transplant",
          y = "Intensity (m/z 734.7471, RT 237.7)"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 12)
        )
  #+ Targeted
  myco <- read_excel("mycophenolate.xlsx", sheet = "MCP") %>%
    filter(ACR!= "Preop")
 rm_plot <- function(data, metabolite) {
  ggplot(data, aes(x = Sample_no, y = {{ metabolite }})) +
    geom_point(aes(color = ACR), size = 2.5, alpha = 0.8) +
    geom_smooth(aes(group = Patient_no), method = "loess", se = FALSE, color = "gray", linewidth = 0.8) +
    scale_color_manual(values = acr_colors) +
    facet_wrap(~Patient_no, scales = "free_x") +
    labs(
      x = "Sample Number",
      y = rlang::as_label(rlang::enquo(metabolite)),
      color = "ACR Grade"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_text(face = "bold"),
      legend.position = "right"  # make sure it's shown
    )
  }
  rm_plot(myco, Mycophenolate)
  # Rename using short names for all four features
  targeted_interest <- HILIC_C18_TFT_steroids %>%
    # filter(!steroids %in% c("Preop", "Washout", "Yes")) %>%
    # filter(!ACR %in% c("Preop")) %>%
    mutate(across("UDP_C00015_M+2Na-H_HILIC":"2-Butyloctyl_sulfate_C22403_M+FA-H_C18", ~ scale(.)[, 1])) %>%
    rename(
      phenethylamine = "Phenethylamine_C05332_M+H_HILIC",
      cholest_5_ene_3beta = "Cholest-5-ene-3beta,26-diol_C15610_M+H_HILIC",
      glucopyranosyl_manose = "4-O-beta-D-Glucopyranosyl-D-mannose_C02964_M-H_C18",
      mannobiose = "Mannobiose_C20861_M-H_C18"
    ) %>%
    select(Patient_no, ACR, Sample_no, paired_POD_EMB_day, everything())

    rm_plot(targeted_interest, phenethylamine)
  rm_plot(targeted_interest, cholest_5_ene_3beta)

