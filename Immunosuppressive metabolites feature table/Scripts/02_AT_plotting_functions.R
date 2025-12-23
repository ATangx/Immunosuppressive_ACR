# Plotting functions for metabolomics analysis

#' Create bar plots comparing metabolite levels between ACR groups
#' @param sig_results Dataframe with significant metabolite results
#' @param group1_data First group data (0R/1R)
#' @param group2_data Second group data (2R+)
#' @param use_median Logical, whether to use median instead of mean
plot_ttest_bars <- function(sig_results, group1_data, group2_data, use_median = FALSE) {
  
  plot_list <- list()
  
  for (i in 1:nrow(sig_results)) {
    metabolite <- sig_results$Metabolite[i]
    
    # Get the metabolite data
    met_col <- which(colnames(group1_data) == metabolite)
    if (length(met_col) == 0) next
    
    # Prepare data for plotting
    group1_values <- group1_data[[met_col]]
    group2_values <- group2_data[[met_col]]
    
    # Create combined dataframe for ggplot
    plot_data <- data.frame(
      values = c(group1_values, group2_values),
      group = c(rep("0R/1R", length(group1_values)), 
                rep("2R+", length(group2_values)))
    )
    
    # Remove NA values
    plot_data <- plot_data[!is.na(plot_data$values), ]
    
    # Calculate summary statistics
    if (use_median) {
      summary_data <- plot_data %>%
        group_by(group) %>%
        summarise(
          center = median(values, na.rm = TRUE),
          .groups = 'drop'
        )
      y_label <- paste("Median", metabolite, "Level")
    } else {
      summary_data <- plot_data %>%
        group_by(group) %>%
        summarise(
          center = mean(values, na.rm = TRUE),
          se = sd(values, na.rm = TRUE) / sqrt(n()),
          .groups = 'drop'
        )
      y_label <- paste("Mean", metabolite, "Level")
    }
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = group, y = values, fill = group)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      stat_summary(fun = ifelse(use_median, median, mean), 
                   geom = "point", 
                   size = 3, 
                   color = "red") +
      labs(
        title = paste("Comparison of", metabolite, "levels"),
        subtitle = paste("p-value =", round(sig_results$p_value[i], 4)),
        x = "ACR Group",
        y = y_label
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none"
      ) +
      scale_fill_manual(values = c("0R/1R" = "lightblue", "2R+" = "lightcoral"))
    
    # Store the plot
    plot_list[[metabolite]] <- p
  }
  
  return(plot_list)
}