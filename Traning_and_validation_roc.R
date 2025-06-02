# R Script: Lung Cancer Immune Therapy Data Analysis
# Author: Shaoqiu Chen
# Date: 2025-06-02 
# Description: This script loads pre-processed data, required packages,
#              and then performs training and validation analyses,
#              including ROC curves, forest plots, and violin plots.

# --- 0. Setup and Recommendations ---
#
# IMPORTANT FOR GITHUB:
# 1. Working Directory:
#    - This script AVOIDS using `setwd()` with absolute paths.
#    - BEST PRACTICE: Use RStudio Projects. Create a .Rproj file in your
#      main project folder. When you open the .Rproj file, R's working
#      directory will automatically be set to that folder.
#    - All file paths below (for loading data and saving plots) should then
#      be relative to this project root.
#
# 2. Data File ("Science_advance.RData"):
#    - Place this file in your project's root directory or in a subdirectory
#      (e.g., "data/Science_advance.RData"). Adjust the `data_file` path below
#      accordingly.
#    - If this file is too large for GitHub, host it elsewhere and provide
#      download instructions in your README.md.
#
# 3. Output Directories:
#    - This script will save plots. It's good practice to save them into
#      subdirectories (e.g., "results/plots/").
#    - You might need to create these directories manually or add R code
#      to create them (e.g., using `dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)`).


###download Science_advance.RData at link:



# --- 1. Load Pre-processed Data/Environment ---


data_file <- "Science_advance.RData" # Or "data/Science_advance.RData" if in a subfolder
if (file.exists(data_file)) {
  load(data_file)
  print(paste("Successfully loaded:", data_file))
  # Optional: List loaded objects to verify
  # print("Objects loaded:")
  # print(ls())
} else {
  warning(paste("Data file not found:", data_file,
                "- Please ensure it's in the correct path relative to your project root."))
}

# --- 2. Create Output Directories (Recommended) ---
# Create directories to store the generated plots if they don't exist.
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/training_plots")) dir.create("results/training_plots")
if (!dir.exists("results/validation_plots")) dir.create("results/validation_plots")
if (!dir.exists("results/forest_plots")) dir.create("results/forest_plots")

# --- 3. Load Required R Packages ---
print("Loading required packages...")
# For meta-analysis
library(MetaIntegrator)

# For data manipulation and general utilities
library(dplyr)    # Data manipulation
library(tidyr)    # Tidying data
library(tibble)   # Modern data frames

# For plotting and visualization
library(ggplot2)  # Creating graphics
library(ggthemes) # Additional themes for ggplot2
library(ggpubr)   # Enhancements for ggplot2
library(survminer) # Survival analysis plots

# For statistical analysis and modeling
library(pROC)     # ROC curve analysis
library(multiROC) # Multi-class ROC curve analysis
library(rstatix)  # Pipe-friendly statistical tests
library(survival) # Core survival analysis functions

print("All required packages have been loaded (or attempted to load).")

# --- 4. Define Reusable Plotting Theme (Optional) ---
# You can define a common theme to make your plots consistent.
common_theme <- theme_bw() + # Or another base theme like theme_classic()
  theme(
    legend.title = element_blank(),
    panel.background = element_blank(), # Often handled by theme_bw/classic
    plot.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5) # ggplot2 uses linewidth
  )

# --- 5. Training Analysis ---
print("Starting Training Analysis...")

# Calculate training ROC (assuming discovery_exampleMetaObj and forwardRes are loaded)
if (exists("discovery_exampleMetaObj") && exists("forwardRes")) {
  training_roc_plot_obj <- summaryROCPlot(metaObject = discovery_exampleMetaObj,
                                          filterObject = forwardRes,
                                          bootstrapReps = 100)
  
  training_roc_plot_obj <- training_roc_plot_obj +
    common_theme # Apply common theme
  
  ggsave("results/training_plots/training_roc_curve.pdf", plot = training_roc_plot_obj, width = 6.73, height = 6.83)
  print("Saved training ROC curve.")
} else {
  warning("Skipping training ROC plot: 'discovery_exampleMetaObj' or 'forwardRes' not found.")
}

# Violin Plots for Training Datasets
# GSE135222
if (exists("forwardRes") && exists("GSE135222")) {
  gse135222_violin_plot <- violinPlot(forwardRes, GSE135222, labelColumn = 'group') +
    common_theme
  ggsave("results/training_plots/GSE135222_violinplot.pdf", plot = gse135222_violin_plot, width = 5.09, height = 6.02)
  print("Saved GSE135222 violin plot.")
} else {
  warning("Skipping GSE135222 violin plot: 'forwardRes' or 'GSE135222' not found.")
}

# POPLAR
if (exists("forwardRes") && exists("POPLAR")) {
  poplar_violin_plot <- violinPlot(forwardRes, POPLAR, labelColumn = 'group') +
    common_theme
  ggsave("results/training_plots/POPLAR_violinplot.pdf", plot = poplar_violin_plot, width = 5.09, height = 6.02)
  print("Saved POPLAR violin plot.")
} else {
  warning("Skipping POPLAR violin plot: 'forwardRes' or 'POPLAR' not found.")
}

# OKA (different type)
if (exists("forwardRes") && exists("OKA")) {
  oka_hist_violin_plot <- violinPlot(forwardRes, OKA, labelColumn = 'HIST') +
    common_theme
  ggsave("results/training_plots/OKA_different_type_violinplot.pdf", plot = oka_hist_violin_plot, width = 5.09, height = 6.02)
  print("Saved OKA different type violin plot.")
} else {
  warning("Skipping OKA (HIST) violin plot: 'forwardRes' or 'OKA' not found.")
}

# --- 6. Validation Analysis ---
print("Starting Validation Analysis...")

# Calculate validation ROC (assuming validation_exampleMetaObj and forwardRes are loaded)
if (exists("validation_exampleMetaObj") && exists("forwardRes")) {
  validation_roc_plot_obj <- summaryROCPlot(metaObject = validation_exampleMetaObj,
                                            filterObject = forwardRes,
                                            bootstrapReps = 100)
  
  validation_roc_plot_obj <- validation_roc_plot_obj +
    common_theme # Apply common theme
  
  ggsave("results/validation_plots/validation_roc_curve.pdf", plot = validation_roc_plot_obj, width = 6.73, height = 6.83)
  print("Saved validation ROC curve.")
} else {
  warning("Skipping validation ROC plot: 'validation_exampleMetaObj' or 'forwardRes' not found.")
}

# Violin Plots for Validation Datasets
# OKA (group)
if (exists("forwardRes") && exists("OKA")) {
  oka_group_violin_plot <- violinPlot(forwardRes, OKA, labelColumn = 'group') +
    common_theme
  ggsave("results/validation_plots/OKA_violinplot.pdf", plot = oka_group_violin_plot, width = 5.09, height = 6.02)
  print("Saved OKA (group) violin plot.")
} else {
  warning("Skipping OKA (group) violin plot: 'forwardRes' or 'OKA' not found.")
}

# GSE111414
if (exists("forwardRes") && exists("GSE111414")) {
  gse111414_violin_plot <- violinPlot(forwardRes, GSE111414, labelColumn = 'group') +
    common_theme
  ggsave("results/validation_plots/GSE111414_violinplot.pdf", plot = gse111414_violin_plot, width = 5.09, height = 6.02)
  print("Saved GSE111414 violin plot.")
} else {
  warning("Skipping GSE111414 violin plot: 'forwardRes' or 'GSE111414' not found.")
}

# WHTJ2 Violin Plot
if (exists("forwardRes") && exists("WHTJ2")) {
  whtj2_violin_plot <- violinPlot(forwardRes, WHTJ2, labelColumn = 'group') +
    common_theme
  ggsave("results/validation_plots/WHTJ2_violinplot.pdf", plot = whtj2_violin_plot, width = 5.09, height = 6.02) # Added ggsave for this
  print("Saved WHTJ2 violin plot.")
  
  # WHTJ2 ROC Plot
  whtj2_roc_plot_obj <- rocPlot(forwardRes, WHTJ2, title = "ROC plot for PSC Validation WHTJ2") + # Added WHTJ2 to title
    common_theme
  ggsave("results/validation_plots/WHTJ2_roc_curve.pdf", plot = whtj2_roc_plot_obj, width = 6.73, height = 6.83) # Added ggsave
  print("Saved WHTJ2 ROC curve.")
  
} else {
  warning("Skipping WHTJ2 plots: 'forwardRes' or 'WHTJ2' not found.")
}


# --- 7. Forest Plots for Genes ---
# This section generates forest plots for genes identified in forwardRes.
# It uses 'exampleMetaObj'. If you have separate metaObjects for discovery/validation
# and want different forest plots, you might need to adapt this section or run it twice
# with the appropriate metaObject.

if (exists("forwardRes") && exists("exampleMetaObj")) {
  all_genes <- c(forwardRes$posGeneNames, forwardRes$negGeneNames)
  
  if (length(all_genes) > 0) {
    print(paste("Generating forest plots for", length(all_genes), "genes..."))
    # Note: The original forestPlot function from MetaIntegrator might draw directly
    # or return a grid object. The pdf()/print()/dev.off() pattern is typical.
    # `graphics.off()` is generally not needed if `dev.off()` is used correctly for each pdf.
    for (gene in all_genes) {
      file_name <- file.path("results", "forest_plots", paste0(gsub("[^A-Za-z0-9_.-]", "_", gene), "_forest_plot.pdf")) # Sanitize gene name for filename
      pdf(file_name, width = 4.48, height = 4.85)
      tryCatch({
        # forestPlot in MetaIntegrator might not return a ggplot object,
        # it might plot directly. Ensure it works with pdf redirection.
        # If it's a base plot, just calling the function is enough.
        # If it returns a grid object (like some complex plots), it might need `grid.draw()`.
        # For now, assuming it plots directly or print() works for its output type.
        current_plot <- forestPlot(
          geneName = gene,
          metaObject = exampleMetaObj,
          boxColor = "black",
          whiskerColor = "black",
          zeroLineColor = "black",
          summaryColor = "black",
          textColor = "black"
        )
        if (!is.null(current_plot)) { # If it returns an object that can be printed (like ggplot or grid)
          print(current_plot)
        }
        print(paste("Saved forest plot for gene:", gene, "to", file_name))
      }, error = function(e) {
        warning(paste("Error generating forest plot for gene:", gene, "-", e$message))
      }, finally = {
        dev.off() # Ensure device is closed even if there's an error
      })
    }
    print("Finished generating forest plots.")
  } else {
    print("No genes found in forwardRes to generate forest plots.")
  }
} else {
  warning("Skipping forest plots: 'forwardRes' or 'exampleMetaObj' not found.")
}

# --- End of Script ---
print("Script execution finished.")

