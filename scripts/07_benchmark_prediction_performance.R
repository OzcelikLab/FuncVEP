library(data.table)
library(pROC)
library(glue)
library(ggplot2)
library(gridExtra)
library(grid)
library(patchwork)
library(RColorBrewer)
library(showtext)

merge_model_predictions <- function(model_name, benchmarking_data, result_dir) {
  prediction_path <- file.path(result_dir, paste0(model_name, ".txt"))
  
  predictions <- fread(prediction_path, sep = "\t")
  
  if (!("ID" %in% colnames(predictions)) || !(model_name %in% colnames(predictions))) {
    stop(glue("Missing 'ID' or '{model_name}' column in {prediction_path}"))
  }
  predictions <- predictions[, .(ID, get(model_name))]
  setnames(predictions, "V2", model_name)
  
  merged_data <- merge(benchmarking_data, predictions, by = "ID", all.x = TRUE)
  return(merged_data)
}

benchmark_tool_performance <- function(data, label_column, id_column = "ID") {
  tool_columns <- setdiff(colnames(data), c(id_column, label_column))
  
  data[, (tool_columns) := lapply(.SD, as.numeric), .SDcols = tool_columns]
  
  tool_columns <- tool_columns[sapply(tool_columns, function(tool) {
    sum(!is.na(data[get(label_column) == 1, get(tool)])) > 0 &&
      sum(!is.na(data[get(label_column) == 0, get(tool)])) > 0
  })]
  
  results <- data.table(
    Tool = character(),
    Accuracy = numeric(),
    AUC = numeric(),
    OptimalCutoff = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    YoudensJ = numeric(),
    HigherIndicatesPathogenicity = logical()
  )
  
  for (tool in tool_columns) {
    roc_obj <- pROC::roc(data[[label_column]], data[[tool]], quiet = TRUE)
    coords <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
    
    auc_value <- pROC::auc(roc_obj)
    youdens_j <- coords$sensitivity + coords$specificity - 1
    threshold <- as.numeric(coords$threshold)
    
    median_path <- median(data[get(label_column) == 1, get(tool)], na.rm = TRUE)
    median_benign <- median(data[get(label_column) == 0, get(tool)], na.rm = TRUE)
    higher_indicates_pathogenicity <- median_path > median_benign
    
    if (!is.na(threshold)) {
      predictions <- ifelse(data[[tool]] >= threshold, 1, 0)
      valid <- !is.na(predictions) & !is.na(data[[label_column]])
      accuracy <- sum(predictions[valid] == data[[label_column]][valid]) / sum(valid)
      
      if (!higher_indicates_pathogenicity) {
        accuracy <- 1 - accuracy
      }
    } else {
      accuracy <- NA
    }
    
    results <- rbind(
      results,
      data.table(
        Tool = tool,
        Accuracy = accuracy,
        AUC = as.numeric(auc_value),
        OptimalCutoff = threshold,
        Sensitivity = as.numeric(coords$sensitivity),
        Specificity = as.numeric(coords$specificity),
        YoudensJ = youdens_j,
        HigherIndicatesPathogenicity = higher_indicates_pathogenicity
      ),
      fill = TRUE
    )
    
    plot(roc_obj, main = paste("ROC Curve -", tool))
  }
  
  return(results[, .(Tool, Accuracy, AUC, OptimalCutoff, Sensitivity, Specificity, YoudensJ, HigherIndicatesPathogenicity)])
}

# Benchmark for functional dataset
variant_features <- fread("../data/intermediate/variant_features_core.txt", sep = "\t")
variant_labels <- fread("../data/intermediate/variant_labels.txt", sep = "\t")
# primateai <- fread("../data/restricted/primateai.txt", sep = "\t")

# variant_features[primateai, on=.(ID), `:=` (PrimateAI_score=PrimateAI_score, PAI3D_score=PAI3D_score)]

variant_labels[, functional_label := fifelse(functional_label == "PS3", 1,
                                             fifelse(functional_label == "BS3", 0, NA_integer_))]
variant_labels[, clinical_label := fifelse(clinical_label == "P", 1,
                                           fifelse(clinical_label == "B", 0, NA_integer_))]

all_columns <- fread("../resources/feature_lists/all_columns.txt", sep = "\t")
prediction_tools <- all_columns[Type == "Variant Effect Predictor" & Name %in% colnames(variant_features), Name]

benchmarking_data_functional <- variant_features[, c("ID", prediction_tools), with = FALSE]
benchmarking_data_functional <- merge(benchmarking_data_functional, variant_labels[, .(ID, functional_label)], by = "ID", all.x = TRUE)
benchmarking_data_functional <- merge(benchmarking_data_functional, variant_labels[, .(ID, weight)], by = "ID", all.x = TRUE)
benchmarking_data_functional <- benchmarking_data_functional[weight == 1 & !is.na(functional_label)]
benchmarking_data_functional[, weight := NULL]  # optional cleanup
benchmarking_data_functional <- unique(benchmarking_data_functional, by = "ID")


models <- c("FuncVEP_CTI", "FuncVEP_CTE", "FuncVEP_SP", "ClinVEP_CTI", "ClinVEP_CTE", "ClinVEP_SP")
for (model in models) {
  benchmarking_data_functional <- merge_model_predictions(model, benchmarking_data_functional, result_dir = "../results/predictions/functional")
  benchmarking_data_functional <- unique(benchmarking_data_functional, by = "ID")
}

results_functional <- benchmark_tool_performance(benchmarking_data_functional, label_column = "functional_label", id_column = "ID")

no_training_set <- fread("../resources/feature_lists/tools_excluded_due_to_unavailable_training_sets.txt", header = FALSE)[[1]]
results_functional_filtered <- results_functional[!Tool %in% no_training_set]

fwrite(results_functional_filtered, "../results/benchmarking/functional.txt", sep = "\t")


# Benchmark for clinical dataset
benchmarking_data_clinical <- variant_features[, c("ID", prediction_tools), with = FALSE]
benchmarking_data_clinical <- merge(benchmarking_data_clinical, variant_labels[, .(ID, clinical_label)], by = "ID", all.x = TRUE)
benchmarking_data_clinical <- benchmarking_data_clinical[!is.na(clinical_label) & clinical_label != ""]

for (model in models) {
  benchmarking_data_clinical <- merge_model_predictions(model, benchmarking_data_clinical, result_dir = "../results/predictions/clinical")
  benchmarking_data_clinical <- unique(benchmarking_data_clinical, by = "ID")
}

results_clinical <- benchmark_tool_performance(benchmarking_data_clinical, label_column = "clinical_label", id_column = "ID")

results_clinical_filtered <- results_clinical[!Tool %in% no_training_set]
fwrite(results_clinical_filtered, "../results/benchmarking/clinical.txt", sep = "\t")



benchmark_tool_performance <- function(data, label_column, id_columns = c("chr", "pos", "ref", "alt", "ensg")) {
  tool_columns <- setdiff(colnames(data), c(id_columns, label_column))
  
  data[, (tool_columns) := lapply(.SD, as.numeric), .SDcols = tool_columns]
  
  tool_columns <- tool_columns[sapply(tool_columns, function(tool) {
    sum(!is.na(data[get(label_column) == 1, get(tool)])) > 0 &&
      sum(!is.na(data[get(label_column) == 0, get(tool)])) > 0
  })]
  
  results <- data.table()
  
  for (tool in tool_columns) {
    roc_obj <- roc(data[[label_column]], data[[tool]], quiet = TRUE)
    coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
    
    auc_value <- auc(roc_obj)
    youdens_j <- coords$sensitivity + coords$specificity - 1
    threshold <- as.numeric(coords$threshold)
    
    median_path <- median(data[get(label_column) == 1, get(tool)], na.rm = TRUE)
    median_benign <- median(data[get(label_column) == 0, get(tool)], na.rm = TRUE)
    higher_indicates_pathogenicity <- median_path > median_benign
    
    predictions <- ifelse(data[[tool]] >= threshold, 1, 0)
    valid <- !is.na(predictions) & !is.na(data[[label_column]])
    accuracy <- sum(predictions[valid] == data[[label_column]][valid]) / sum(valid)
    
    if (!higher_indicates_pathogenicity) {
      accuracy <- 1 - accuracy
    }
    
    results <- rbind(
      results,
      data.table(
        Tool = tool,
        Accuracy = accuracy,
        AUC = as.numeric(auc_value),
        OptimalCutoff = threshold,
        Sensitivity = as.numeric(coords$sensitivity),
        Specificity = as.numeric(coords$specificity),
        YoudensJ = youdens_j,
        HigherIndicatesPathogenicity = higher_indicates_pathogenicity
      ),
      fill = TRUE
    )
  }
  
  return(results)
}


mart <- fread("../resources/gene_ids/mart_export_02_april_2025.txt", sep="\t")
setnames(mart, old = "Gene", new = "gene")
mart[, enst := NULL]
mart <- mart[!duplicated(gene) & !duplicated(gene, fromLast = TRUE)]

no_training_set <- fread("../resources/feature_lists/tools_excluded_due_to_unavailable_training_sets.txt", header = FALSE)[[1]]

test_variants <- fread("../data/intermediate/denovo_test_variants.txt", sep = "\t")

all_columns <- fread("../resources/feature_lists/all_columns.txt", sep = "\t")
prediction_tools <- all_columns[Type == "Variant Effect Predictor" & Name %in% colnames(test_variants), Name]

id_cols <- c("chr", "pos", "ref", "alt", "ensg")

# Benchmark for DD
dd_dataset <- fread("../data/datasets/testing/de_novo/DD_dataset.txt", sep = "\t")
dd_dataset <- dd_dataset[phenotype %in% c("DD", "DD_control")]
dd_dataset[, label := fifelse(phenotype == "DD", 1, 0)]
dd_dataset[mart, on=.(gene), ensg:=i.ensg]
dd_dataset <- dd_dataset[!is.na(ensg)]

dd_with_scores <- merge(dd_dataset[, c(id_cols, "label"), with = FALSE], test_variants, by = id_cols, all.x = TRUE)
dd_results <- benchmark_tool_performance(dd_with_scores, label_column = "label", id_columns = id_cols)
dd_results <- dd_results[!Tool %in% no_training_set]
fwrite(dd_results, "../results/benchmarking/dd.txt", sep = "\t")

# Benchmark for NDD
ndd_dataset <- fread("../data/datasets/testing/de_novo/NDD_dataset.txt", sep = "\t")
ndd_dataset <- ndd_dataset[phenotype %in% c("NDD", "NDD_control")]
ndd_dataset[, label := fifelse(phenotype == "NDD", 1, 0)]
ndd_dataset[mart, on=.(gene), ensg:=i.ensg]
ndd_dataset <- ndd_dataset[!is.na(ensg)]

ndd_with_scores <- merge(ndd_dataset[, c(id_cols, "label"), with = FALSE], test_variants, by = id_cols, all.x = TRUE)
ndd_results <- benchmark_tool_performance(ndd_with_scores, label_column = "label", id_columns = id_cols)
ndd_results <- ndd_results[!Tool %in% no_training_set]
fwrite(ndd_results, "../results/benchmarking/ndd.txt", sep = "\t")

# Benchmark for cancer
cancer_dataset <- fread("../data/datasets/testing/cancer/cancer_hotspot_dataset.txt", sep = "\t")
cancer_dataset <- cancer_dataset[!is.na(hotspot)]
cancer_dataset[, label := as.integer(hotspot)]
cancer_dataset[mart, on=.(gene), ensg:=i.ensg]
cancer_dataset <- cancer_dataset[!is.na(ensg)]

cancer_with_scores <- merge(cancer_dataset[, c(id_cols, "label"), with = FALSE], test_variants, by = id_cols, all.x = TRUE)
cancer_results <- benchmark_tool_performance(cancer_with_scores, label_column = "label", id_columns = id_cols)
cancer_results <- cancer_results[!Tool %in% no_training_set]
fwrite(cancer_results, "../results/benchmarking/cancer.txt", sep = "\t")



save_combined_functional_clinical_roc <- function(data_functional,
                                                  data_clinical,
                                                  tool_names,
                                                  label_column_functional,
                                                  label_column_clinical) {
  palette_colors <- brewer.pal(min(length(tool_names), 8), "Set1")
  
  build_roc_plot <- function(data, label_column, panel_title) {
    roc_dt <- rbindlist(lapply(tool_names, function(tool) {
      display_name <- gsub("_", "-", tool)
      d <- data[!is.na(get(tool)) & !is.na(get(label_column))]
      roc_obj <- pROC::roc(d[[label_column]], d[[tool]], quiet = TRUE)
      auc_val <- as.numeric(pROC::auc(roc_obj))
      
      coords <- as.data.table(
        pROC::coords(roc_obj, "all", ret = c("specificity", "sensitivity"))
      )
      coords[, `:=`(
        FPR  = 1 - specificity,
        Tool = glue("{display_name} (AUC = {format(auc_val, digits = 3, nsmall = 3)})")
      )]
      coords
    }))
    
    ggplot(roc_dt, aes(FPR, sensitivity, colour = Tool)) +
      geom_line(size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
      coord_fixed() +
      scale_colour_brewer(palette = "Set1") +
      labs(title = panel_title,
           x = "1 – Specificity",
           y = "Sensitivity",
           colour = NULL) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position.inside = c(0.98, 0.02),
        legend.justification = c("right", "bottom"),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        legend.text = element_text(size = 10)
      ) +
      guides(colour = guide_legend(override.aes = list(size = 1.5)))
  }
  
  p_func <- build_roc_plot(data_functional, label_column_functional, "Functional dataset")
  p_clin <- build_roc_plot(data_clinical,  label_column_clinical,  "Clinical dataset")
  
  combined <- p_func + p_clin +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = "A",
                    theme = theme(plot.tag = element_text(face = "bold",
                                                          size = 16,
                                                          hjust = -0.1,
                                                          vjust = 1.2)))
  
  out_dir  <- "../results/figures/auc"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_file <- file.path(out_dir, "combined_functional_clinical_roc.png")
  ggsave(out_file, plot = combined, width = 14, height = 7, dpi = 300)
  
  message(glue("✔ Combined ROC figure saved to {out_file}"))
}

tools_to_plot <- c("FuncVEP_CTI", "FuncVEP_CTE", "FuncVEP_SP", "ClinVEP_CTI", "ClinVEP_CTE", "ClinVEP_SP")
save_combined_functional_clinical_roc(
  data_functional = benchmarking_data_functional,
  data_clinical = benchmarking_data_clinical,
  tool_names = tools_to_plot,
  label_column_functional = "functional_label",
  label_column_clinical = "clinical_label"
)



font_add(family = "Aptos",
         regular = "../resources/fonts/Aptos.ttf")
showtext_auto()
save_combined_distribution_plot <- function(data,
                                            tool_names,
                                            label_column,
                                            label_type = c("functional", "clinical"),
                                            show_threshold = TRUE,
                                            save_individual = TRUE,
                                            font_size_mult = 3,
                                            show_legends_individual = FALSE) {
  label_type <- match.arg(label_type)
  plots <- vector("list", length(tool_names))
  
  # --- shared theme with square legend keys and no extra margins ---
  theme_mpl_match <- function(size_mult = font_size_mult, base_family = "Aptos") {
    theme_minimal(base_family = base_family) +
      theme(
        plot.title      = element_text(size = 15 * size_mult, face = "plain", hjust = 0.5,
                                       margin = margin(b = 10 * size_mult)),
        axis.title.x    = element_text(size = 14 * size_mult, face = "plain",
                                       margin = margin(t = 8 * size_mult)),
        axis.title.y    = element_text(size = 14 * size_mult, face = "plain",
                                       margin = margin(r = 8 * size_mult)),
        axis.text.x     = element_text(size = 14 * size_mult, face = "plain"),
        axis.text.y     = element_text(size = 12 * size_mult, face = "plain"),
        legend.title    = element_text(size = 15 * size_mult, face = "plain"),
        legend.text     = element_text(size = 14 * size_mult, face = "plain"),
        legend.position = "bottom",
        legend.key.width  = unit(26, "pt"),
        legend.key.height = unit(12, "pt"),
        legend.margin = margin(0, 0, 0, 0)
      )
  }
  
  out_dir_root       <- "../results/figures/distribution_plots"
  out_dir_combined   <- file.path(out_dir_root, label_type)
  out_dir_individual <- out_dir_combined
  if (!dir.exists(out_dir_combined)) dir.create(out_dir_combined, recursive = TRUE)
  if (save_individual && !dir.exists(out_dir_individual)) dir.create(out_dir_individual, recursive = TRUE)
  
  legend_saved <- FALSE
  
  for (i in seq_along(tool_names)) {
    tool         <- tool_names[i]
    display_name <- gsub("_", "-", tool)
    
    plot_data <- data[!is.na(get(tool)) & !is.na(get(label_column)),
                      .(Score = get(tool),
                        Label = factor(get(label_column),
                                       levels = c(0, 1),
                                       labels = c("Benign", "Pathogenic")))]
    
    thr_val <- NA
    if (show_threshold) {
      roc_obj <- pROC::roc(as.numeric(plot_data$Label) - 1,
                           plot_data$Score, quiet = TRUE)
      thr_val <- as.numeric(pROC::coords(roc_obj, "best", ret = "threshold"))
    }
    
    p <- ggplot(plot_data, aes(Score, fill = Label)) +
      geom_histogram(position = "identity", bins = 50,
                     alpha = 0.6, colour = "black") +
      scale_fill_manual(values = c("Benign"     = "#3CB371",
                                   "Pathogenic" = "#E74C3C")) +
      labs(title = display_name,
           x     = "Prediction score",
           y     = "Variant count",
           fill  = NULL) +
      coord_cartesian(clip = "off") +
      theme_mpl_match()
    
    if (!is.na(thr_val)) {
      p <- p +
        geom_vline(xintercept = thr_val, linetype = "dashed",
                   colour = "blue", linewidth = 0.6) +
        annotate("text",
                 x      = thr_val,
                 y      = Inf,
                 vjust  = 1.5,
                 hjust  = 1.1,
                 label  = glue::glue("Thr = {format(thr_val, digits = 3, nsmall = 3)}"),
                 colour = "blue",
                 size   = (12 * font_size_mult) / ggplot2::.pt)
    }
    
    if (!legend_saved) {
      gtab <- ggplotGrob(p)
      legend_index <- which(sapply(gtab$grobs, function(x) x$name) == "guide-box")
      if (length(legend_index) > 0) {
        legend_grob <- gtab$grobs[[legend_index]]
        
        tmp_file <- tempfile(fileext = ".png")
        png(filename = tmp_file, width = 800, height = 200, units = "px", res = 150)
        grid::grid.newpage()
        grid::grid.draw(legend_grob)
        dev.off()
        
        library(magick)
        img <- image_read(tmp_file)
        img_trimmed <- image_trim(img)
        
        legend_out <- file.path(out_dir_root, "legend_only.png")
        image_write(img_trimmed, path = legend_out, format = "png")
        message(glue::glue("Legend saved (auto-trimmed) to {legend_out}"))
        
        
        legend_saved <- TRUE
      }
    }
    
    plots[[i]] <- p
    
    if (save_individual) {
      p_indiv <- if (!show_legends_individual) p + theme(legend.position = "none") else p
      out_file_individual <- file.path(out_dir_individual, paste0(display_name, "_distribution.png"))
      ggsave(out_file_individual, plot = p_indiv, width = 6, height = 5, dpi = 300, units = "in")
      message(glue::glue("Individual plot saved to {out_file_individual}"))
    }
  }
  
  combined <- patchwork::wrap_plots(plots, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    theme_mpl_match()
  
  out_file_combined <- file.path(out_dir_combined, "combined_score_distributions.png")
  ggsave(out_file_combined, plot = combined, width = 10, height = 8, dpi = 300, units = "in")
  message(glue::glue("Combined distribution plot saved to {out_file_combined}"))
}



save_combined_distribution_plot(
  benchmarking_data_functional,
  models,
  label_column = "functional_label",
  label_type = "functional",
  font_size_mult = 3,
  show_legends_individual = FALSE,
)  

save_combined_distribution_plot(
  benchmarking_data_clinical,
  models,
  label_column = "clinical_label",
  label_type = "clinical",
  font_size_mult = 3,
  show_legends_individual = FALSE,
)
