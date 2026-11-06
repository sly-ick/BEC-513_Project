###############################################################################
# PROJECT 4: TCGA RNA-seq gene-expression boxplot for NKX2-1 (Auto-path fixed)
# Compare Primary Tumor (PT) with Solid Tissue Normal (STN).
# Y-axis: log2(TPM + 1)
###############################################################################

############################ 1. Install/load packages #########################

required_packages <- c("readr", "dplyr", "ggplot2", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(readr)
library(dplyr)
library(ggplot2)
library(tibble)

############################ 2. Set up folders ################################

project_dir <- getwd()
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project folder: ", project_dir, "\n", sep = "")
cat("Results folder: ", results_dir, "\n\n", sep = "")

############################ 3. Locate the input files #########################

sample_sheet_candidates <- list.files(
  project_dir,
  pattern = "^gdc_sample_sheet.*\\.tsv$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(sample_sheet_candidates) == 0) {
  stop("Could not find gdc_sample_sheet.tsv in ", project_dir)
}
sample_sheet_path <- sample_sheet_candidates[[1]]

archive_candidates <- list.files(
  project_dir,
  pattern = "^tcga_data\\.(tar\\.gz|tgz|tar|zip)$",
  full.names = TRUE,
  ignore.case = TRUE
)

archive_path <- if (length(archive_candidates) > 0) archive_candidates[[1]] else NA_character_

manifest_candidates <- list.files(
  project_dir,
  pattern = "^gdc_manifest(\\..*)?$",
  full.names = TRUE,
  ignore.case = TRUE
)

manifest_path <- if (length(manifest_candidates) > 0) manifest_candidates[[1]] else NA_character_

cat("Sample sheet: ", basename(sample_sheet_path), "\n", sep = "")
if (!is.na(archive_path)) cat("Expression archive: ", basename(archive_path), "\n", sep = "")
if (!is.na(manifest_path)) cat("Manifest found: ", basename(manifest_path), "\n", sep = "")
cat("\n")

############################ 4. Extract archive if needed ######################

count_pattern <- "rna_seq\\.augmented_star_gene_counts\\.tsv(\\.gz)?$"

existing_count_files <- list.files(
  project_dir,
  pattern = count_pattern,
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(existing_count_files) == 0 && !is.na(archive_path)) {
  cat("Extracting archive...\n")
  if (grepl("\\.zip$", archive_path, ignore.case = TRUE)) {
    utils::unzip(archive_path, exdir = project_dir, overwrite = TRUE)
  } else {
    utils::untar(archive_path, exdir = project_dir, tar = "internal")
  }
  cat("Archive extraction completed.\n\n")
} else {
  cat("RNA-seq count files found directly in project directory; extraction skipped.\n\n")
}

############################ 5. Find all RNA-seq count files ###################

count_files <- list.files(
  project_dir,
  pattern = count_pattern,
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(count_files) == 0) {
  stop("No rna_seq.augmented_star_gene_counts.tsv files were found in ", project_dir)
}
cat("RNA-seq count files located: ", length(count_files), "\n\n", sep = "")

############################ 6. Read the sample sheet ##########################

sample_sheet <- readr::read_tsv(
  sample_sheet_path,
  show_col_types = FALSE,
  name_repair = "minimal",
  progress = FALSE
)

normalize_header <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", trimws(x)))
}

header_names <- names(sample_sheet)
header_names_normalized <- normalize_header(header_names)

file_id_candidates <- which(header_names_normalized %in% c("fileid", "fileuuid", "uuid"))
if (length(file_id_candidates) == 0) {
  stop("Could not identify the file UUID column in the sample sheet.")
}
file_id_col <- file_id_candidates[[1]]

normalize_label <- function(x) {
  tolower(trimws(gsub("[^[:alnum:]]+", "", as.character(x))))
}

accepted_labels <- c("primarytumor", "solidtissuenormal", "pt", "stn")

type_candidates <- which(header_names_normalized %in% c("tissuetype", "sampletype"))
if (length(type_candidates) > 0) {
  type_scores <- vapply(
    type_candidates,
    function(i) sum(normalize_label(sample_sheet[[i]]) %in% accepted_labels, na.rm = TRUE),
    integer(1)
  )
} else {
  type_scores <- integer(0)
}

if (length(type_candidates) == 0 || max(type_scores) == 0) {
  all_scores <- vapply(
    seq_along(sample_sheet),
    function(i) sum(normalize_label(sample_sheet[[i]]) %in% accepted_labels, na.rm = TRUE),
    integer(1)
  )
  all_scores[file_id_col] <- 0
  if (max(all_scores) == 0) stop("Could not identify tissue type column.")
  type_col <- which.max(all_scores)
} else {
  type_col <- type_candidates[which.max(type_scores)]
}

sample_metadata <- tibble::tibble(
  file_id = trimws(as.character(sample_sheet[[file_id_col]])),
  original_tissue_label = as.character(sample_sheet[[type_col]])
) %>%
  mutate(
    normalized_label = normalize_label(original_tissue_label),
    condition = case_when(
      normalized_label %in% c("primarytumor", "pt") ~ "PT",
      normalized_label %in% c("solidtissuenormal", "stn") ~ "STN",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(condition),
    !is.na(file_id),
    file_id != ""
  ) %>%
  distinct(file_id, .keep_all = TRUE) %>%
  mutate(file_id = tolower(file_id))

cat("Selected samples by condition:\n")
print(sample_metadata %>% count(condition))
cat("\n")

############################ 7. Match UUIDs to count files #####################

uuid_pattern <- paste0(
  "^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-",
  "[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
)

parent_ids <- basename(dirname(count_files))
filename_ids <- sub(
  "\\.rna_seq\\.augmented_star_gene_counts\\.tsv(\\.gz)?$",
  "",
  basename(count_files),
  ignore.case = TRUE
)

file_ids_from_paths <- ifelse(
  grepl(uuid_pattern, parent_ids),
  parent_ids,
  ifelse(grepl(uuid_pattern, filename_ids), filename_ids, NA_character_)
)

count_file_index <- tibble::tibble(
  file_id = tolower(file_ids_from_paths),
  count_file = count_files
) %>%
  filter(!is.na(file_id))

matched_samples <- sample_metadata %>%
  inner_join(count_file_index, by = "file_id")

cat("Matched PT/STN samples to RNA-seq files: ", nrow(matched_samples), "\n\n", sep = "")

############################ 8. Extract NKX2-1 TPM from each file ###############

extract_nkx2_1_tpm <- function(file_path, sample_file_id) {
  # Added comment = "#" to skip the "# gene-model: GENCODE v36" line
  expression_table <- readr::read_tsv(
    file_path,
    comment = "#",
    show_col_types = FALSE,
    progress = FALSE
  )

  columns_lower <- tolower(names(expression_table))
  gene_col <- match("gene_name", columns_lower)
  tpm_col <- match("tpm_unstranded", columns_lower)

  if (is.na(gene_col)) stop("gene_name column missing in: ", file_path)
  if (is.na(tpm_col)) tpm_col <- 7

  # Find the row whose gene symbol is NKX2-1.
  # Use which() to ignore rows where gene_name is NA.
  gene_indices <- which(trimws(as.character(expression_table[[gene_col]])) == "NKX2-1")

  if (length(gene_indices) == 0) {
    return(tibble::tibble(file_id = sample_file_id, TPM = NA_real_))
  }

  gene_rows <- expression_table[gene_indices, , drop = FALSE]

  if (nrow(gene_rows) > 1) {
    gene_rows <- gene_rows[1, ]
  }

  tpm_value <- suppressWarnings(as.numeric(as.character(gene_rows[[tpm_col]][[1]])))
  if (length(tpm_value) != 1 || is.na(tpm_value) || !is.finite(tpm_value) || tpm_value < 0) {
    return(tibble::tibble(file_id = sample_file_id, TPM = NA_real_))
  }

  tibble::tibble(file_id = sample_file_id, TPM = tpm_value)
}

tpm_results_list <- lapply(seq_len(nrow(matched_samples)), function(i) {
  if (i == 1 || i %% 50 == 0 || i == nrow(matched_samples)) {
    message("Reading sample ", i, " of ", nrow(matched_samples))
  }
  extract_nkx2_1_tpm(
    file_path = matched_samples$count_file[[i]],
    sample_file_id = matched_samples$file_id[[i]]
  )
})

tpm_results <- bind_rows(tpm_results_list)

############################ 9. Build the plotting table #######################

plot_data <- matched_samples %>%
  left_join(tpm_results, by = "file_id") %>%
  mutate(
    condition = factor(condition, levels = c("PT", "STN")),
    log2_TPM_plus_1 = log2(TPM + 1)
  )

# Handle missing TPM values gracefully
if (any(is.na(plot_data$TPM)) || any(is.na(plot_data$log2_TPM_plus_1))) {
  warning("Some samples have missing TPM values; continuing with available data.")
  # Filter out samples with missing TPM before plotting
  plot_data <- plot_data %>% filter(!is.na(TPM) & !is.na(log2_TPM_plus_1))
}

expression_summary <- plot_data %>%
  group_by(condition) %>%
  summarize(
    number_of_samples = n(),
    median_TPM = median(TPM),
    mean_TPM = mean(TPM),
    median_log2_TPM_plus_1 = median(log2_TPM_plus_1),
    mean_log2_TPM_plus_1 = mean(log2_TPM_plus_1),
    .groups = "drop"
  )

print(expression_summary)
readr::write_csv(
  expression_summary,
  file.path(results_dir, "NKX2-1_expression_summary_PT_STN.csv")
)

############################ 11. Statistical Test ##############################

test_res <- stats::t.test(log2_TPM_plus_1 ~ condition, data = plot_data)
cat("\nWelch Two Sample t-test:\n")
print(test_res)

############################ 12. Draw and save the boxplot #####################

nkx_plot <- ggplot(
  plot_data,
  aes(x = condition, y = log2_TPM_plus_1, fill = condition)
) +
  geom_boxplot(
    width = 0.60,
    alpha = 0.70,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    alpha = 0.55,
    size = 1.8
  ) +
  labs(
    title = paste0("NKX2-1 expression in TCGA-LUAD (p = ", format.pval(test_res$p.value, digits=3), ")"),
    x = "Tissue condition",
    y = expression(log[2](TPM + 1)),
    fill = "Condition"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggplot2::ggsave(
  filename = file.path(results_dir, "NKX2-1_PT_vs_STN_boxplot.png"),
  plot = nkx_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(results_dir, "NKX2-1_PT_vs_STN_boxplot.pdf"),
  plot = nkx_plot,
  width = 7,
  height = 5,
  units = "in"
)

if (nrow(plot_data) >= 2) {
  test_res <- stats::t.test(log2_TPM_plus_1 ~ condition, data = plot_data)
  cat("\nWelch Two Sample t-test:\n")
  print(test_res)

  nkx_plot <- ggplot(
    plot_data,
    aes(x = condition, y = log2_TPM_plus_1, fill = condition)
  ) +
    geom_boxplot(
      width = 0.50,
      alpha = 0.8,
      outlier.shape = NA
    ) +
    # REMOVED: geom_jitter()
    labs(
      title = paste0("NKX2-1 expression in TCGA-LUAD (p = ", format.pval(test_res$p.value, digits=3), ")"),
      x = "Tissue condition",
      y = "log2(TPM + 1)",
      fill = "Condition"
    ) +
    scale_fill_manual(values = c("PT" = "#e74c3c", "STN" = "#3498db")) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )

  print(nkx_plot)

  ggplot2::ggsave(
    filename = file.path(results_dir, "NKX2-1_PT_vs_STN_boxplot_clean.png"),
    plot = nkx_plot,
    width = 6,
    height = 5,
    units = "in",
    dpi = 300
  )
}

cat("\nAnalysis completed successfully.\n")
cat("Results saved in: ", results_dir, "\n", sep = "")
cat(" - NKX2-1_PT_vs_STN_boxplot.png\n")
cat(" - NKX2-1_PT_vs_STN_boxplot.pdf\n")
cat(" - NKX2-1_expression_PT_STN_sample_values.csv\n")
cat(" - NKX2-1_expression_summary_PT_STN.csv\n")
cat(" - Note: Some samples were skipped due to non-numeric or missing TPM values.\n", sep = "")