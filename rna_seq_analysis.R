#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(tximport)
  library(DESeq2)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

# -------------------------
# Argument Parser
# -------------------------
parser <- ArgumentParser(description = "RNA-seq downstream analysis script.")
parser$add_argument("--mode", type = "character", required = TRUE,
                    help = "Analysis mode: 'align' or 'quant'.")
parser$add_argument("--config", type = "character", required = TRUE,
                    help = "Path to the config.yaml file.")
parser$add_argument("--contrast", type = "character", required = TRUE,
                    help = "Contrast to test, format: 'treat,control'.")
parser$add_argument("--output_prefix", type = "character", required = TRUE,
                    help = "Prefix for result files.")
parser$add_argument("--input_counts", type = "character",
                    help = "Path to featureCounts file (align mode).")
parser$add_argument("--input_rsem_files", type = "character",
                    help = "Comma-separated list of RSEM .genes.results files (align mode).")
parser$add_argument("--input_quants", type = "character",
                    help = "Comma-separated list of Salmon quant.sf files (quant mode).")
parser$add_argument("--tx2gene", type = "character",
                    help = "Path to tx2gene mapping file (quant mode).")
parser$add_argument("--output_counts_matrix", type = "character", required = TRUE,
                    help = "Path to save the final gene symbol count matrix.")
parser$add_argument("--output_tpm_matrix", type = "character", required = TRUE,
                    help = "Path to save the final gene symbol TPM matrix.")

args <- parser$parse_args()

# -------------------------
# Helper: ID conversion + aggregate to gene symbols
# -------------------------
convert_ids_and_aggregate <- function(matrix_in, annot_db) {
  cat("Converting ENSEMBL IDs to Gene Symbols...\n")

  if (is.null(rownames(matrix_in)) || anyNA(rownames(matrix_in))) {
    stop("Input matrix has no rownames. Expected ENSEMBL IDs as rownames.")
  }

  ensembl_ids <- gsub("\\..*$", "", rownames(matrix_in))
  symbols <- mapIds(
    annot_db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  df_to_agg <- as.data.frame(matrix_in) %>%
    rownames_to_column("ensembl_version") %>%
    mutate(
      ensembl = gsub("\\..*$", "", ensembl_version),
      symbol = symbols[ensembl]
    ) %>%
    filter(!is.na(symbol) & symbol != "") %>%
    select(symbol, everything(), -ensembl_version, -ensembl)

  agg_df <- df_to_agg %>%
    group_by(symbol) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("symbol")

  cat("ID conversion complete. Original IDs:", nrow(matrix_in),
      "Mapped & Aggregated Symbols:", nrow(agg_df), "\n")
  agg_df
}

# -------------------------
# 1) Load Config and Prepare Metadata
# -------------------------
config <- yaml::read_yaml(args$config)

if (is.null(config$samples) || length(config$samples) == 0) {
  stop("config.yaml missing 'samples' section or it is empty.")
}

samples_df <- tibble(sample = names(config$samples)) %>%
  mutate(condition = sapply(sample, function(s) config$samples[[s]]$condition))

if (any(is.na(samples_df$condition)) || any(samples_df$condition == "")) {
  stop("Some samples have missing/empty condition in config.yaml.")
}

sample_names <- samples_df$sample

# -------------------------
# 2) Load Annotation DB
# -------------------------
cat("Loading Annotation Database...\n")
annot_db_name <- config$analysis$r_annotation_db
if (is.null(annot_db_name) || annot_db_name == "") {
  stop("config$analysis$r_annotation_db is missing/empty (e.g., org.Hs.eg.db).")
}
suppressPackageStartupMessages(library(annot_db_name, character.only = TRUE))
annot_db <- get(annot_db_name)

# -------------------------
# 3) Load and Process Data
# -------------------------
counts_data_raw <- NULL
tpm_data_raw <- NULL

if (args$mode == "align") {
  cat("Mode: align. Loading featureCounts and RSEM data.\n")

  if (is.null(args$input_counts) || args$input_counts == "") {
    stop("--input_counts is required in align mode.")
  }
  if (is.null(args$input_rsem_files) || args$input_rsem_files == "") {
    stop("--input_rsem_files is required in align mode.")
  }

  # featureCounts
  counts_data_raw <- read.table(args$input_counts, header = TRUE, row.names = 1,
                                sep = "\t", check.names = FALSE)
  # typical featureCounts: first 5 columns are annotation fields, then samples
  if (ncol(counts_data_raw) < 7) {
    stop("featureCounts file seems to have too few columns; expected >= 7.")
  }
  counts_data_raw <- counts_data_raw[, 6:ncol(counts_data_raw), drop = FALSE]

  clean_colnames <- basename(colnames(counts_data_raw))
  colnames(counts_data_raw) <- sub("\\.unique\\.dedup\\.bam$", "", clean_colnames)

  # keep only samples in config (and in that order)
  missing_in_counts <- setdiff(sample_names, colnames(counts_data_raw))
  if (length(missing_in_counts) > 0) {
    stop(paste0("These samples in config are missing in featureCounts columns: ",
                paste(missing_in_counts, collapse = ", ")))
  }
  counts_data_raw <- counts_data_raw[, sample_names, drop = FALSE]

  # RSEM TPM files
  cat("Loading and merging RSEM TPM files...\n")
  rsem_files <- strsplit(args$input_rsem_files, ",")[[1]]
  rsem_files <- trimws(rsem_files)

  if (length(rsem_files) != length(sample_names)) {
    stop(paste0("Number of --input_rsem_files (", length(rsem_files),
                ") must equal number of samples in config (", length(sample_names), ")."))
  }

  tpm_list <- map2(rsem_files, sample_names, function(file_path, sample_name) {
    read_tsv(file_path, show_col_types = FALSE) %>%
      select(gene_id, TPM) %>%
      rename(!!sample_name := TPM)
  })

  tpm_df_merged <- reduce(tpm_list, full_join, by = "gene_id")
  tpm_data_raw <- tpm_df_merged %>%
    column_to_rownames("gene_id") %>%
    as.matrix()

  cat("RSEM TPM data successfully merged.\n")

} else if (args$mode == "quant") {
  cat("Mode: quant. Loading Salmon data with tximport.\n")

  if (is.null(args$input_quants) || args$input_quants == "") {
    stop("--input_quants is required in quant mode.")
  }
  if (is.null(args$tx2gene) || args$tx2gene == "") {
    stop("--tx2gene is required in quant mode.")
  }

  quant_files <- strsplit(args$input_quants, ",")[[1]]
  quant_files <- trimws(quant_files)

  if (length(quant_files) != length(sample_names)) {
    stop(paste0("Number of --input_quants (", length(quant_files),
                ") must equal number of samples in config (", length(sample_names), ")."))
  }

  names(quant_files) <- sample_names

  tx2gene_df <- read_tsv(args$tx2gene, col_names = c("TXNAME", "GENEID"), show_col_types = FALSE)

  txi <- tximport(
    quant_files,
    type = "salmon",
    tx2gene = tx2gene_df,
    ignoreTxVersion = TRUE,
    ignoreAfterBar = TRUE
  )

  counts_data_raw <- round(txi$counts)
  tpm_data_raw <- txi$abundance

} else {
  stop("--mode must be either 'align' or 'quant'.")
}

# Ensure matrices have colnames and match config
if (is.null(colnames(counts_data_raw))) stop("counts_data_raw has no colnames.")
if (is.null(colnames(tpm_data_raw))) stop("tpm_data_raw has no colnames.")

# -------------------------
# 6) Save Aggregated Matrices (always)
# -------------------------
counts_symbol <- convert_ids_and_aggregate(counts_data_raw, annot_db)
tpm_symbol <- convert_ids_and_aggregate(tpm_data_raw, annot_db)

write.table(counts_symbol, file = args$output_counts_matrix,
            sep = "\t", quote = FALSE, col.names = NA)
write.table(tpm_symbol, file = args$output_tpm_matrix,
            sep = "\t", quote = FALSE, col.names = NA)
cat("Final symbol-keyed count and TPM matrices saved.\n")

# -------------------------
# 4) DESeq2 Analysis (run only if replicates exist)
# -------------------------
# Build coldata aligned to counts columns
coldata <- samples_df %>%
  filter(sample %in% colnames(counts_data_raw)) %>%
  distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

coldata <- coldata[colnames(counts_data_raw), , drop = FALSE]
coldata$condition <- factor(coldata$condition)

# Parse contrast
contrast_vec <- strsplit(args$contrast, ",")[[1]]
contrast_vec <- trimws(contrast_vec)
if (length(contrast_vec) != 2) {
  stop("--contrast must be in format 'treat,control'.")
}
treat <- contrast_vec[1]
control <- contrast_vec[2]

if (!(treat %in% levels(coldata$condition)) || !(control %in% levels(coldata$condition))) {
  stop(paste0(
    "Contrast levels not found in condition factor. Requested: ",
    treat, " vs ", control,
    ". Available: ", paste(levels(coldata$condition), collapse = ", ")
  ))
}

rep_table <- table(coldata$condition)
min_rep <- min(rep_table)

if (min_rep < 2) {
  msg <- paste0(
    "DESeq2 skipped: not enough biological replicates.\n",
    "Replicates per condition: ",
    paste(names(rep_table), rep_table, sep = "=", collapse = ", "),
    "\nNeed >= 2 per condition to estimate dispersion (DESeq2 no longer supports treating samples as replicates)."
  )
  cat(msg, "\n")

  # Write placeholder DESeq2 results
  placeholder <- data.frame(
    baseMean = rep(NA_real_, nrow(counts_data_raw)),
    log2FoldChange = rep(NA_real_, nrow(counts_data_raw)),
    lfcSE = rep(NA_real_, nrow(counts_data_raw)),
    stat = rep(NA_real_, nrow(counts_data_raw)),
    pvalue = rep(NA_real_, nrow(counts_data_raw)),
    padj = rep(NA_real_, nrow(counts_data_raw)),
    symbol = rep(NA_character_, nrow(counts_data_raw)),
    entrez = rep(NA_character_, nrow(counts_data_raw)),
    row.names = rownames(counts_data_raw)
  )
  write.csv(placeholder, file = paste0(args$output_prefix, ".deseq2_results.csv"))

  # Placeholders for plot/enrichment
  file.create(paste0(args$output_prefix, ".volcano_plot.png"))
  file.create(paste0(args$output_prefix, ".GO_enrichment.csv"))
  file.create(paste0(args$output_prefix, ".GO_dotplot.png"))

  cat("Analysis complete (counts/TPM done; DESeq2/Volcano/GO skipped).\n")
  quit(save = "no", status = 0)
}

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts_data_raw,
  colData = coldata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds, fitType = "local") # avoid some warnings

res <- results(dds, contrast = c("condition", treat, control))
res_df <- as.data.frame(res)

# -------------------------
# 5) Annotate Results
# -------------------------
res_df$symbol <- mapIds(
  annot_db,
  keys = gsub("\\..*", "", rownames(res_df)),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res_df$entrez <- mapIds(
  annot_db,
  keys = gsub("\\..*", "", rownames(res_df)),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

res_df_ordered <- res_df[order(res_df$padj), ]
write.csv(res_df_ordered, file = paste0(args$output_prefix, ".deseq2_results.csv"))

# -------------------------
# 7) Volcano Plot
# -------------------------
cat("Generating Volcano Plot...\n")

keyvals <- ifelse(!is.na(res_df$padj) & res_df$log2FoldChange < -1 & res_df$padj < 0.05, "blue",
                  ifelse(!is.na(res_df$padj) & res_df$log2FoldChange > 1 & res_df$padj < 0.05, "red", "grey"))
keyvals[is.na(keyvals)] <- "grey"
names(keyvals)[keyvals == "red"] <- "Up-regulated"
names(keyvals)[keyvals == "grey"] <- "Non-significant"
names(keyvals)[keyvals == "blue"] <- "Down-regulated"

selected_genes <- res_df_ordered %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  head(5) %>%
  pull(symbol)

png(paste0(args$output_prefix, ".volcano_plot.png"), width = 10, height = 10, units = "in", res = 300)
print(
  EnhancedVolcano(
    res_df,
    lab = res_df$symbol,
    x = "log2FoldChange",
    y = "padj",
    colCustom = keyvals,
    colAlpha = 0.8,
    selectLab = selected_genes,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    max.overlaps = Inf,
    pCutoff = 0.05,
    FCcutoff = 1.0,
    title = paste(treat, "vs", control),
    subtitle = "padj < 0.05 & |log2FC| > 1 (color); y-axis uses pvalue",
    legendPosition = "right",
    legendLabSize = 12,
    legendIconSize = 4.0,
    pointSize = 3.0,
    labSize = 4.0
  )
)
dev.off()

# -------------------------
# 8) GO Enrichment Analysis
# -------------------------
cat("Running GO Enrichment Analysis...\n")

fdr_thr <- config$analysis$fdr_threshold
if (is.null(fdr_thr) || !is.numeric(fdr_thr)) {
  fdr_thr <- 0.05
  cat("config$analysis$fdr_threshold missing/invalid; using default 0.05\n")
}

sig_genes <- res_df_ordered %>%
  filter(!is.na(padj), padj < fdr_thr, !is.na(entrez)) %>%
  pull(entrez)

all_genes <- res_df_ordered %>%
  filter(!is.na(entrez)) %>%
  pull(entrez)

if (length(sig_genes) > 10) {
  ego <- enrichGO(
    gene = sig_genes,
    universe = all_genes,
    OrgDb = annot_db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file = paste0(args$output_prefix, ".GO_enrichment.csv"))
    p <- dotplot(ego, showCategory = 20) + ggtitle("GO Enrichment Analysis")
    ggsave(paste0(args$output_prefix, ".GO_dotplot.png"), plot = p, width = 10, height = 8, dpi = 300)
  } else {
    file.create(paste0(args$output_prefix, ".GO_enrichment.csv"))
    file.create(paste0(args$output_prefix, ".GO_dotplot.png"))
  }
} else {
  cat("Not enough significant genes for GO analysis. Skipping.\n")
  file.create(paste0(args$output_prefix, ".GO_enrichment.csv"))
  file.create(paste0(args$output_prefix, ".GO_dotplot.png"))
}

cat("Analysis complete.\n")