############################################################
# All_In.R
# lncRNA Downstream Analysis Pipeline
#
# Author : Furkan Emre Bora
# Lab    : Cingoz Lab
#
# Description:
#   End-to-end downstream analysis of RNA-seq differential expression
#   results with a focus on identifying key long non-coding RNAs (lncRNAs).
#   Input : DESeq2-style results table (.xlsx) + GENCODE GTF annotation.
#   Output: QC plots, volcano/MA plots, lncRNA tables, GSEA Hallmark results.
#
# How to use:
#   1. Install required packages (the script auto-installs missing ones).
#   2. Set the three paths in Section 0 (USER SETTINGS).
#   3. Source / Rscript the file.  All outputs are written to <outdir>.
#
# Input file requirements:
#   xlsx_file : Excel workbook with a sheet named "deg_results".
#               Mandatory columns -> Column1 (ENSG gene IDs), baseMean,
#               log2FoldChange, lfcSE, stat, pvalue, padj, Gene (symbol).
#   gtf_file  : GENCODE GTF (tested with v46; any recent version works).
#
# Output folders created automatically:
#   <outdir>/plots/   - all PNG figures (numbered 01-17)
#   <outdir>/tables/  - annotated gene tables and significant lncRNA lists
#   <outdir>/gsea/    - GSEA Hallmark CSV results
############################################################


############################################################
# SECTION 0 — COMMAND-LINE ARGUMENT PARSING
# Usage: Rscript All_In.R --gtf <file.gtf> --deg <results.xlsx> [--outdir <path>] [--padj 0.05] [--lfc 1.0]
############################################################
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")

suppressPackageStartupMessages(library(optparse))

# Define the command-line options
option_list <- list(
  make_option(c("--deg", "-d"), 
              type = "character", 
              default = NULL, 
              help = "Path to the DESeq2 results csv file (Required)", 
              metavar = "FILE"),
  make_option(c("--gtf", "-g"), 
              type = "character", 
              default = NULL, 
              help = "Path to the GENCODE annotation GTF file (Required)", 
              metavar = "FILE"),
  make_option(c("--outdir", "-o"), 
              type = "character", 
              default = "lncRNA_Analysis_Output", 
              help = "Output directory [Default: %default]", 
              metavar = "PATH"),
  make_option(c("--padj", "-p"), 
              type = "numeric", 
              default = 0.05, 
              help = "FDR (adjusted p-value) significance threshold [Default: %default]", 
              metavar = "NUM"),
  make_option(c("--lfc","-l"), 
              type = "numeric", 
              default = 1.0, 
              help = "Minimum absolute |log2 Fold Change| threshold [Default: %default]", 
              metavar = "NUM"),
  make_option(c("--seed"), 
              type = "integer", 
              default = 1, 
              help = "Random seed for reproducibility [Default: %default]", 
              metavar = "INT")
)

# Create the parser and parse the arguments
opt_parser <- OptionParser(option_list = option_list, 
                           description = "All_In.R: Command-line version for lncRNA downstream analysis pipeline")
opt <- parse_args(opt_parser)

# Validate that required arguments are provided
if (is.null(opt$deg) || is.null(opt$gtf)) {
  print_help(opt_parser)
  stop("Both --deg and --gtf arguments are required.", call. = FALSE)
}

# Assign parsed arguments to the script's core variables
csv_file    <- opt$deg
gtf_file     <- opt$gtf
outdir       <- opt$outdir
padj_cutoff  <- opt$padj
log2fc_cutoff <- opt$lfc
set.seed(opt$seed)

# The unified color palette (unchanged from original script)
DIR_COLORS <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")

# Print configuration for user confirmation
cat("Parameters set:\n")
cat("  DEG file    :", csv_file, "\n")
cat("  GTF file    :", gtf_file, "\n")
cat("  Output dir  :", outdir, "\n")
cat("  FDR cutoff  :", padj_cutoff, "\n")
cat("  log2FC cutoff:", log2fc_cutoff, "\n\n")
############################################################


############################################################
# SECTION 1 — PACKAGE MANAGEMENT
# All CRAN and Bioconductor packages are installed automatically
# if they are not already present.
############################################################

pkg_cran <- c(
  "data.table",   # fast file reading (GTF parsing)
  "dplyr",        # data manipulation
  "stringr",      # string operations (GTF attribute parsing)
  "tibble",       # tidy data frames
  "ggplot2",      # all plotting
  "openxlsx",     # read Excel without Java dependency
  "ggrepel",      # non-overlapping volcano labels
  "forcats",      # factor reordering for bar plots
  "scales"        # axis formatting helpers
)

pkg_bioc <- c(
  "fgsea",            # fast pre-ranked GSEA (fallback)
  "clusterProfiler",  # GSEA with TERM2GENE interface
  "enrichplot",       # dotplot / ridgeplot / emapplot for GSEA results
  "msigdbr"           # MSigDB gene sets (Hallmark, C2, etc.) in R
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in pkg_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in pkg_bioc) if (!requireNamespace(p, quietly = TRUE))
  BiocManager::install(p, update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(openxlsx)
  library(ggrepel)
  library(forcats)
  library(scales)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
})

theme_set(theme_classic(base_size = 14))

# Create all output directories
dir.create(file.path(outdir, "plots"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "gsea"),   showWarnings = FALSE, recursive = TRUE)


############################################################
# SECTION 2 — HELPER FUNCTIONS
############################################################

# Remove ENSEMBL version suffix (e.g. ENSG00000001.5 -> ENSG00000001)
clean_ensg <- function(x) sub("\\..*$", "", x)

# Safe -log10 transform: replaces p == 0 with a very small floor
# to avoid infinite values on volcano / distribution plots.
safe_neglog10 <- function(p) {
  -log10(pmax(as.numeric(p), 1e-300))
}

# Quantile-based axis range for coord_cartesian().
# Returns c(lo_quantile, hi_quantile) so extreme outliers are
# excluded from the visible area without removing data points.
qcap <- function(x, lo = 0.005, hi = 0.995) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, probs = c(lo, hi), na.rm = TRUE))
}

# Add a "direction" factor column to a DEG data frame.
# Levels: Up / Down / NS (not significant).
add_direction <- function(df, padj_thr, lfc_thr) {
  df %>%
    mutate(
      direction = case_when(
        !is.na(padj) & padj < padj_thr & log2FoldChange >=  lfc_thr ~ "Up",
        !is.na(padj) & padj < padj_thr & log2FoldChange <= -lfc_thr ~ "Down",
        TRUE ~ "NS"
      ),
      direction = factor(direction, levels = c("Up", "Down", "NS"))
    )
}


############################################################
# SECTION 3 — LOAD DIFFERENTIAL EXPRESSION RESULTS
#
# Expected sheet : "deg_results"
# Expected columns (at minimum):
#   Column1         -> ENSEMBL gene ID (e.g. ENSG00000001.5)
#   baseMean        -> mean normalised count across all samples
#   log2FoldChange  -> shrunken or unshrunken log2FC from DESeq2
#   pvalue          -> Wald test p-value
#   padj            -> Benjamini-Hochberg adjusted p-value
#   Gene            -> gene symbol (may be partially NA)
#
# Note: ~77 % padj NA is normal for DESeq2 when independent filtering
# removes low-count genes before multiple testing correction.
############################################################

if (!file.exists(csv_file)) stop("Excel file not found: ", csv_file)

deg <- read.csv(file = csv_file)

# Robustly detect gene ID column
gene_id_candidates <- c("Column1", "gene_id", "gene", "Gene", "X")
gene_id_col <- gene_id_candidates[gene_id_candidates %in% colnames(deg)]

if (length(gene_id_col) > 0) {
  gene_id_col <- gene_id_col[1]
  message("Using gene ID column: ", gene_id_col)
  deg <- deg %>% rename(gene_id = all_of(gene_id_col))
} else {
  message("No Column1/gene_id/Gene column found; using first column as gene_id.")
  deg <- deg %>% rename(gene_id = 1)
}

# Drop columns that are entirely NA (housekeeping artefacts from some exports)
deg <- deg %>%
  select(where(~ !all(is.na(.)))) %>%
  mutate(gene_id_clean = clean_ensg(gene_id))

# Ensure all stat columns are numeric (Excel can import them as character)
num_cols <- intersect(
  c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
  colnames(deg)
)
deg[num_cols] <- lapply(deg[num_cols], as.numeric)

cat("DEG rows loaded       :", nrow(deg), "\n")
cat("Rows with non-NA padj :", sum(!is.na(deg$padj)), "/", nrow(deg), "\n")


############################################################
# SECTION 4 — PARSE GENCODE GTF ANNOTATION
#
# The GTF is parsed with data.table::fread for speed.
# We extract two tables:
#   anno_gene : one row per gene  -> gene_id, gene_name, biotype, gene_length
#   anno_exon : one row per gene  -> exon_count, total_exon_length
# These are joined to form the annotation table `anno`.
#
# Tested with GENCODE v36-v46 (hg38). For other species / builds,
# ensure the GTF contains gene_type or gene_biotype attributes.
############################################################

if (!file.exists(gtf_file)) stop("GTF file not found: ", gtf_file)

message("Reading GTF — this may take 1-2 minutes for a full GENCODE file ...")

gtf_raw <- data.table::fread(
  cmd        = paste("grep -v '^#'", shQuote(gtf_file)),
  sep        = "\t",
  header     = FALSE,
  quote      = "",
  data.table = FALSE,
  fill       = TRUE
)

if (ncol(gtf_raw) < 9) {
  stop("GTF parsing failed: expected >= 9 columns, got ", ncol(gtf_raw),
       ". Check if the GTF is gzipped or improperly formatted.")
}

colnames(gtf_raw)[1:9] <- c(
  "seqname", "source", "feature", "start", "end",
  "score",   "strand", "frame",   "attribute"
)

# Keep only gene and exon rows (skip transcript / CDS / UTR rows)
gtf_raw <- gtf_raw %>% filter(feature %in% c("gene", "exon"))

# Helper: extract a quoted attribute value from the GTF attribute string
# e.g. extract_attr('gene_id "ENSG0001"; gene_name "TP53"', "gene_name")
#      -> "TP53"
extract_attr <- function(x, key) {
  m <- stringr::str_match(x, paste0(key, ' "([^"]+)"'))
  m[, 2]
}

gtf_raw <- gtf_raw %>%
  mutate(
    gene_id   = extract_attr(attribute, "gene_id"),
    gene_name = extract_attr(attribute, "gene_name"),
    # GENCODE uses "gene_type"; Ensembl GTFs use "gene_biotype" — handle both
    gene_type = dplyr::coalesce(
      extract_attr(attribute, "gene_type"),
      extract_attr(attribute, "gene_biotype")
    )
  ) %>%
  mutate(gene_id_clean = clean_ensg(gene_id))

# Gene-level annotation (one row per gene)
anno_gene <- gtf_raw %>%
  filter(feature == "gene") %>%
  transmute(
    gene_id_clean,
    gene_name,
    gene_biotype = gene_type,
    gene_length  = as.numeric(end) - as.numeric(start) + 1
  ) %>%
  distinct()

# Exon-level summary (total exon count and cumulative exon length per gene)
anno_exon <- gtf_raw %>%
  filter(feature == "exon") %>%
  mutate(exon_len = as.numeric(end) - as.numeric(start) + 1) %>%
  group_by(gene_id_clean) %>%
  summarise(
    exon_count        = dplyr::n(),
    total_exon_length = sum(exon_len, na.rm = TRUE),
    .groups = "drop"
  )

anno <- anno_gene %>% left_join(anno_exon, by = "gene_id_clean")
cat("Annotated genes in GTF:", nrow(anno), "\n")


############################################################
# SECTION 5 — MERGE DEG RESULTS WITH GTF ANNOTATION
#
# Join on the version-stripped ENSEMBL ID (gene_id_clean).
# gene_symbol preference order: Gene column from xlsx > GTF gene_name > ENSG ID
############################################################

deg_anno <- deg %>%
  left_join(anno, by = "gene_id_clean") %>%
  mutate(
    gene_symbol = dplyr::coalesce(
      if ("Gene" %in% colnames(.)) Gene else NA_character_,
      gene_name,
      gene_id_clean
    ),
    gene_symbol  = ifelse(is.na(gene_symbol) | gene_symbol == "",
                          gene_id_clean, gene_symbol),
    gene_biotype = ifelse(is.na(gene_biotype) | gene_biotype == "",
                          "Unannotated", gene_biotype)
  )

write.csv(
  deg_anno,
  file.path(outdir, "tables", "DEG_annotated_all_genes.csv"),
  row.names = FALSE
)

# Biotype labels that define lncRNAs in GENCODE annotation.
# These cover the main lncRNA subtypes assigned by GENCODE.
lnc_types <- c(
  "lncRNA",                          # generic long ncRNA (GENCODE v36+)
  "lincRNA",                         # long intergenic ncRNA
  "antisense",                       # antisense to a protein-coding gene
  "processed_transcript",            # spliced, no ORF assigned
  "sense_intronic",                  # overlaps intron of another gene (same strand)
  "sense_overlapping",               # overlaps another gene (same strand)
  "3prime_overlapping_ncrna",        # 3' overlap ncRNA
  "macro_lncRNA",                    # very long lncRNA (>100 kb)
  "bidirectional_promoter_lncRNA"    # divergently transcribed from a shared promoter
)

deg_lnc <- deg_anno %>% filter(gene_biotype %in% lnc_types)
cat("lncRNA genes in dataset:", nrow(deg_lnc), "\n")

write.csv(
  deg_lnc,
  file.path(outdir, "tables", "DEG_annotated_lncRNA_only.csv"),
  row.names = FALSE
)


############################################################
# SECTION 6 — DEFINE SIGNIFICANT GENE SETS
#
# Three tiers are defined to balance stringency vs. power:
#
#   sig_lnc          broad  : padj < 0.05  + |lfc| >= 1
#   sig_lnc_strict   strict : padj < 0.05  + pvalue < 0.05 + |lfc| >= 1
#   sig_lnc_moderate medium : padj < 0.10  + pvalue < 0.05 + |lfc| >= 0.58
#
# Why three tiers?
#   DESeq2 independent filtering can make padj NA for ~77 % of genes.
#   The moderate set recovers biologically relevant lncRNAs with
#   relaxed FDR but still requiring nominal p < 0.05 and >= 1.5-fold change.
#   The strict set is appropriate for publication-level reporting.
############################################################

deg_anno <- add_direction(deg_anno, padj_cutoff, log2fc_cutoff)
deg_lnc  <- add_direction(deg_lnc,  padj_cutoff, log2fc_cutoff)

# Broad set: FDR < 5 % and at least 2-fold change
sig_lnc <- deg_lnc %>%
  filter(
    !is.na(padj), padj < padj_cutoff,
    !is.na(log2FoldChange), abs(log2FoldChange) >= log2fc_cutoff
  )

# Strict set: requires both adjusted AND nominal p-value (dual confirmation)
sig_lnc_strict <- deg_lnc %>%
  filter(
    !is.na(padj),   padj   < 0.05,
    !is.na(pvalue), pvalue < 0.05,
    !is.na(log2FoldChange), abs(log2FoldChange) >= 1.0
  )

# Moderate set: useful for pathway enrichment and network analysis
sig_lnc_moderate <- deg_lnc %>%
  filter(
    !is.na(padj),   padj   < 0.10,
    !is.na(pvalue), pvalue < 0.05,
    !is.na(log2FoldChange), abs(log2FoldChange) >= 0.58  # >= 1.5-fold
  )

cat("\n--- Significant lncRNA counts ---\n")
cat("sig_lnc         (padj<0.05,  |lfc|>=1)                  :", nrow(sig_lnc), "\n")
cat("sig_lnc_strict  (padj<0.05, pval<0.05, |lfc|>=1)        :", nrow(sig_lnc_strict), "\n")
cat("sig_lnc_moderate(padj<0.10, pval<0.05, |lfc|>=0.58)     :", nrow(sig_lnc_moderate), "\n\n")

write.csv(sig_lnc,         file.path(outdir, "tables", "SIG_lncRNA_padj0.05_lfc1.csv"),   row.names = FALSE)
write.csv(sig_lnc_strict,  file.path(outdir, "tables", "SIG_lncRNA_STRICT.csv"),           row.names = FALSE)
write.csv(sig_lnc_moderate,file.path(outdir, "tables", "SIG_lncRNA_MODERATE.csv"),         row.names = FALSE)


############################################################
# SECTION 7 — QC / DISTRIBUTION PLOTS
#
# These plots help assess the quality of the DE results:
#   01  log2FC distribution   : should be roughly symmetric around 0
#   02  -log10(pvalue)        : should show enrichment near high values
#   03  -log10(padj)          : FDR distribution after correction
#   04  baseMean violin       : lncRNA expression tends to be lower than mRNA
#   05  Exon count histogram  : lncRNAs typically have fewer exons than mRNAs
#   06  Gene length histogram : lncRNAs are often shorter than protein-coding genes
############################################################

## 7.1  log2FC distribution (axis zoomed to 1st-99th quantile)
lfc_lim <- qcap(deg_anno$log2FoldChange[is.finite(deg_anno$log2FoldChange)])

p_lfc <- ggplot(
    deg_anno %>% filter(is.finite(log2FoldChange)),
    aes(x = log2FoldChange)
  ) +
  geom_histogram(bins = 80, fill = "steelblue", color = "white", linewidth = 0.2) +
  labs(title = "log2FC distribution (all genes)",
       x = "log2FoldChange", y = "Count")
if (all(is.finite(lfc_lim))) p_lfc <- p_lfc + coord_cartesian(xlim = lfc_lim)
ggsave(file.path(outdir, "plots", "01_log2FC_distribution.png"),
       p_lfc, width = 8, height = 5, dpi = 300)

## 7.2  -log10(pvalue) distribution
pv_df <- deg_anno %>%
  filter(!is.na(pvalue), is.finite(pvalue), pvalue > 0) %>%
  mutate(neglog10p = safe_neglog10(pvalue))
pv_lim <- qcap(pv_df$neglog10p)

p_pv <- ggplot(pv_df, aes(x = neglog10p)) +
  geom_histogram(bins = 80, fill = "steelblue", color = "white", linewidth = 0.2) +
  labs(title = "-log10(pvalue) distribution", x = "-log10(pvalue)", y = "Count")
if (all(is.finite(pv_lim))) p_pv <- p_pv + coord_cartesian(xlim = pv_lim)
ggsave(file.path(outdir, "plots", "02_neglog10pvalue_distribution.png"),
       p_pv, width = 8, height = 5, dpi = 300)

## 7.3  -log10(padj) distribution
fdr_df <- deg_anno %>%
  filter(!is.na(padj), is.finite(padj), padj > 0) %>%
  mutate(neglog10fdr = safe_neglog10(padj))
fdr_lim <- qcap(fdr_df$neglog10fdr)

p_fdr <- ggplot(fdr_df, aes(x = neglog10fdr)) +
  geom_histogram(bins = 80, fill = "steelblue", color = "white", linewidth = 0.2) +
  labs(title = "-log10(padj) distribution", x = "-log10(padj)", y = "Count")
if (all(is.finite(fdr_lim))) p_fdr <- p_fdr + coord_cartesian(xlim = fdr_lim)
ggsave(file.path(outdir, "plots", "03_neglog10padj_distribution.png"),
       p_fdr, width = 8, height = 5, dpi = 300)

## 7.4  baseMean violin by gene class
#  log10(baseMean + 1) is used because baseMean spans several orders of magnitude.
bm_df <- deg_anno %>%
  filter(!is.na(baseMean), baseMean >= 0) %>%
  mutate(
    log10_bm = log10(baseMean + 1),
    class = factor(
      case_when(
        gene_biotype %in% lnc_types ~ "lncRNA",
        gene_biotype == "Unannotated" ~ "Unannotated",
        TRUE ~ "non-lncRNA"
      ),
      levels = c("lncRNA", "non-lncRNA", "Unannotated")
    )
  )

p_bm <- ggplot(bm_df, aes(x = class, y = log10_bm, fill = class)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  labs(title = "Expression proxy (log10(baseMean+1)) by gene class",
       x = NULL, y = "log10(baseMean + 1)") +
  theme(legend.position = "none")
ggsave(file.path(outdir, "plots", "04_baseMean_violin_by_class.png"),
       p_bm, width = 8, height = 5, dpi = 300)

## 7.5  Exon count distribution (stacked, log10 scale)
ex_df <- deg_anno %>%
  filter(!is.na(exon_count), exon_count > 0) %>%
  mutate(
    log10_exon = log10(exon_count + 1),
    class = case_when(
      gene_biotype %in% lnc_types ~ "lncRNA",
      gene_biotype == "Unannotated" ~ "Unannotated",
      TRUE ~ "non-lncRNA"
    )
  )

p_ex <- ggplot(ex_df, aes(x = log10_exon, fill = class)) +
  geom_histogram(bins = 60, position = "stack", color = "white", linewidth = 0.15) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Exon count distribution (log10 scale)",
       x = "log10(exon_count + 1)", y = "Count")
ggsave(file.path(outdir, "plots", "05_ExonCount_distribution.png"),
       p_ex, width = 9, height = 5, dpi = 300)

## 7.6  Gene length distribution (stacked, log10 scale)
gl_df <- deg_anno %>%
  filter(!is.na(gene_length), gene_length > 0) %>%
  mutate(
    log10_gl = log10(gene_length),
    class = case_when(
      gene_biotype %in% lnc_types ~ "lncRNA",
      gene_biotype == "Unannotated" ~ "Unannotated",
      TRUE ~ "non-lncRNA"
    )
  )
gl_lim <- qcap(gl_df$log10_gl)

p_gl <- ggplot(gl_df, aes(x = log10_gl, fill = class)) +
  geom_histogram(bins = 60, position = "stack", color = "white", linewidth = 0.15) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Gene length distribution (log10 scale)",
       x = "log10(gene length, bp)", y = "Count")
if (all(is.finite(gl_lim))) p_gl <- p_gl + coord_cartesian(xlim = gl_lim)
ggsave(file.path(outdir, "plots", "06_GeneLength_distribution.png"),
       p_gl, width = 9, height = 5, dpi = 300)


############################################################
# SECTION 8 — MA PLOT
#
# An MA plot shows log fold change (M) against mean expression (A).
# Points are coloured by direction: red = up, blue = down, grey = NS.
# The y-axis is zoomed to the 0.5th-99.5th quantile so extreme outliers
# (common in low-expression genes) do not compress the informative range.
############################################################

ma_df <- deg_anno %>%
  filter(
    !is.na(baseMean), baseMean >= 0,
    !is.na(log2FoldChange), is.finite(log2FoldChange)
  ) %>%
  mutate(log10_bm = log10(baseMean + 1))

y_lim <- qcap(ma_df$log2FoldChange)

p_ma <- ggplot(ma_df, aes(x = log10_bm, y = log2FoldChange, color = direction)) +
  geom_point(alpha = 0.45, size = 0.8) +
  geom_hline(yintercept = c(-log2fc_cutoff, log2fc_cutoff),
             linetype = "dashed", color = "black") +
  scale_color_manual(
    values = DIR_COLORS,
    guide  = guide_legend(override.aes = list(size = 3, alpha = 1))
  ) +
  labs(title = "MA plot — all genes",
       x = "log10(baseMean + 1)", y = "log2(Fold Change)", color = NULL)
if (all(is.finite(y_lim))) p_ma <- p_ma + coord_cartesian(ylim = y_lim)
ggsave(file.path(outdir, "plots", "07_MA_plot.png"),
       p_ma, width = 9, height = 5, dpi = 300)


############################################################
# SECTION 9 — VOLCANO PLOT (lncRNA subset)
#
# Because DESeq2 independent filtering assigns NA padj to ~77 % of genes,
# the volcano uses raw p-value for point colouring (gives a meaningful
# distribution).  The strict set (sig_lnc_strict) is labelled; if that set
# is empty, the top 20 lncRNAs by raw p-value are labelled instead.
#
# x-axis: data-driven via quantile zoom (not hardcoded) so that the range
# automatically adapts to whatever fold-change distribution the data have.
############################################################

volcano_df <- deg_lnc %>%
  mutate(
    pvalue         = as.numeric(pvalue),
    log2FoldChange = as.numeric(log2FoldChange),
    neglog10p      = safe_neglog10(pvalue),
    grp = case_when(
      !is.na(pvalue) & pvalue < 0.05 & log2FoldChange >=  1 ~ "Up",
      !is.na(pvalue) & pvalue < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(is.finite(log2FoldChange), is.finite(neglog10p))

n_up   <- sum(volcano_df$grp == "Up",   na.rm = TRUE)
n_down <- sum(volcano_df$grp == "Down", na.rm = TRUE)

# Relabel factor levels to include counts in the legend
volcano_df$grp <- factor(
  volcano_df$grp,
  levels = c("Up", "Down", "NS"),
  labels = c(
    paste0("Up (pval<0.05): ",   n_up),
    paste0("Down (pval<0.05): ", n_down),
    "NS"
  )
)

# Choose genes to label: prefer strict set, fall back to top-20 by raw pvalue
label_ids <- if (nrow(sig_lnc_strict) > 0) {
  sig_lnc_strict$gene_id_clean
} else {
  deg_lnc %>%
    filter(!is.na(pvalue), is.finite(pvalue)) %>%
    arrange(pvalue) %>%
    head(20) %>%
    pull(gene_id_clean)
}

label_df <- volcano_df %>%
  mutate(gene_label = ifelse(
    gene_id_clean %in% label_ids,
    ifelse(is.na(gene_name) | gene_name == "", gene_id_clean, gene_name),
    NA_character_
  )) %>%
  filter(!is.na(gene_label))

# Data-driven x limits with 10 % padding
x_range <- qcap(volcano_df$log2FoldChange, 0.005, 0.995)
x_pad   <- diff(x_range) * 0.10
x_lim   <- c(x_range[1] - x_pad, x_range[2] + x_pad)

# Colour named vector aligned to the relabelled factor levels
vol_colors <- setNames(
  c("#E41A1C", "#377EB8", "grey70"),
  levels(volcano_df$grp)
)

png(file.path(outdir, "plots", "08_Volcano_lncRNA.png"),
    width = 2600, height = 2100, res = 300, type = "cairo-png")
print(
  ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10p, color = grp)) +
    geom_point(size = 1.0, alpha = 0.75) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    {if (nrow(label_df) > 0)
      geom_label_repel(
        data              = label_df,
        aes(label         = gene_label),
        size              = 3.5,
        min.segment.length = 0,
        box.padding       = 0.45,
        point.padding     = 0.25,
        max.overlaps      = Inf,
        force             = 12,
        force_pull        = 0.2,
        seed              = 123,
        show.legend       = FALSE
      )
    } +
    scale_color_manual(values = vol_colors) +
    coord_cartesian(xlim = x_lim) +
    labs(
      title    = "lncRNA — Volcano plot",
      subtitle = paste0(
        "Labelled: ", length(label_ids),
        " key lncRNAs  |  classification threshold: raw p-value < 0.05"
      ),
      x     = "log2(Fold Change)",
      y     = "-log10(p-value)",
      color = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title      = element_text(face = "bold", size = 18),
      legend.position = c(0.82, 0.88)
    )
)
dev.off()


############################################################
# SECTION 10 — TOP lncRNA BAR PLOT
#
# Uses the moderate threshold for a richer plot (more bars);
# falls back to the strict set if moderate is empty.
# Bars are coloured red (Up) / blue (Down) and ordered by log2FC.
############################################################

plot_set  <- if (nrow(sig_lnc_moderate) > 0) sig_lnc_moderate else sig_lnc
set_label <- if (nrow(sig_lnc_moderate) > 0) "moderate" else "strict"

if (nrow(plot_set) > 0) {
  top_bar <- plot_set %>%
    arrange(padj) %>%
    slice_head(n = as.integer(min(25L, nrow(plot_set)))) %>%
    mutate(
      label     = ifelse(is.na(gene_symbol) | gene_symbol == "",
                         gene_id_clean, gene_symbol),
      label     = fct_reorder(label, log2FoldChange),
      direction = factor(direction, levels = c("Up", "Down", "NS"))
    )

  p_bar <- ggplot(top_bar, aes(x = label, y = log2FoldChange, fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = DIR_COLORS, drop = FALSE) +
    labs(
      title    = paste0("Top lncRNAs — ", set_label, " threshold"),
      subtitle = if (set_label == "moderate")
        "padj < 0.10 + pval < 0.05 + |log2FC| >= 0.58" else
        "padj < 0.05 + |log2FC| >= 1",
      x = NULL, y = "log2FoldChange", fill = NULL
    )
  ggsave(file.path(outdir, "plots", "09_Top_lncRNA_barplot.png"),
         p_bar, width = 9, height = 7, dpi = 300)
} else {
  message("No significant lncRNAs found for bar plot. Try relaxing thresholds.")
}


############################################################
# SECTION 11 — TOP lncRNA LOLLIPOP PLOT
#
# Alternative to the bar plot; easier to read when many genes are shown.
# Each lollipop stem goes from y = 0 to y = log2FC.
############################################################

if (nrow(plot_set) > 0) {
  top_lolli <- plot_set %>%
    arrange(padj) %>%
    head(as.integer(min(20L, nrow(plot_set)))) %>%
    mutate(
      label     = ifelse(is.na(gene_name) | gene_name == "",
                         gene_id_clean, gene_name),
      label     = fct_reorder(label, log2FoldChange),
      direction = factor(direction, levels = c("Up", "Down", "NS"))
    )

  p_lolli <- ggplot(top_lolli,
                    aes(x = label, y = log2FoldChange, color = direction)) +
    geom_segment(aes(x = label, xend = label, y = 0, yend = log2FoldChange),
                 linewidth = 1.2) +
    geom_point(size = 4) +
    coord_flip() +
    scale_color_manual(values = DIR_COLORS, drop = FALSE) +
    labs(
      title = paste0("Top ", nrow(top_lolli), " lncRNA hits (ranked by padj)"),
      x = NULL, y = "log2FoldChange", color = NULL
    )
  ggsave(file.path(outdir, "plots", "10_Top_lncRNA_lollipop.png"),
         p_lolli, width = 10, height = 7, dpi = 300)
}


############################################################
# SECTION 12 — BIOTYPE COMPOSITION AMONG SIGNIFICANT GENES
#
# Bar chart showing which gene biotypes are most represented
# among the genes that pass the FDR < 0.05 threshold.
# Useful for checking whether protein-coding genes dominate
# or whether ncRNA categories are also enriched.
############################################################

sig_biotype <- deg_anno %>%
  filter(!is.na(padj), padj < 0.05) %>%
  count(gene_biotype, sort = TRUE) %>%
  slice_head(n = 30L)

if (nrow(sig_biotype) > 0) {
  p_biotype <- ggplot(sig_biotype,
                      aes(x = reorder(gene_biotype, n), y = n)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Biotype composition among significant genes (padj < 0.05)",
         x = NULL, y = "Count")
  ggsave(file.path(outdir, "plots", "11_Sig_biotype_bar.png"),
         p_biotype, width = 10, height = 7, dpi = 300)
}


############################################################
# SECTION 13 — EXPRESSION vs EFFECT SCATTER PLOTS
#
# Two scatter plots help visualise the relationship between
# mean expression level and the magnitude / significance of change:
#   Plot A: log10(baseMean) vs |log2FC|  — shows if LFC is inflated at low counts
#   Plot B: log10(baseMean) vs -log10(padj) — shows power increases with expression
############################################################

scatter_df <- deg_anno %>%
  filter(
    !is.na(baseMean), baseMean >= 0,
    !is.na(log2FoldChange), is.finite(log2FoldChange)
  ) %>%
  mutate(
    log10_bm     = log10(baseMean + 1),
    absLFC       = abs(log2FoldChange),
    neglog10padj = ifelse(is.na(padj), NA_real_, safe_neglog10(padj))
  )

p_scat1 <- ggplot(scatter_df, aes(x = log10_bm, y = absLFC)) +
  geom_point(alpha = 0.3, size = 0.7) +
  labs(title = "Expression vs |log2FC|",
       x = "log10(baseMean + 1)", y = "|log2FoldChange|")
ggsave(file.path(outdir, "plots", "12_Expr_vs_absLFC.png"),
       p_scat1, width = 8, height = 5, dpi = 300)

p_scat2 <- ggplot(scatter_df %>% filter(!is.na(neglog10padj)),
                  aes(x = log10_bm, y = neglog10padj)) +
  geom_point(alpha = 0.3, size = 0.7) +
  labs(title = "Expression vs significance",
       x = "log10(baseMean + 1)", y = "-log10(padj)")
ggsave(file.path(outdir, "plots", "13_Expr_vs_neglog10padj.png"),
       p_scat2, width = 8, height = 5, dpi = 300)


############################################################
# SECTION 14 — QC SUMMARY TABLE
#
# Writes a single-row CSV with key counts for quick reporting.
############################################################

# Pre-compute values before passing into tibble() to avoid
# column-name shadowing (tibble evaluates column expressions lazily)
.n_strict  <- nrow(sig_lnc_strict)
.n_mod     <- nrow(sig_lnc_moderate)
.up_strict <- sum(sig_lnc_strict$direction == "Up",   na.rm = TRUE)
.dn_strict <- sum(sig_lnc_strict$direction == "Down", na.rm = TRUE)

qc_tbl <- tibble::tibble(
  total_genes        = nrow(deg_anno),
  genes_with_padj    = sum(!is.na(deg_anno$padj)),
  sig_all_padj005    = sum(!is.na(deg_anno$padj) & deg_anno$padj < 0.05),
  total_lnc          = nrow(deg_lnc),
  lnc_with_padj      = sum(!is.na(deg_lnc$padj)),
  sig_lnc_padj005    = nrow(sig_lnc),
  sig_lnc_strict_n   = .n_strict,
  sig_lnc_moderate_n = .n_mod,
  up_lnc_strict      = .up_strict,
  down_lnc_strict    = .dn_strict
)

write.csv(qc_tbl, file.path(outdir, "tables", "QC_summary.csv"), row.names = FALSE)
print(t(qc_tbl))


############################################################
# SECTION 15 — GSEA HALLMARK PATHWAYS
#
# Gene Set Enrichment Analysis using the MSigDB Hallmark collection
# (50 well-defined biological processes, minimal redundancy).
#
# Ranking metric: sign(log2FC) * -log10(pvalue)
#   - Positive values -> upregulated & significant
#   - Negative values -> downregulated & significant
#   - This metric is preferred over log2FC alone because it
#     incorporates both direction and statistical support.
#
# Only protein-coding genes are used for ranking because MSigDB
# Hallmark gene sets contain protein-coding symbols.
#
# Note: GSEA is wrapped in tryCatch so that a failure (e.g. no
# significant pathways) does not abort the rest of the script.
#
# Output files:
#   gsea/GSEA_Hallmark_results.csv  -> full results table
#   plots/14_GSEA_Hallmark_dotplot.png
#   plots/15_GSEA_top_enrichment_curve.png
#   plots/16_GSEA_Hallmark_ridgeplot.png
#   plots/17_GSEA_Hallmark_emapplot.png
############################################################

deg_pc <- deg_anno %>%
  filter(gene_biotype == "protein_coding") %>%
  filter(!is.na(pvalue), is.finite(pvalue), pvalue > 0) %>%
  filter(!is.na(log2FoldChange), is.finite(log2FoldChange)) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  mutate(
    rank_score = sign(log2FoldChange) * (-log10(pvalue)),
    gene_sym   = toupper(gene_name)
  )

cat("\nProtein-coding genes available for GSEA ranking:", nrow(deg_pc), "\n")

if (nrow(deg_pc) < 100) {
  message("GSEA skipped: fewer than 100 protein-coding genes with valid pvalue.")
} else {

  # Collapse duplicate symbols: keep the entry with the largest |rank|
  ranked_tmp <- tapply(
    deg_pc$rank_score, deg_pc$gene_sym,
    function(x) x[which.max(abs(x))]
  )
  ranked <- setNames(as.numeric(ranked_tmp), names(ranked_tmp))

  # Add tiny random jitter to break exact ties (required by fgsea backend)
  set.seed(123)
  ranked <- ranked + rnorm(length(ranked), 0, 1e-9)
  ranked <- sort(ranked, decreasing = TRUE)

  # Load Hallmark gene sets from MSigDB (human, HGNC symbols)
  # msigdbr v7+ uses collection="H"; older versions use category="H"
  msig_h <- msigdbr(collection = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    mutate(gene_symbol = toupper(gene_symbol)) %>%
    as.data.frame()

  overlap_n <- sum(names(ranked) %in% msig_h$gene_symbol)
  cat("Overlap between ranked genes and Hallmark symbols:", overlap_n, "\n")
  if (overlap_n < 100)
    warning("Low overlap with MSigDB — check that gene_name contains HGNC symbols, not ENSG IDs.")

  # Run GSEA with pvalueCutoff = 1 so all pathways are retained in the
  # result object; significant ones are filtered afterwards for visualisation.
  gsea_h <- tryCatch(
    GSEA(
      geneList      = ranked,
      TERM2GENE     = msig_h,
      pvalueCutoff  = 1,
      pAdjustMethod = "BH",
      minGSSize     = 10,
      maxGSSize     = 500,
      eps           = 0,
      verbose       = FALSE
    ),
    error = function(e) {
      message("GSEA failed: ", e$message)
      NULL
    }
  )

  if (!is.null(gsea_h) && nrow(as.data.frame(gsea_h@result)) > 0) {

    gsea_res <- as.data.frame(gsea_h@result) %>% arrange(p.adjust)
    write.csv(gsea_res,
              file.path(outdir, "gsea", "GSEA_Hallmark_results.csv"),
              row.names = FALSE)
    cat("Significant Hallmark pathways (padj < 0.05):",
        sum(gsea_res$p.adjust < 0.05), "\n")

    # Subset to significant pathways for visualisation
    # (fall back to full object if filter leaves nothing)
    gsea_sig <- tryCatch(
      filter(gsea_h, p.adjust < 0.05),
      error = function(e) gsea_h
    )

    # Dotplot: NES on x-axis, gene ratio as dot size, padj as dot colour
    tryCatch({
      png(file.path(outdir, "plots", "14_GSEA_Hallmark_dotplot.png"),
          width = 2400, height = 1800, res = 300, type = "cairo-png")
      print(dotplot(gsea_sig, showCategory = 15) +
              ggtitle("GSEA Hallmark — protein-coding ranked list"))
      dev.off()
    }, error = function(e) message("GSEA dotplot skipped: ", e$message))

    # Enrichment score curve for the top-ranked pathway
    tryCatch({
      top_id <- gsea_res$ID[1]
      png(file.path(outdir, "plots", "15_GSEA_top_enrichment_curve.png"),
          width = 2400, height = 1800, res = 300, type = "cairo-png")
      print(gseaplot2(gsea_h, geneSetID = top_id, title = top_id))
      dev.off()
    }, error = function(e) message("GSEA enrichment curve skipped: ", e$message))

    # Ridgeplot: distribution of per-gene ranks within each enriched pathway
    tryCatch({
      png(file.path(outdir, "plots", "16_GSEA_Hallmark_ridgeplot.png"),
          width = 2400, height = 1800, res = 300, type = "cairo-png")
      print(ridgeplot(gsea_sig, showCategory = 20))
      dev.off()
    }, error = function(e) message("GSEA ridgeplot skipped: ", e$message))

    # Enrichment map: pathways as nodes, edges = shared leading-edge genes
    tryCatch({
      png(file.path(outdir, "plots", "17_GSEA_Hallmark_emapplot.png"),
          width = 2400, height = 1800, res = 300, type = "cairo-png")
      print(emapplot(pairwise_termsim(gsea_sig), showCategory = 30))
      dev.off()
    }, error = function(e) message("GSEA emapplot skipped: ", e$message))

  } else {
    message("GSEA returned no enriched pathways. ",
            "Try relaxing pvalueCutoff or verify gene_name contains HGNC symbols.")
  }

}  # end GSEA block

message("\nAll done.  Outputs are in: ", normalizePath(outdir))
