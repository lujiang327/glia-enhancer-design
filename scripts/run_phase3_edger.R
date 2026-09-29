#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) stop("Missing argument: ", flag)
  args[[position + 1L]]
}

matrix_path <- arg_value("--count-matrix")
cell_type_path <- arg_value("--cell-types")
results_dir <- arg_value("--results-dir")
report_dir <- arg_value("--report-dir")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "contrasts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "leave_one_donor_out"), recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

donors <- c("LGS1", "LGS2", "LGS3", "LVG1")
cell_type_table <- read.delim(cell_type_path, stringsAsFactors = FALSE, check.names = FALSE)
groups <- cell_type_table$analysis_group
if (length(groups) != 13L || anyDuplicated(groups) || groups[[1]] != "MG") {
  stop("Expected 13 unique analysis groups with MG first")
}
off_targets <- groups[groups != "MG"]

input <- read.delim(gzfile(matrix_path), stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("peak_id", "chrom", "start", "end")
if (!identical(names(input)[seq_along(required_columns)], required_columns)) {
  stop("Unexpected count-matrix coordinate columns")
}
coordinates <- input[, required_columns]
counts <- as.matrix(input[, -(seq_along(required_columns))])
storage.mode(counts) <- "integer"
rownames(counts) <- coordinates$peak_id
if (anyNA(counts) || any(counts < 0L)) stop("Counts must be nonnegative integers")
if (anyDuplicated(coordinates$peak_id)) stop("Duplicate peak IDs in count matrix")
sample_ids <- colnames(counts)
parts <- strsplit(sample_ids, "__", fixed = TRUE)
if (any(lengths(parts) != 2L)) stop("Malformed sample ID")
metadata <- data.frame(
  sample_id = sample_ids,
  cell_type = vapply(parts, `[[`, character(1), 1L),
  donor = vapply(parts, `[[`, character(1), 2L),
  stringsAsFactors = FALSE
)
rownames(metadata) <- sample_ids
expected_ids <- as.vector(t(outer(groups, donors, paste, sep = "__")))
if (!identical(sample_ids, expected_ids)) stop("Count-matrix sample order does not match contract")
metadata$cell_type <- factor(metadata$cell_type, levels = groups)
metadata$donor <- factor(metadata$donor, levels = donors)

design <- model.matrix(~ 0 + cell_type + donor, data = metadata)
colnames(design) <- sub("^cell_type", "", colnames(design))
rownames(design) <- sample_ids
if (qr(design)$rank != ncol(design)) stop("Differential design is not full rank")

y0 <- DGEList(counts = counts, samples = metadata, genes = coordinates)
keep <- filterByExpr(y0, design = design)
if (sum(keep) < 2L) stop("Fewer than two peaks passed filterByExpr")
y <- y0[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

contrast_names <- c("MG_vs_rest", paste0("MG_vs_", off_targets))
contrast_matrix <- matrix(
  0,
  nrow = ncol(design),
  ncol = length(contrast_names),
  dimnames = list(colnames(design), contrast_names)
)
contrast_matrix["MG", "MG_vs_rest"] <- 1
contrast_matrix[off_targets, "MG_vs_rest"] <- -1 / length(off_targets)
for (group in off_targets) {
  name <- paste0("MG_vs_", group)
  contrast_matrix["MG", name] <- 1
  contrast_matrix[group, name] <- -1
}
if (any(abs(colSums(contrast_matrix[groups, , drop = FALSE])) > 1e-12)) {
  stop("Cell-type contrast weights do not sum to zero")
}
if (any(contrast_matrix[setdiff(rownames(contrast_matrix), groups), ] != 0)) {
  stop("Donor coefficients must be zero in all contrasts")
}

write.table(
  cbind(sample_id = rownames(design), as.data.frame(design, check.names = FALSE)),
  file.path(report_dir, "design_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  cbind(coefficient = rownames(contrast_matrix), as.data.frame(contrast_matrix, check.names = FALSE)),
  file.path(report_dir, "contrast_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

raw_library <- colSums(counts)
sample_qc <- data.frame(
  metadata,
  raw_input_library_size = as.numeric(raw_library),
  filtered_library_size = y$samples$lib.size,
  norm_factor = y$samples$norm.factors,
  effective_library_size = y$samples$lib.size * y$samples$norm.factors,
  stringsAsFactors = FALSE
)
write.table(
  sample_qc,
  file.path(report_dir, "normalization_factors.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

log_cpm <- cpm(y, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
donor_effects <- function(comparator) {
  effects <- sapply(donors, function(donor) {
    target_column <- which(metadata$cell_type == "MG" & metadata$donor == donor)
    if (comparator == "rest") {
      comparison_columns <- which(metadata$cell_type != "MG" & metadata$donor == donor)
      log_cpm[, target_column] - rowMeans(log_cpm[, comparison_columns, drop = FALSE])
    } else {
      comparison_column <- which(metadata$cell_type == comparator & metadata$donor == donor)
      log_cpm[, target_column] - log_cpm[, comparison_column]
    }
  })
  colnames(effects) <- paste0("donor_logFC_", donors)
  effects
}

write_gzip_table <- function(table, path) {
  connection <- gzfile(path, open = "wt", compression = 1)
  on.exit(close(connection))
  write.table(table, connection, sep = "\t", quote = FALSE, row.names = FALSE)
}

contrast_summaries <- list()
pvalue_matrix <- matrix(NA_real_, nrow = nrow(y), ncol = length(contrast_names), dimnames = list(rownames(y), contrast_names))
full_glial_stats <- list()
primary_test <- NULL
for (contrast_name in contrast_names) {
  qlf <- glmQLFTest(fit, contrast = contrast_matrix[, contrast_name])
  table <- topTags(qlf, n = Inf, sort.by = "none")$table
  if (!identical(rownames(table), rownames(y))) stop("edgeR result order changed")
  comparator <- if (contrast_name == "MG_vs_rest") "rest" else sub("^MG_vs_", "", contrast_name)
  effects <- donor_effects(comparator)
  output <- cbind(
    coordinates[match(rownames(table), coordinates$peak_id), , drop = FALSE],
    table,
    effects,
    donors_positive = rowSums(effects > 0),
    minimum_donor_logFC = apply(effects, 1, min)
  )
  output_path <- file.path(results_dir, "contrasts", paste0(contrast_name, ".tsv.gz"))
  write_gzip_table(output, output_path)
  pvalue_matrix[, contrast_name] <- table$PValue
  significant <- table$FDR < 0.05
  positive <- significant & table$logFC > 0
  negative <- significant & table$logFC < 0
  contrast_summaries[[contrast_name]] <- list(
    tested_peaks = nrow(table),
    fdr_lt_0_05 = sum(significant),
    positive_fdr_lt_0_05 = sum(positive),
    negative_fdr_lt_0_05 = sum(negative),
    positive_fdr_lt_0_05_logFC_ge_1 = sum(positive & table$logFC >= 1),
    positive_all_four_donors = sum(positive & rowSums(effects > 0) == 4L),
    result = output_path
  )
  if (contrast_name %in% c("MG_vs_Astrocyte", "MG_vs_Microglia")) {
    full_glial_stats[[contrast_name]] <- data.frame(
      peak_id = rownames(table), logFC = table$logFC, FDR = table$FDR,
      stringsAsFactors = FALSE
    )
  }
  if (contrast_name == "MG_vs_rest") primary_test <- qlf
}

summary_rows <- do.call(rbind, lapply(names(contrast_summaries), function(name) {
  values <- contrast_summaries[[name]]
  data.frame(
    contrast = name,
    tested_peaks = values$tested_peaks,
    fdr_lt_0_05 = values$fdr_lt_0_05,
    positive_fdr_lt_0_05 = values$positive_fdr_lt_0_05,
    negative_fdr_lt_0_05 = values$negative_fdr_lt_0_05,
    positive_fdr_lt_0_05_logFC_ge_1 = values$positive_fdr_lt_0_05_logFC_ge_1,
    positive_all_four_donors = values$positive_all_four_donors,
    stringsAsFactors = FALSE
  )
}))
write.table(summary_rows, file.path(report_dir, "contrast_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

run_lodo <- function(excluded_donor, comparator) {
  selected <- metadata$donor != excluded_donor
  lodo_metadata <- droplevels(metadata[selected, , drop = FALSE])
  lodo_design <- model.matrix(~ 0 + cell_type + donor, data = lodo_metadata)
  colnames(lodo_design) <- sub("^cell_type", "", colnames(lodo_design))
  if (qr(lodo_design)$rank != ncol(lodo_design)) stop("LODO design is not full rank")
  lodo_y <- DGEList(counts = y$counts[, selected, drop = FALSE], samples = lodo_metadata)
  lodo_y <- calcNormFactors(lodo_y, method = "TMM")
  lodo_y <- estimateDisp(lodo_y, lodo_design, robust = TRUE)
  lodo_fit <- glmQLFit(lodo_y, lodo_design, robust = TRUE)
  contrast <- rep(0, ncol(lodo_design)); names(contrast) <- colnames(lodo_design)
  contrast["MG"] <- 1; contrast[comparator] <- -1
  lodo_test <- glmQLFTest(lodo_fit, contrast = contrast)
  table <- topTags(lodo_test, n = Inf, sort.by = "none")$table
  name <- paste0("MG_vs_", comparator)
  output <- data.frame(
    peak_id = rownames(table), table, stringsAsFactors = FALSE, check.names = FALSE
  )
  output_path <- file.path(results_dir, "leave_one_donor_out", paste0(name, "_without_", excluded_donor, ".tsv.gz"))
  write_gzip_table(output, output_path)
  full <- full_glial_stats[[name]]
  stopifnot(identical(full$peak_id, rownames(table)))
  full_positive <- full$FDR < 0.05 & full$logFC > 0
  if (!any(full_positive)) stop("No full-model positive peaks for ", name)
  data.frame(
    contrast = name,
    excluded_donor = excluded_donor,
    pearson_logFC = cor(full$logFC, table$logFC, method = "pearson"),
    spearman_logFC = cor(full$logFC, table$logFC, method = "spearman"),
    sign_concordance_all_peaks = mean(sign(full$logFC) == sign(table$logFC)),
    full_positive_peaks = sum(full_positive),
    fraction_full_positive_remaining_positive = mean(table$logFC[full_positive] > 0),
    fraction_full_positive_remaining_fdr = mean(table$FDR[full_positive] < 0.05 & table$logFC[full_positive] > 0),
    result = output_path,
    stringsAsFactors = FALSE
  )
}

lodo_rows <- do.call(rbind, lapply(c("Astrocyte", "Microglia"), function(comparator) {
  do.call(rbind, lapply(donors, run_lodo, comparator = comparator))
}))
write.table(lodo_rows, file.path(report_dir, "leave_one_donor_out_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

group_colors <- setNames(grDevices::hcl.colors(length(groups), "Dynamic"), groups)
donor_pch <- setNames(c(16, 15, 17, 18), donors)
png(file.path(report_dir, "mds.png"), width = 1800, height = 1400, res = 180)
plotMDS(y, top = min(20000, nrow(y)), labels = as.character(metadata$donor),
        col = group_colors[as.character(metadata$cell_type)], pch = donor_pch[as.character(metadata$donor)])
legend("topright", legend = groups, col = group_colors, pch = 16, cex = 0.65, ncol = 2)
dev.off()

png(file.path(report_dir, "normalization.png"), width = 2200, height = 1300, res = 180)
par(mfrow = c(2, 1), mar = c(9, 4, 2, 1))
barplot(sample_qc$raw_input_library_size / 1e6, names.arg = sample_ids, las = 2, cex.names = 0.45,
        col = group_colors[as.character(metadata$cell_type)], ylab = "Raw peak counts (millions)")
barplot(sample_qc$norm_factor, names.arg = sample_ids, las = 2, cex.names = 0.45,
        col = group_colors[as.character(metadata$cell_type)], ylab = "TMM normalization factor")
dev.off()

png(file.path(report_dir, "bcv.png"), width = 1600, height = 1300, res = 180)
plotBCV(y)
dev.off()
png(file.path(report_dir, "ql_dispersion.png"), width = 1600, height = 1300, res = 180)
plotQLDisp(fit)
dev.off()
png(file.path(report_dir, "primary_md.png"), width = 1600, height = 1300, res = 180)
plotMD(primary_test, status = decideTestsDGE(primary_test, adjust.method = "BH", p.value = 0.05),
       main = "MG vs equal-weight off-target mean")
abline(h = c(-1, 1), col = "grey50", lty = 2)
dev.off()
png(file.path(report_dir, "pvalue_histograms.png"), width = 2200, height = 1800, res = 180)
par(mfrow = c(4, 4), mar = c(3, 3, 2, 1))
for (name in contrast_names) hist(pvalue_matrix[, name], breaks = 40, main = name, xlab = "P-value")
dev.off()

writeLines(capture.output(sessionInfo()), file.path(report_dir, "sessionInfo.txt"))
result <- list(
  status = "DIFFERENTIAL_MODEL_DIAGNOSTICS_READY_FOR_REVIEW",
  candidate_selection_status = "BLOCKED_PENDING_DIFFERENTIAL_MODEL_QC_REVIEW",
  samples = ncol(y),
  peaks_before_filtering = nrow(y0),
  peaks_after_filterByExpr = nrow(y),
  design_rows = nrow(design),
  design_columns = ncol(design),
  design_rank = qr(design)$rank,
  design_condition_number = kappa(design),
  normalization = list(
    method = "TMM",
    factor_range = range(y$samples$norm.factors),
    product = prod(y$samples$norm.factors)
  ),
  dispersion = list(
    common_dispersion = y$common.dispersion,
    common_biological_coefficient_of_variation = sqrt(y$common.dispersion),
    tagwise_dispersion_quantiles = as.list(setNames(as.numeric(quantile(y$tagwise.dispersion)), names(quantile(y$tagwise.dispersion)))),
    posterior_ql_dispersion_quantiles = as.list(setNames(as.numeric(quantile(fit$var.post)), names(quantile(fit$var.post))))
  ),
  contrast_summaries = contrast_summaries,
  leave_one_donor_out = split(lodo_rows, seq_len(nrow(lodo_rows))),
  package_versions = list(
    R = R.version.string,
    edgeR = as.character(packageVersion("edgeR")),
    limma = as.character(packageVersion("limma"))
  ),
  scope = "Model fitting and diagnostics only. No candidate ranking or selection was performed."
)
write_json(result, file.path(report_dir, "differential_model_qc.json"), auto_unbox = TRUE, pretty = TRUE, digits = 10)
cat(toJSON(result, auto_unbox = TRUE, pretty = TRUE, digits = 10), "\n")
