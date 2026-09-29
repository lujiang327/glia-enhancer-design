#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) stop("Missing argument: ", flag)
  args[[position + 1L]]
}

library_path <- normalizePath(arg_value("--library"), mustWork = FALSE)
manifest_dir <- normalizePath(arg_value("--manifest-dir"), mustWork = FALSE)
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_path, .libPaths()))

cran <- "https://cloud.r-project.org"
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = cran, lib = library_path)
}
BiocManager::install(version = "3.19", ask = FALSE, update = FALSE, lib = library_path)
BiocManager::install(
  "edgeR",
  version = "3.19",
  ask = FALSE,
  update = FALSE,
  dependencies = c("Depends", "Imports", "LinkingTo"),
  lib = library_path
)

required <- c("BiocManager", "BiocVersion", "edgeR", "limma", "locfit", "Rcpp")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing installed packages: ", paste(missing, collapse = ", "))

versions <- data.frame(
  package = required,
  version = vapply(required, function(x) as.character(packageVersion(x)), character(1)),
  stringsAsFactors = FALSE
)
write.table(
  versions,
  file = file.path(manifest_dir, "installed_packages.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(manifest_dir, "sessionInfo.txt"))
jsonlite::write_json(
  list(
    status = "PHASE3_EDGER_RUNTIME_READY",
    R = R.version.string,
    Bioconductor = as.character(BiocManager::version()),
    library = "library",
    packages = setNames(as.list(versions$version), versions$package)
  ),
  file.path(manifest_dir, "runtime_manifest.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
cat("Installed isolated Phase 3 edgeR runtime.\n")
