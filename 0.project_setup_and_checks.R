# TBI paper-order analysis stage 0
# See the script index for figure coverage and execution dependencies.

# Run from Rscript, RStudio Source, or this notebook; no setwd is required.
# UTF-8 is required for Greek subtype labels and Chinese file paths.
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  .locale_candidates <- if (.Platform$OS.type == "windows") c("English_United States.utf8", ".UTF-8") else c("C.UTF-8", "en_US.UTF-8")
  for (.locale_candidate in .locale_candidates) {
    suppressWarnings(try(Sys.setlocale("LC_CTYPE", .locale_candidate), silent = TRUE))
    if (isTRUE(l10n_info()[["UTF-8"]])) break
  }
  if (!isTRUE(l10n_info()[["UTF-8"]])) stop("A UTF-8 LC_CTYPE locale is required. Configure an available UTF-8 locale before running this script.")
}
if (!exists("TBI_ROOT", inherits = FALSE)) {
  .root_env <- Sys.getenv("TBI_PROJECT_ROOT", unset = "")
  .file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
  .this_file <- if (length(.source_files)) tail(.source_files, 1)[[1]] else if (length(.file_arg)) sub("^--file=", "", .file_arg[1]) else ""
  TBI_ROOT <- if (nzchar(.root_env)) normalizePath(.root_env, mustWork = TRUE) else
    if (nzchar(.this_file)) dirname(normalizePath(.this_file, mustWork = TRUE)) else normalizePath(getwd(), mustWork = TRUE)
}
if (!file.exists(file.path(TBI_ROOT, "config", "settings.R"))) stop("Set TBI_PROJECT_ROOT to the folder containing the numbered scripts.")
source(file.path(TBI_ROOT, "config", "settings.R"), encoding = "UTF-8")
for (.module in c("common", "preprocessing", "plots", "tdeg", "cellchat", "trajectory", "validation"))
  source(file.path(TBI_ROOT, "functions", paste0(.module, ".R")), encoding = "UTF-8")


for (.d in c("inputs/objects", "inputs/reference", "inputs/matrices", "results/audit")) ensure_dir(project_path(.d))
.packages <- c("Seurat", "DoubletFinder", "harmony", "clusterProfiler", "org.Mm.eg.db", "mgcv", "CellChat", "pheatmap", "dplyr", "tidyr", "ggplot2", "patchwork", "scales", "openxlsx", "ComplexUpset", "slingshot", "monocle", "SingleCellExperiment", "S4Vectors", "Biobase", "VGAM", "Matrix", "future")
.check <- data.frame(package = .packages, installed = vapply(.packages, requireNamespace, logical(1), quietly = TRUE))
write_csv(.check, file.path(outdir("audit"), "dependency_preflight.csv"))
print(.check)
# Metadata from an existing integrated object is authoritative; do not guess sham Day from Lib.
.main_path <- input_path(cfg$input_main)
if (file.exists(.main_path) && requireNamespace("Seurat", quietly = TRUE)) {
  .obj <- normalize_metadata(readRDS(.main_path))
  .m <- unique(.obj@meta.data[, c("Lib", "Day", "Group", "Source", "Mice_per_library", "Control")])
  .template <- read_csv(input_path(cfg$sample_manifest))
  .extra <- setdiff(.template$Lib, .m$Lib)
  if (length(.extra)) warning("Template libraries absent from object: ", paste(.extra, collapse = ", "))
  .m$matrix_dir <- .template$matrix_dir[match(.m$Lib, .template$Lib)]
  .m$gene_column <- .template$gene_column[match(.m$Lib, .template$Lib)]
  .m$matrix_dir[is.na(.m$matrix_dir)] <- paste0("inputs/matrices/", .m$Lib[is.na(.m$matrix_dir)])
  .m$gene_column[is.na(.m$gene_column)] <- ifelse(.m$Source[is.na(.m$gene_column)] == "Public", 2L, 1L)
  write_csv(.m, file.path(outdir("audit"), "sample_manifest_from_actual_object.csv"))
  message("Review results/audit/sample_manifest_from_actual_object.csv and use it as config/sample_manifest.csv for raw-matrix QC.")
}
message("No DNBC4tools/SoupX run was reconstructed. For existing objects start at script 3 after configuring inputs; for matrices run 1 and 2 first.")

