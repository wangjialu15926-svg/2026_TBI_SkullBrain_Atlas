# TBI paper-order analysis stage 7
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

with_stage("07", {
  run_lineage_panels("Macrophage", "Figure2_S2")
})
