run_release_checks <- function() {
  need("dplyr"); obj <- get_main(); d <- outdir("audit")
  checks <- list()
  record <- function(test, pass, detail) { checks[[length(checks) + 1L]] <<- data.frame(test = test, pass = pass, detail = detail) }
  record("deleted cells absent", !any(obj@misc$tbi_deleted_ids %in% colnames(obj)), "Cell IDs, not renamed labels")
  record("no level4", !"level4" %in% names(obj@meta.data), "Only level1/2/3 are used")
  record("no banned level3", !any(obj$level3 %in% c("Macrophage13", paste0("Neutrophil", 10:13), "CD8-ODG", "Microglia9")), "Old Mphi12 is tracked by cell IDs")
  libs <- unique(obj@meta.data[, c("Lib", "Day", "Group")])
  record("31 libraries", nrow(libs) == cfg$expected_libraries, paste(nrow(libs), "observed"))
  observed <- data.frame(level2 = names(table(obj$level2)), cells = as.integer(table(obj$level2)))
  write_csv(observed, file.path(d, "final_level2_counts.csv"))
  record("16 level2 types", length(unique(obj$level2)) == 16, paste(length(unique(obj$level2)), "observed"))
  record("58 confirmed global level3 labels", setequal(unique(as.character(obj$level3)), global_level3_labels()), paste(length(unique(obj$level3)), "observed; erythroblast states remain local"))
  confirmed <- read_csv(project_path("config", "confirmed_level3_counts.csv"))
  observed3 <- as.data.frame(table(obj$level3), stringsAsFactors = FALSE); names(observed3) <- c("level3", "observed")
  compare3 <- merge(confirmed, observed3, by = "level3", all = TRUE)
  compare3$difference <- compare3$observed - compare3$confirmed_cells
  write_csv(compare3, file.path(d, "confirmed_level3_count_comparison.csv"))
  for (lineage in names(subtype_orders)) {
    sub <- get_lineage(lineage)
    record(paste(lineage, "state coverage"), setequal(unique(sub$level3), subtype_orders[[lineage]]), paste(sort(unique(sub$level3)), collapse = ";"))
  }
  cc <- get_cellchat(); record("7 full CellChat objects", length(cc) == 7L, paste(names(cc), collapse = ";"))
  record("no CellChat downsampling", identical(attr(cc, "downsampled"), FALSE), "Only min50 subtype exclusion is permitted")
  tdeg_files <- list.files(project_path("results", "tDEG"), pattern = "tDEG_final.rds$", recursive = TRUE, full.names = TRUE)
  record("9 tDEG results", length(tdeg_files) == 9L, paste(length(tdeg_files), "observed"))
  keys <- names(cfg$tdeg_module_order)
  observed_keys <- vapply(tdeg_files, function(file) readRDS(file)$key, character(1))
  record("exact nine tDEG identities", !anyDuplicated(observed_keys) && setequal(observed_keys, keys), paste(observed_keys, collapse = ";"))
  for (file in tdeg_files) {
    value <- readRDS(file)
    parts <- strsplit(value$key, "_", fixed = TRUE)[[1]]
    current <- tryCatch({ read_tdeg_result(parts[1], parts[2]); TRUE }, error = function(e) FALSE)
    record(paste(value$key, "current subtype and settings"), current, "Main, subtype-object and settings fingerprints")
    record(paste(value$key, "candidate containment"), all(value$final_genes %in% value$candidates), "All final tDEGs must be pairwise DEGs")
    record(paste(value$key, "current annotation"), identical(value$main_md5, hash_input(main_final_path())), "Input fingerprint")
  }
  tables <- list.files(project_path("results", "Supplementary_tables"), pattern = "^TableS[0-9]+.*\\.xlsx$", full.names = TRUE)
  record("11 supplement workbooks", length(tables) == 11L, paste(length(tables), "observed"))
  table_ids <- sub("^(TableS[0-9]+).*", "\\1", basename(tables))
  record("exact TableS1-S11 identities", !anyDuplicated(table_ids) && setequal(table_ids, paste0("TableS", 1:11)), paste(table_ids, collapse = ";"))
  for (stage in as.character(3:36)) {
    path <- project_path("results", "audit", paste0(stage, ".status.rds"))
    passed <- file.exists(path) && identical(readRDS(path)$status, "success")
    record(paste("stage", stage, "completed"), passed, "Latest recorded execution must succeed")
  }
  # Manuscript counts are comparison targets, never forced totals after subtype exclusions.
  manuscript <- data.frame(tissue = c("All", "Brain", "Skull"), manuscript = c(256134L, 226647L, 29487L),
                             observed = c(ncol(obj), sum(obj$Group == "Brain"), sum(obj$Group == "Skull")))
  manuscript$difference <- manuscript$observed - manuscript$manuscript
  write_csv(manuscript, file.path(d, "manuscript_cell_count_comparison.csv"))
  result <- do.call(rbind, checks); write_csv(result, file.path(d, "release_checks.csv"))
  writeLines(capture.output(utils::sessionInfo()), file.path(d, "final_sessionInfo.txt"))
  versions <- vapply(c("Seurat", "DoubletFinder", "harmony", "clusterProfiler", "org.Mm.eg.db", "mgcv", "CellChat", "pheatmap", "dplyr", "tidyr", "data.table", "stringr", "ggplot2", "patchwork", "ComplexUpset", "slingshot", "monocle", "SingleCellExperiment", "openxlsx"), function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "NOT_INSTALLED"
  }, character(1))
  write_csv(data.frame(package = names(versions), version = versions), file.path(d, "actual_package_versions.csv"))
  if (any(!result$pass)) stop("Release checks failed; see results/audit/release_checks.csv.")
  invisible(result)
}
