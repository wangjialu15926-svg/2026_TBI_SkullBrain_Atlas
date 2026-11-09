options(stringsAsFactors = FALSE)
project_path <- function(...) file.path(TBI_ROOT, ...)
input_path <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) stop("An input path must be one nonempty string.")
  x <- path.expand(x)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", x) || startsWith(x, "\\\\")) x else project_path(x)
}
need <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". Install them in the analysis environment; scripts do not install packages.")
}
ensure_dir <- function(path) { dir.create(path, recursive = TRUE, showWarnings = FALSE); path }
outdir <- function(...) ensure_dir(project_path("results", ...))
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
read_csv <- function(path, ...) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", ...)
read_input <- function(path) {
  path <- input_path(path)
  if (!file.exists(path)) stop("Required input is missing: ", path)
  readRDS(path)
}
check_columns <- function(x, required, context) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(context, " lacks columns: ", paste(missing, collapse = ", "))
}
hash_input <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) stop("Cannot fingerprint missing input: ", path)
  value <- unname(tools::md5sum(path))
  if (is.na(value)) stop("Cannot read input for fingerprinting: ", path)
  value
}
audit_run <- function(stage, inputs = character(), settings = list()) {
  d <- outdir("audit")
  saveRDS(list(stage = stage, time = Sys.time(), settings = settings,
               inputs = data.frame(path = inputs, md5 = vapply(inputs, hash_input, character(1))),
               session = utils::sessionInfo()), file.path(d, paste0(stage, ".provenance.rds")))
  writeLines(capture.output(utils::sessionInfo()), file.path(d, paste0(stage, ".sessionInfo.txt")))
}
day_order <- c("D0", "D1", "D7", "M6")
canonical_labels <- function(x) {
  x <- as.character(x)
  map <- c("B1-like MBC" = "B1 Cell", "B1-like_MBC" = "B1 Cell",
           "erythroblast1" = "Erythroblast1", "erythroblast2" = "Erythroblast2",
           "erythrocyte" = "Erythrocyte", "CD4 T" = "CD4 T Cell", "CD8 T" = "CD8 T Cell",
           "NK" = "NK Cell", "γδ T" = "γδ T Cell", "VLMC" = "VLMCs")
  hit <- !is.na(x) & x %in% names(map); x[hit] <- map[x[hit]]
  x <- sub("^Mφ([0-9]+)$", "Macrophage\\1", x)
  x <- sub("^Neu([0-9]+)$", "Neutrophil\\1", x)
  x <- sub("^MG([0-9]+)$", "Microglia\\1", x)
  x[x == "Mono"] <- "Monocyte"
  x
}
level2_order <- c("Promyelocyte", "Macrophage", "Microglia", "Neutrophil", "Basophil", "DC",
                  "B Cell", "T Cell", "Erythrocyte", "Astrocyte", "Endothelial Cell", "Pericyte",
                  "Oligodendrocyte", "Neuron", "CPEC", "Ependymal Cell")
level2_colors <- setNames(c("#a30059", "#960018", "#e60010", "#ff2f80", "#ffae97", "#cc5500",
                            "#1CE6FF", "#377eb8", "#FFFF00", "#556b2f", "#228b22", "#66ff00",
                            "#a5d854", "#fdbe02", "#66c1a5", "#007f66"), level2_order)
subtype_orders <- list(
  Macrophage = c("Monocyte", paste0("Macrophage", 1:12)),
  Neutrophil = paste0("Neutrophil", 1:9), Microglia = paste0("Microglia", 1:8),
  Bcell = c("B progenitors-1", "B progenitors-2", "Pro-B Cell", "Pre-B Cell", "B Cell", "B1 Cell", "Plasma Cell"),
  Tcell = c("CD4 T Cell", "CD8 T Cell", "CD8-IFN", "NK Cell", "T progenitors", "γδ T Cell"),
  Erythroid = c("Erythroblast1", "Erythroblast2", "Erythrocyte"),
  NVU = c("homeostatic_Astro", "DAA", "VLMCs", "SMC", "homeostatic_pericyte", "homeostatic_endo", "Angiogenic-like_endo")
)
subtype_colors <- c(
  setNames(c("#FFD6D6", "#FFBABA", "#FFA0A0", "#FF8585", "#FF6B6B", "#FF5252", "#E63946",
             "#D62828", "#C1121F", "#A4161A", "#960018", "#800014", "#4D000C"), subtype_orders$Macrophage),
  setNames(c("#FFE0EE", "#FFCCE0", "#FFB3D9", "#FF99CC", "#FF80BF", "#FF66B3", "#FF4DA6", "#FF2F80", "#E62973"), subtype_orders$Neutrophil),
  setNames(c("#B32D00", "#CC3D00", "#E65016", "#FF6B3D", "#FF9364", "#FFB08B", "#FFCDB2", "#FFE5D9"), subtype_orders$Microglia),
  setNames(c("#97fdfd", "#3ce2fe", "#00c5cd", "#018cc5", "#2a52be", "#1e3a8a", "#0000ff"), subtype_orders$Bcell),
  setNames(c("#00B060", "#0A7A35", "#28A745", "#90E87A", "#C4F291", "#E8F9B8"), subtype_orders$Tcell),
  setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628"), subtype_orders$NVU),
  setNames(c("#9c93e5", "#ff83a5", "#ec3f64"), subtype_orders$Erythroid))
module_colors <- c(C1 = "#E74C3C", C2 = "#3498DB", C3 = "#2ECC71", C4 = "#9B59B6",
                   C5 = "#E67E22", C6 = "#F1C40F", C7 = "#A0522D", C8 = "#FF69B4")
marker_level2 <- list(Promyelocyte = c("Elane", "Mpo"), Macrophage = c("Fn1", "Lyz2", "Arg1"),
  Microglia = c("Cx3cr1", "P2ry12", "Tmem119"), Neutrophil = c("S100a9", "S100a8", "Ngp"),
  Basophil = c("Ms4a2", "Fcer1a", "Itga2b"), DC = c("Bst2", "Cd74"), `B Cell` = c("Cd79a", "Cd19"),
  `T Cell` = c("Cd3d", "Trbc2"), Erythrocyte = c("Hba-a1", "Hbb-bs", "Aqp1"),
  Astrocyte = c("Slc1a2", "Slc1a3"), `Endothelial Cell` = c("Cdh5", "Kdr", "Flt1"),
  Pericyte = c("Rgs5", "Pdgfrb", "Vtn"), Oligodendrocyte = c("Mbp", "Plp1", "Mal"),
  Neuron = c("Stmn2", "Ccnd2", "Sox4"), CPEC = c("Ttr", "Otx2"), `Ependymal Cell` = c("Foxj1", "Tmem212"))
lineage_parent <- list(Macrophage = "Macrophage", Neutrophil = "Neutrophil", Microglia = "Microglia",
                       Bcell = "B Cell", Tcell = "T Cell", Erythroid = "Erythrocyte",
                       NVU = c("Astrocyte", "Endothelial Cell", "Pericyte"))
global_level3_labels <- function() unique(c(
  setdiff(unlist(subtype_orders, use.names = FALSE), c("Erythroblast1", "Erythroblast2")),
  setdiff(level2_order, unlist(lineage_parent, use.names = FALSE))))
global_erythroid_labels <- function(labels) {
  labels <- canonical_labels(labels)
  labels[labels %in% subtype_orders$Erythroid] <- "Erythrocyte"
  labels
}
assert_no_incomplete_stages <- function(stages) {
  for (stage in as.character(stages)) {
    path <- project_path("results", "audit", paste0(stage, ".status.rds"))
    # Older imported results may have no status file; a recorded failed/running
    # execution must never be bypassed by accepting the preceding cached file.
    if (file.exists(path) && !identical(readRDS(path)$status, "success"))
      stop("Script ", stage, " did not complete successfully; rerun it before consuming its results.")
  }
  invisible(TRUE)
}
main_final_path <- function() project_path("results", "objects", "main_final.rds")
get_main <- function() {
  need("Seurat")
  assert_no_incomplete_stages("3")
  obj <- read_input(main_final_path())
  if (!identical(obj@misc$tbi_schema, "TBI_final_level3_v1")) stop("Run script 3 to construct the audited final object.")
  if (any(!as.character(obj$level3) %in% global_level3_labels())) stop("Final global level3 labels do not match the confirmed 58-state vocabulary; rerun script 3. Erythroblast states belong only to the local erythroid object.")
  # Imported objects can retain integrated/SCT as their active assay.
  # Expression panels and marker tests in this project use normalized RNA.
  Seurat::DefaultAssay(obj) <- "RNA"
  obj
}
get_lineage <- function(name, tissue = NULL) {
  need("Seurat")
  obj <- read_input(project_path("results", "objects", paste0(name, "_final.rds")))
  main <- get_main()
  if (!identical(obj@misc$tbi_main_md5, hash_input(main_final_path()))) stop("Subtype input fingerprint is stale or missing; rerun script 3: ", name)
  if (!all(colnames(obj) %in% colnames(main))) stop("Stale subtype object: ", name)
  for (field in c("Lib", "Day", "Group", "level2", "level3")) {
    expected <- as.character(main@meta.data[colnames(obj), field])
    observed <- as.character(obj@meta.data[[field]])
    if (name == "Erythroid" && field == "level3") {
      if (any(!observed %in% subtype_orders$Erythroid)) stop("Invalid local erythroid subtype labels.")
      observed <- global_erythroid_labels(observed)
    }
    if (!identical(observed, expected)) stop("Subtype ", field, " metadata disagree with main: ", name)
  }
  if (!is.null(tissue)) obj <- subset(obj, cells = colnames(obj)[obj$Group == tissue])
  if (!is.null(tissue)) obj@misc$tbi_scope_tissues <- tissue
  if (!ncol(obj)) stop("No cells for ", name, " ", tissue)
  Seurat::DefaultAssay(obj) <- "RNA"
  obj
}
normalize_metadata <- function(obj, require_level2 = TRUE) {
  meta <- obj@meta.data
  ids <- colnames(obj)
  if (is.null(ids) || anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids))) stop("Main object has invalid or duplicated cell IDs.")
  if (is.null(rownames(meta)) || anyDuplicated(rownames(meta)) || !setequal(rownames(meta), ids)) stop("Cell IDs and metadata row names do not match.")
  meta <- meta[ids, , drop = FALSE]
  if (!"Lib" %in% names(meta) && "sample_batch_library" %in% names(meta)) meta$Lib <- as.character(meta$sample_batch_library)
  check_columns(meta, c("Lib", "Day", "Group"), "Object metadata")
  meta$Lib <- as.character(meta$Lib); meta$Day <- as.character(meta$Day); meta$Group <- as.character(meta$Group)
  if (anyNA(meta[, c("Lib", "Day", "Group")]) || any(!nzchar(meta$Lib))) stop("Missing Lib/Day/Group metadata.")
  if (any(!meta$Day %in% day_order) || any(!meta$Group %in% c("Brain", "Skull"))) stop("Unexpected Day/Group; do not infer from library names.")
  perlib <- unique(meta[, c("Lib", "Day", "Group")]); if (anyDuplicated(perlib$Lib)) stop("A library maps to multiple Day/Group values.")
  meta$Source <- ifelse(grepl("^TBI", meta$Lib), "In-house", "Public")
  meta$Dataset <- ifelse(meta$Source == "Public", "Public Brain", paste("In-house", meta$Group))
  meta$Mice_per_library <- ifelse(grepl("^TBI", meta$Lib), 4L, NA_integer_)
  if (!"Control" %in% names(meta) || anyNA(meta$Control) || any(!as.character(meta$Control) %in% c("sham", "naive", "injured"))) {
    meta$Control <- ifelse(meta$Day == "D0", ifelse(meta$Source == "In-house", "sham", "naive"), "injured")
  }
  if (any(meta$Group == "Skull" & meta$Day == "M6") || any(meta$Source == "In-house" & meta$Day == "M6")) stop("M6 must be public Brain data in this study.")
  if (require_level2) {
    check_columns(meta, "level2", "Object metadata")
    meta$level2 <- as.character(meta$level2)
    if (anyNA(meta$level2) || any(!meta$level2 %in% level2_order)) stop("Unrecognized level2 labels.")
    meta$level1 <- ifelse(meta$level2 %in% c("Promyelocyte", "Macrophage", "Microglia", "Neutrophil", "Basophil", "DC"), "Myeloid",
                           ifelse(meta$level2 %in% c("B Cell", "T Cell"), "Lymphoid", "Non-immune"))
  }
  obj@meta.data <- meta
  obj
}
save_plot <- function(p, path, width = 8, height = 6) {
  need("ggplot2")
  device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(path, plot = p, width = width, height = height, device = device, limitsize = FALSE)
}
with_pdf_device <- function(path, width, height, code) {
  device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  device(path, width = width, height = height)
  opened <- grDevices::dev.cur()
  on.exit({ if (opened %in% grDevices::dev.list()) grDevices::dev.off(opened) }, add = TRUE)
  force(code)
}
build_library_proportions <- function(meta, label, allowed = NULL, library_manifest = NULL) {
  required <- c("Lib", "Day", "Group", "Source", "Mice_per_library")
  check_columns(meta, c(required, label), "Proportion metadata")
  if (!nrow(meta)) stop("Cannot calculate proportions for an empty cell set.")
  meta$cell_type <- as.character(meta[[label]])
  if (anyNA(meta[, c("Lib", "Day", "Group", "cell_type")])) stop("Missing library, day, tissue or cell-type metadata.")
  if (is.null(allowed)) allowed <- sort(unique(meta$cell_type))
  allowed <- as.character(allowed)
  if (!length(allowed) || anyNA(allowed) || anyDuplicated(allowed) || any(!meta$cell_type %in% allowed)) stop("Proportion state list must cover every observed state exactly once.")
  libraries <- if (is.null(library_manifest)) unique(meta[, required, drop = FALSE]) else {
    check_columns(library_manifest, required, "Library reference manifest")
    unique(library_manifest[, required, drop = FALSE])
  }
  libraries$Lib <- as.character(libraries$Lib)
  if (anyNA(libraries[, c("Lib", "Day", "Group")]) || anyDuplicated(libraries$Lib)) stop("Each reference library must map to exactly one day and tissue.")
  if (any(!meta$Lib %in% libraries$Lib)) stop("Observed cells have a library absent from the reference manifest.")
  matched <- libraries[match(meta$Lib, libraries$Lib), , drop = FALSE]
  if (any(as.character(meta$Day) != as.character(matched$Day)) || any(as.character(meta$Group) != as.character(matched$Group))) stop("Cell Day/Group values disagree with the library reference.")
  grid <- merge(libraries, data.frame(cell_type = allowed), by = NULL)
  counts <- as.data.frame(table(factor(meta$Lib, levels = libraries$Lib), factor(meta$cell_type, levels = allowed)), stringsAsFactors = FALSE)
  names(counts) <- c("Lib", "cell_type", "n")
  grid <- merge(grid, counts, by = c("Lib", "cell_type"), all.x = TRUE, sort = FALSE)
  denominators <- setNames(as.integer(table(factor(meta$Lib, levels = libraries$Lib))), libraries$Lib)
  grid$denominator <- unname(denominators[grid$Lib])
  grid$percentage <- ifelse(grid$denominator > 0, 100 * grid$n / grid$denominator, NA_real_)
  grid$status <- ifelse(grid$denominator > 0, "included", "no_parent_lineage_cells")
  grid
}
summarize_library_coverage <- function(grid) {
  check_columns(grid, c("Lib", "Day", "Group", "denominator", "status"), "Library denominator grid")
  coverage <- unique(grid[, c("Lib", "Day", "Group", "denominator", "status")])
  if (anyDuplicated(coverage$Lib)) stop("Inconsistent denominator coverage within a library.")
  if (!nrow(coverage)) stop("No library coverage to summarize.")
  keys <- unique(coverage[, c("Group", "Day")]); rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    tab <- coverage[coverage$Group == keys$Group[i] & coverage$Day == keys$Day[i], , drop = FALSE]
    valid <- is.finite(tab$denominator) & tab$denominator > 0 & tab$status == "included"
    rows[[i]] <- data.frame(Group = keys$Group[i], Day = keys$Day[i], n_reference_libraries = nrow(tab),
      n_libraries = sum(valid), n_unavailable = sum(!valid),
      included_Lib = paste(sort(tab$Lib[valid]), collapse = ";"),
      unavailable_Lib = paste(sort(tab$Lib[!valid]), collapse = ";"),
      statistical_unit = "independent library/pool", test = "two-sided unpaired Student t-test; var.equal=TRUE")
  }
  result <- do.call(rbind, rows)
  result[order(result$Group, match(result$Day, day_order)), , drop = FALSE]
}
with_stage <- function(stage, code) {
  if (!is.character(stage) || length(stage) != 1L || is.na(stage) || !nzchar(stage)) stop("Stage must be one nonempty string.")
  # Numbered entry points use 03/04/...; consumers and release checks use 3/4/....
  if (grepl("^[0-9]+$", stage)) stage <- sub("^0+(?=[0-9])", "", stage, perl = TRUE)
  set.seed(cfg$seed)
  ensure_dir(project_path("results"))
  status_path <- file.path(outdir("audit"), paste0(stage, ".status.rds"))
  saveRDS(list(stage = stage, status = "running", time = Sys.time()), status_path)
  result <- tryCatch({
    value <- force(code)
    audit_run(stage, settings = cfg)
    saveRDS(list(stage = stage, status = "success", time = Sys.time()), status_path)
    value
  }, error = function(e) {
    saveRDS(list(stage = stage, status = "failed", time = Sys.time(), message = conditionMessage(e)), status_path)
    stop(e)
  })
  invisible(result)
}
