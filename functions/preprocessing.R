run_qc <- function() {
  need(c("Seurat", "DoubletFinder"))
  manifest_path <- input_path(cfg$sample_manifest)
  m <- read_csv(manifest_path)
  check_columns(m, c("Lib", "matrix_dir", "gene_column", "Day", "Group", "Source"), "sample_manifest")
  if (nrow(m) != cfg$expected_libraries || anyDuplicated(m$Lib)) stop("Expected 31 unique library rows.")
  if (sum(m$Source == "Public") != cfg$expected_public_libraries) stop("Expected 25 public libraries.")
  if (anyNA(m[, c("Lib", "matrix_dir", "Day", "Group", "Source", "gene_column")]) || any(!nzchar(trimws(m$Lib))) || any(!nzchar(trimws(m$matrix_dir))) || any(!m$Source %in% c("Public", "In-house")) || any(!m$gene_column %in% c(1L, 2L)) || any(!m$Day %in% day_order) ||
      any(!m$Group %in% c("Brain", "Skull"))) stop("Complete the manifest from actual metadata before QC.")
  get_df_function <- function(old, new) {
    ns <- asNamespace("DoubletFinder")
    if (exists(old, envir = ns, inherits = FALSE)) get(old, ns) else get(new, ns)
  }
  sweep_fun <- get_df_function("paramSweep_v3", "paramSweep")
  doublet_fun <- get_df_function("doubletFinder_v3", "doubletFinder")
  objects <- vector("list", nrow(m)); names(objects) <- m$Lib
  counts <- list()
  for (i in seq_along(objects)) {
    matrix_dir <- input_path(m$matrix_dir[i])
    files <- list.files(matrix_dir, full.names = TRUE)
    pick <- function(pattern) {
      hit <- files[grepl(pattern, basename(files))]
      if (length(hit) != 1L) stop("Expected exactly one ", pattern, " file in ", matrix_dir)
      hit
    }
    # Read prefixed GEO files directly; never rename the raw download in place.
    counts_in <- Seurat::ReadMtx(mtx = pick("matrix\\.mtx(\\.gz)?$"),
      cells = pick("barcodes\\.tsv(\\.gz)?$"), features = pick("(features|genes)\\.tsv(\\.gz)?$"),
      feature.column = as.integer(m$gene_column[i]))
    if (is.list(counts_in)) counts_in <- counts_in[["Gene Expression"]]
    if (is.null(counts_in)) stop("No gene expression matrix: ", m$Lib[i])
    obj <- Seurat::CreateSeuratObject(counts_in, min.cells = 1L, project = m$Lib[i])
    obj <- Seurat::RenameCells(obj, add.cell.id = m$Lib[i])
    obj$Lib <- m$Lib[i]; obj$Day <- m$Day[i]; obj$Group <- m$Group[i]
    obj$Source <- m$Source[i]; obj$sample_batch_library <- m$Lib[i]
    obj$percent.mt <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-|^MT-|^Mt-")
    before <- ncol(obj)
    # Preserve the original sorted-position rule rather than changing to interpolated quantile.
    index <- max(1L, floor(cfg$upper_feature_fraction * before))
    upper <- sort(obj$nFeature_RNA)[index]
    obj <- subset(obj, subset = nFeature_RNA > cfg$min_feature & nFeature_RNA < upper & percent.mt < cfg$max_mt)
    after_qc <- ncol(obj)
    if (after_qc < 20L) stop("Too few cells after QC in ", m$Lib[i])
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = cfg$nfeatures, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = min(cfg$npcs, ncol(obj) - 1L), seed.use = cfg$seed, verbose = FALSE)
    set.seed(cfg$seed)
    sweep <- sweep_fun(obj, PCs = cfg$doublet_pcs, sct = FALSE)
    stats <- DoubletFinder::summarizeSweep(sweep, GT = FALSE)
    bcm <- DoubletFinder::find.pK(stats)
    pK <- as.numeric(as.character(bcm$pK[which.max(bcm$BCmetric)]))
    nExp <- round(8e-6 * ncol(obj)^2)
    if (nExp < 1L || nExp >= ncol(obj)) stop("Invalid expected doublet count in ", m$Lib[i])
    set.seed(cfg$seed)
    obj <- doublet_fun(obj, PCs = cfg$doublet_pcs, pN = cfg$doublet_pN, pK = pK,
                       nExp = nExp, reuse.pANN = FALSE, sct = FALSE)
    label <- grep("^DF.classifications", names(obj@meta.data), value = TRUE)
    if (length(label) != 1L) stop("Ambiguous DoubletFinder classification column.")
    obj$doublet_info <- obj@meta.data[[label]]
    obj <- subset(obj, subset = doublet_info == "Singlet")
    objects[[i]] <- obj
    counts[[i]] <- data.frame(Lib = m$Lib[i], before = before, after_qc = after_qc,
                              after_doublets = ncol(obj), upper_feature = upper, pK = pK, nExp = nExp)
    saveRDS(obj, file.path(outdir("qc", "libraries"), paste0(m$Lib[i], ".rds")))
  }
  write_csv(do.call(rbind, counts), file.path(outdir("qc"), "library_cell_counts.csv"))
  saveRDS(objects, file.path(outdir("objects"), "qc_libraries.rds"))
  audit_run("01_QC", manifest_path, cfg)
}
run_integration <- function() {
  need(c("Seurat", "harmony"))
  objects <- read_input(project_path("results", "objects", "qc_libraries.rds"))
  if (length(objects) != cfg$expected_libraries) stop("QC object list does not contain 31 libraries.")
  obj <- merge(objects[[1]], y = objects[-1])
  if (anyDuplicated(colnames(obj))) stop("Cell IDs must be globally unique.")
  obj <- Seurat::NormalizeData(obj, scale.factor = 10000, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = cfg$nfeatures, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = cfg$npcs, seed.use = cfg$seed, verbose = FALSE)
  set.seed(cfg$seed)
  obj <- harmony::RunHarmony(obj, group.by.vars = "Lib")
  obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = cfg$dims)
  obj <- Seurat::FindClusters(obj, resolution = cfg$resolution, random.seed = cfg$seed)
  obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = cfg$dims, seed.use = cfg$seed)
  saveRDS(obj, file.path(outdir("objects"), "integrated_unannotated.rds"))
  # Recomputed cluster numbers are not guaranteed to match the historical 0-33 map.
  features <- intersect(unique(unlist(marker_level2)), rownames(obj))
  save_plot(Seurat::DotPlot(obj, features = features) + Seurat::RotatedAxis(),
            file.path(outdir("integration"), "review_cluster_markers.pdf"), 15, 8)
  mapfile <- project_path("config", "level2_cluster_map.csv")
  map <- read_csv(mapfile)
  check_columns(map, c("cluster", "level2"), "level2 cluster map")
  if (anyNA(map[, c("cluster", "level2")]) || anyDuplicated(map$cluster) || any(!map$level2 %in% c(level2_order, "DROP"))) stop("Invalid or duplicated level2 cluster map.")
  if (!isTRUE(cfg$accept_recomputed_cluster_map)) {
    stop("Integration saved. Review marker plot and level2_cluster_map.csv, then set accept_recomputed_cluster_map=TRUE in settings.R.")
  }
  ids <- as.character(Seurat::Idents(obj))
  if (any(!ids %in% as.character(map$cluster))) stop("Missing new clusters in the approved level2 map.")
  labels <- map$level2[match(ids, as.character(map$cluster))]
  keep <- !is.na(labels) & labels != "DROP"
  obj <- subset(obj, cells = colnames(obj)[keep]); obj$level2 <- labels[keep]
  obj <- normalize_metadata(obj)
  saveRDS(obj, file.path(outdir("objects"), "integrated_level2.rds"))
}
clean_subtype <- function(obj, lineage, label_column = "level3", state = "auto") {
  if (!lineage %in% names(subtype_orders)) stop("Unknown subtype lineage: ", lineage)
  if (length(state) != 1L || is.na(state) || !state %in% c("auto", "legacy", "final", "raw_clusters")) stop("Invalid annotation_state for ", lineage)
  if (!label_column %in% names(obj@meta.data)) stop(lineage, ": missing label column ", label_column)
  if (state == "raw_clusters") {
    if (!isTRUE(cfg$accept_subtype_cluster_maps)) stop("Review subtype markers/maps before setting accept_subtype_cluster_maps=TRUE.")
    mapping <- read_csv(project_path("config", "subtype_cluster_maps.csv"))
    mapping <- mapping[mapping$lineage == lineage, ]
    if (!nrow(mapping) || anyDuplicated(mapping$cluster) || anyNA(mapping$level3)) stop("Invalid or duplicated subtype cluster map: ", lineage)
    raw <- as.character(obj@meta.data[[label_column]])
    if (any(!raw %in% as.character(mapping$cluster))) stop("Unmapped raw clusters in ", lineage)
    obj$level3 <- mapping$level3[match(raw, as.character(mapping$cluster))]
    label_column <- "level3"; state <- "final"
  }
  labels <- canonical_labels(obj@meta.data[[label_column]])
  if (anyNA(labels) || any(!nzchar(trimws(labels)))) stop(lineage, ": missing labels are not approved exclusions; correct the annotation or supply explicit manual cell exclusions.")
  if (lineage == "Macrophage") {
    if (state == "auto") {
      if ("Macrophage13" %in% labels) state <- "legacy" else {
        if ("Macrophage12" %in% labels) stop("Macrophage12 is ambiguous. Set annotation_state=final or legacy in subtype_inputs.csv after checking the source object.")
        state <- "final"
      }
    }
    if (!state %in% c("legacy", "final")) stop("Macrophage annotation_state must be legacy/final/auto.")
    if (state == "legacy") {
      labels[labels == "Macrophage12"] <- NA_character_
      labels[labels == "Macrophage13"] <- "Macrophage12"
    }
  }
  if (lineage == "Neutrophil") labels[labels %in% paste0("Neutrophil", 10:12)] <- NA_character_
  if (lineage == "Tcell") labels[labels == "CD8-ODG"] <- NA_character_
  if (lineage == "Microglia") labels[labels == "Microglia9"] <- NA_character_
  # Only explicit legacy exclusions are dropped; unexpected labels are errors.
  labels[labels %in% c("Del", "DROP", "Discard")] <- NA_character_
  unexpected <- setdiff(unique(stats::na.omit(labels)), subtype_orders[[lineage]])
  if (length(unexpected)) stop(lineage, ": unexpected labels: ", paste(unexpected, collapse = ", "))
  if (any(!is.na(labels) & !nzchar(labels))) stop(lineage, ": empty subtype annotations.")
  removed <- colnames(obj)[is.na(labels)]
  if (all(is.na(labels))) stop("No retained cells remain in ", lineage)
  obj <- subset(obj, cells = colnames(obj)[!is.na(labels)])
  obj$level3 <- labels[!is.na(labels)]
  obj@meta.data$level4 <- NULL
  obj@misc$tbi_subtype <- list(lineage = lineage, cleaned = TRUE, annotation_state = state)
  list(object = obj, removed = removed)
}
run_annotation <- function() {
  need("Seurat")
  generated <- project_path("results", "objects", "integrated_level2.rds")
  if (identical(cfg$main_source, "generated")) {
    path <- generated
  } else if (identical(cfg$main_source, "configured")) {
    path <- if (file.exists(input_path(cfg$input_main))) input_path(cfg$input_main) else input_path(cfg$input_level2)
  } else stop("Set main_source to configured or generated; implicit precedence is not allowed.")
  main <- normalize_metadata(read_input(path))
  if (cfg$subtype_mode == "recluster_review") {
    run_subtype_reclustering_review(main)
    stop("Subtype candidates and marker tables saved. Review maps, then select reviewed candidates in subtype_inputs.csv using raw_clusters mode.")
  }
  before <- colnames(main)
  manual_path <- project_path("config", "manual_excluded_cells.csv")
  manual <- read_csv(manual_path)
  check_columns(manual, c("cell_id", "reason"), "manual exclusions")
  if (anyNA(manual$cell_id) || anyDuplicated(manual$cell_id) || any(!nzchar(manual$cell_id))) stop("Manual exclusions contain missing, empty or duplicate cell IDs.")
  if (length(intersect(manual$cell_id, before))) main <- subset(main, cells = setdiff(before, manual$cell_id))
  main$level3 <- if ("level3" %in% names(main@meta.data)) canonical_labels(main$level3) else main$level2
  excluded <- list(); subobjects <- list()
  if (cfg$input_mode == "subtype_objects") {
    manifest_path <- input_path(cfg$subtype_manifest)
    manifest <- read_csv(manifest_path)
    check_columns(manifest, c("lineage", "path", "label_column", "annotation_state"), "subtype_inputs")
    if (!setequal(manifest$lineage, names(subtype_orders)) || anyDuplicated(manifest$lineage)) stop("Provide all seven unique subtype inputs.")
    for (lineage in names(subtype_orders)) {
      row <- manifest[manifest$lineage == lineage, ]
      sub <- read_input(row$path)
      cells <- colnames(sub)
      if (anyDuplicated(cells) || any(!cells %in% before)) stop(lineage, ": cell IDs do not match the main object.")
      retained_input <- setdiff(cells, manual$cell_id)
      if (!length(retained_input)) stop("Manual exclusions removed all cells in ", lineage)
      sub <- subset(sub, cells = retained_input); cells <- colnames(sub)
      parents <- lineage_parent[[lineage]]
      if (any(!main@meta.data[cells, "level2"] %in% parents)) stop(lineage, ": cells belong to the wrong parent lineage.")
      if (lineage %in% c("Erythroid", "NVU")) {
        tissue <- if (lineage == "Erythroid") "Skull" else "Brain"
        if (any(main@meta.data[cells, "Group"] != tissue)) stop(lineage, ": wrong tissue compartment.")
      }
      cleaned <- clean_subtype(sub, lineage, row$label_column, row$annotation_state)
      sub <- cleaned$object
      scope <- main$level2 %in% parents
      if (lineage == "Erythroid") scope <- scope & main$Group == "Skull"
      if (lineage == "NVU") scope <- scope & main$Group == "Brain"
      discarded <- if (lineage == "Erythroid") character() else setdiff(colnames(main)[scope], colnames(sub))
      excluded[[lineage]] <- discarded
      main <- subset(main, cells = setdiff(colnames(main), discarded))
      if (lineage == "Erythroid") {
        # Erythroblast states belong to the local Figure6 trajectory, not the
        # confirmed 58-label global CellChat/atlas annotation.
        main$level3[scope] <- "Erythrocyte"
        main$erythroid_subtype <- NA_character_
        main$erythroid_subtype[match(colnames(sub), colnames(main))] <- as.character(sub$level3)
      } else main$level3[match(colnames(sub), colnames(main))] <- as.character(sub$level3)
      subobjects[[lineage]] <- sub
    }
  } else if (cfg$input_mode == "final_main") {
    # Explicitly choosing final_main asserts that old macrophage12 was already removed.
    for (lineage in names(subtype_orders)) {
      if (lineage == "Erythroid") {
        cells <- colnames(main)[main$level2 == "Erythrocyte" & main$Group == "Skull"]
        local_labels <- if ("erythroid_subtype" %in% names(main@meta.data)) canonical_labels(main@meta.data[cells, "erythroid_subtype"]) else canonical_labels(main@meta.data[cells, "level3"])
        if (all(!is.na(local_labels) & local_labels %in% subtype_orders$Erythroid) && all(subtype_orders$Erythroid %in% local_labels)) {
          sub <- subset(main, cells = cells); sub$level3 <- local_labels
        } else {
          manifest <- read_csv(input_path(cfg$subtype_manifest))
          row <- manifest[manifest$lineage == "Erythroid", , drop = FALSE]
          if (nrow(row) != 1L) stop("Provide the local erythroid subtype object for Figure6 trajectories.")
          sub <- read_input(row$path)
          if (anyDuplicated(colnames(sub)) || any(!colnames(sub) %in% cells)) stop("Local erythroid cells must match final Skull erythrocytes by ID.")
          sub <- clean_subtype(sub, "Erythroid", row$label_column, row$annotation_state)$object
        }
        if (!all(subtype_orders$Erythroid %in% sub$level3)) stop("Local erythroid object lacks required trajectory states.")
        main$level3[main$level2 == "Erythrocyte"] <- "Erythrocyte"
        main$erythroid_subtype <- NA_character_
        main$erythroid_subtype[match(colnames(sub), colnames(main))] <- as.character(sub$level3)
        subobjects[[lineage]] <- sub
        next
      }
      scope <- main$level3 %in% subtype_orders[[lineage]] & main$level2 %in% lineage_parent[[lineage]]
      if (lineage == "Erythroid") scope <- scope & main$Group == "Skull"
      if (lineage == "NVU") scope <- scope & main$Group == "Brain"
      cells <- colnames(main)[scope]
      if (!length(cells)) stop("Final main lacks ", lineage, " subtypes.")
      subobjects[[lineage]] <- subset(main, cells = cells)
    }
    if (any(main$level3 %in% c("Macrophage13", paste0("Neutrophil", 10:13), "CD8-ODG", "Microglia9"))) stop("final_main still contains legacy unwanted subtypes.")
  } else stop("input_mode must be subtype_objects or final_main.")
  main$level3 <- global_erythroid_labels(main$level3)
  main@meta.data$level4 <- NULL
  allowed <- unique(c(global_level3_labels(), level2_order))
  if (anyNA(main$level3) || any(!main$level3 %in% allowed)) stop("Invalid final level3 values.")
  if (any(main$level3 == "Erythrocyte" & main$level2 != "Erythrocyte")) stop("Global Erythrocyte labels occur under the wrong parent.")
  for (lineage in names(subtype_orders)) {
    if (lineage == "Erythroid") next
    mask <- main$level2 %in% lineage_parent[[lineage]]
    if (lineage == "Erythroid") mask <- mask & main$Group == "Skull"
    if (lineage == "NVU") mask <- mask & main$Group == "Brain"
    if (any(!main$level3[mask] %in% subtype_orders[[lineage]])) stop("Unannotated cells remain in ", lineage)
  }
  for (lineage in names(subtype_orders)) {
    if (lineage == "Erythroid") {
      if (any(main$level3[main$level2 == "Erythrocyte"] != "Erythrocyte")) stop("Global erythroid labels must remain Erythrocyte.")
      next
    }
    wrongly_parented <- main$level3 %in% subtype_orders[[lineage]] & !main$level2 %in% lineage_parent[[lineage]]
    if (any(wrongly_parented)) stop("Subtype labels occur under the wrong level2 parent: ", lineage)
    if (lineage == "NVU" && any(main$level3 %in% subtype_orders$NVU & main$Group != "Brain")) stop("NVU subtype labels occur outside Brain.")
    if (lineage == "Erythroid" && any(main$level3 %in% c("Erythroblast1", "Erythroblast2") & main$Group != "Skull")) stop("Erythroblast subtype labels occur outside Skull.")
    scope <- main$level2 %in% lineage_parent[[lineage]]
    if (lineage == "Erythroid") scope <- scope & main$Group == "Skull"
    if (lineage == "NVU") scope <- scope & main$Group == "Brain"
    missing_states <- setdiff(subtype_orders[[lineage]], unique(main$level3[scope]))
    if (length(missing_states)) stop("Missing expected final states in ", lineage, ": ", paste(missing_states, collapse = ", "))
  }
  if (any(!main$level3 %in% global_level3_labels())) stop("Unresolved global level3 labels remain outside the confirmed 58-state vocabulary.")
  main <- normalize_metadata(main)
  main$functional_annotation <- as.character(main$level3)
  functional_names <- c(Macrophage8 = "Chemotactic Surveillance", Macrophage10 = "Respiratory Burst", Neutrophil4 = "Parenchyma-Adapted")
  hit <- main$level3 %in% names(functional_names)
  main$functional_annotation[hit] <- unname(functional_names[main$level3[hit]])
  main@misc$tbi_schema <- "TBI_final_level3_v1"
  main@misc$tbi_source_md5 <- hash_input(path)
  main@misc$tbi_deleted_ids <- setdiff(before, colnames(main))
  main@misc$tbi_library_manifest <- unique(main@meta.data[, c("Lib", "Day", "Group", "Source", "Mice_per_library", "Control")])
  ensure_dir(dirname(main_final_path()))
  saveRDS(main, main_final_path())
  main_fingerprint <- hash_input(main_final_path())
  caption_counts <- list()
  for (lineage in names(subobjects)) {
    sub <- subobjects[[lineage]]
    local_erythroid <- if (lineage == "Erythroid") as.character(sub$level3) else NULL
    for (field in c("Lib", "Day", "Group", "Source", "Control", "Mice_per_library", "level1", "level2", "level3", "functional_annotation"))
      sub@meta.data[[field]] <- main@meta.data[colnames(sub), field]
    if (lineage == "Erythroid") {
      sub$global_level3 <- as.character(sub$level3)
      sub$level3 <- local_erythroid; sub$erythroid_subtype <- local_erythroid
      sub$functional_annotation <- local_erythroid
    }
    sub@meta.data$level4 <- NULL
    sub@misc$tbi_main_md5 <- main_fingerprint
    sub@misc$tbi_library_manifest <- main@misc$tbi_library_manifest
    sub@misc$tbi_scope_tissues <- if (lineage == "Erythroid") "Skull" else if (lineage == "NVU") "Brain" else c("Brain", "Skull")
    saveRDS(sub, file.path(outdir("objects"), paste0(lineage, "_final.rds")))
    if (lineage %in% c("Macrophage", "Neutrophil", "Microglia")) {
      meta <- sub@meta.data
      # S2C and S5B are Brain panels; S3C may contain both compartments.
      if (lineage != "Neutrophil") meta <- meta[meta$Group == "Brain", , drop = FALSE]
      tissues <- if (lineage == "Neutrophil") c("Brain", "Skull") else "Brain"
      reference <- main@misc$tbi_library_manifest
      reference <- reference[reference$Group %in% tissues, , drop = FALSE]
      grid <- build_library_proportions(meta, "level3", subtype_orders[[lineage]], reference)
      n <- summarize_library_coverage(grid)
      n$panel <- c(Macrophage = "S2C", Neutrophil = "S3C", Microglia = "S5B")[[lineage]]
      n$lineage <- lineage; caption_counts[[lineage]] <- n
    }
  }
  write_csv(do.call(rbind, caption_counts), file.path(outdir("audit"), "S2C_S3C_S5B_library_n.csv"))
  write_csv(data.frame(cell_id = main@misc$tbi_deleted_ids), file.path(outdir("audit"), "deleted_cell_ids.csv"))
  write_csv(as.data.frame(table(main$Group, main$Day, main$level2)), file.path(outdir("audit"), "final_cell_counts.csv"))
  libs <- unique(main@meta.data[, c("Lib", "Day", "Group", "Source", "Mice_per_library", "Control")])
  libs <- libs[order(libs$Group, match(libs$Day, day_order), libs$Lib), ]
  write_csv(libs, file.path(outdir("audit"), "observed_library_manifest.csv"))
  write_csv(data.frame(total_cells = ncol(main), Brain = sum(main$Group == "Brain"), Skull = sum(main$Group == "Skull"),
                        level2_types = length(unique(main$level2)), level3_types = length(unique(main$level3))),
             file.path(outdir("audit"), "actual_totals_not_manuscript_constants.csv"))
  audit_run("03_final_annotation", path, cfg)
}
run_subtype_reclustering_review <- function(main) {
  need(c("Seurat", "harmony"))
  resolutions <- c(Macrophage = 0.5, Neutrophil = 0.5, Microglia = 0.3, Bcell = 0.4, Tcell = 0.2, Erythroid = 0.3, NVU = 0.3)
  exclusions <- read_csv(project_path("config", "manual_excluded_cells.csv"))
  check_columns(exclusions, c("cell_id", "reason"), "manual exclusions")
  for (lineage in names(subtype_orders)) {
    keep <- main$level2 %in% lineage_parent[[lineage]]
    if (lineage == "Erythroid") keep <- keep & main$Group == "Skull"
    if (lineage == "NVU") keep <- keep & main$Group == "Brain"
    obj <- subset(main, cells = setdiff(colnames(main)[keep], exclusions$cell_id))
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = cfg$nfeatures, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = cfg$npcs, seed.use = cfg$seed, verbose = FALSE)
    set.seed(cfg$seed)
    obj <- harmony::RunHarmony(obj, group.by.vars = "Lib", max.iter.harmony = 40)
    obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = cfg$dims)
    obj <- Seurat::FindClusters(obj, resolution = resolutions[[lineage]], random.seed = cfg$seed)
    obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = cfg$dims, seed.use = cfg$seed)
    d <- outdir("subtype_review", lineage)
    saveRDS(obj, file.path(d, "reclustered_unannotated.rds"))
    save_plot(Seurat::DimPlot(obj, label = TRUE, raster = TRUE), file.path(d, "review_UMAP_raw_clusters.pdf"), 10, 8)
    markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE, min.pct = 0, logfc.threshold = 0, test.use = "wilcox", random.seed = cfg$seed)
    write_csv(markers, file.path(d, "review_raw_cluster_markers.csv"))
    write_csv(data.frame(cluster = sort(unique(as.character(Seurat::Idents(obj))))), file.path(d, "observed_raw_clusters.csv"))
  }
}
