palette_for <- function(labels) {
  all_colors <- c(subtype_colors, level2_colors)
  colors <- all_colors[match(labels, names(all_colors))]
  missing <- is.na(colors)
  if (any(missing)) colors[missing] <- grDevices::hcl.colors(sum(missing), "Dark 3")
  setNames(colors, labels)
}
normalize_marker_result <- function(markers) {
  if (!is.data.frame(markers)) stop("Marker result must be a data frame.")
  if (!nrow(markers)) {
    return(data.frame(gene = character(), cluster = character(), p_val = numeric(),
      avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric()))
  }
  check_columns(markers, c("gene", "cluster", "p_val_adj", "avg_log2FC"), "Marker result")
  markers$gene <- as.character(markers$gene)
  markers$cluster <- as.character(markers$cluster)
  markers
}
run_markers_go <- function(obj, directory, label = "level3") {
  need(c("Seurat", "dplyr", "clusterProfiler", "org.Mm.eg.db"))
  Seurat::Idents(obj) <- obj@meta.data[[label]]
  if (length(unique(as.character(Seurat::Idents(obj)))) < 2L) stop("At least two groups are needed for marker analysis.")
  markers <- Seurat::FindAllMarkers(obj, test.use = "wilcox", only.pos = TRUE, min.pct = 0,
                                     logfc.threshold = 0, random.seed = cfg$seed)
  markers <- normalize_marker_result(markers)
  if (!nrow(markers)) writeLines("FindAllMarkers returned no markers; empty result tables are intentional.", file.path(directory, "markers_status.txt"))
  write_csv(markers, file.path(directory, "markers_all.csv"))
  significant <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.25, ]
  write_csv(significant, file.path(directory, "markers_padj0.05_log2FC0.25.csv"))
  universe <- rownames(obj)
  go <- enrich_gene_groups(split(significant$gene, significant$cluster), directory, universe = NULL)
  # No explicit universe was supplied in the source subtype GO pipeline; preserve that convention.
  top <- significant |>
    dplyr::group_by(cluster) |>
    dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC), .by_group = TRUE) |>
    dplyr::slice_head(n = 10) |>
    dplyr::ungroup()
  write_csv(top, file.path(directory, "top10_markers.csv"))
  if (nrow(top)) {
    genes <- unique(top$gene)
    save_plot(Seurat::DotPlot(obj, features = genes, group.by = label) + Seurat::RotatedAxis(),
               file.path(directory, "top10_marker_dotplot.pdf"), max(9, length(genes) * 0.15), 6)
  }
  invisible(list(markers = significant, go = go))
}
enrich_gene_groups <- function(groups, directory, universe = NULL) {
  need(c("clusterProfiler", "org.Mm.eg.db", "ggplot2", "dplyr"))
  results <- list()
  for (group in names(groups)) {
    symbols <- unique(groups[[group]])
    symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
    if (!length(symbols)) next
    ids <- suppressMessages(clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID",
                                                   OrgDb = org.Mm.eg.db::org.Mm.eg.db))
    if (!nrow(ids)) next
    args <- list(gene = unique(ids$ENTREZID), OrgDb = org.Mm.eg.db::org.Mm.eg.db, ont = "BP",
                 pAdjustMethod = "BH", pvalueCutoff = cfg$go_p, qvalueCutoff = cfg$go_q, readable = TRUE)
    if (!is.null(universe)) {
      bg <- suppressMessages(clusterProfiler::bitr(unique(universe), fromType = "SYMBOL", toType = "ENTREZID",
                                                    OrgDb = org.Mm.eg.db::org.Mm.eg.db))
      args$universe <- unique(bg$ENTREZID)
    }
    fit <- do.call(clusterProfiler::enrichGO, args)
    if (is.null(fit) || !nrow(as.data.frame(fit))) next
    tab <- as.data.frame(fit); tab$group <- group; results[[group]] <- tab
  }
  if (!length(results)) {
    write_csv(data.frame(status = "No enriched GO BP terms at configured thresholds"), file.path(directory, "GO_status.csv"))
    return(data.frame())
  }
  tab <- do.call(rbind, results)
  write_csv(tab, file.path(directory, "GO_BP_all.csv"))
  top <- tab |> dplyr::group_by(group) |> dplyr::slice_min(p.adjust, n = 3, with_ties = FALSE) |> dplyr::ungroup()
  representative <- tab |> dplyr::group_by(group) |> dplyr::slice_min(p.adjust, n = 1, with_ties = FALSE) |> dplyr::ungroup()
  write_csv(top, file.path(directory, "GO_BP_top3.csv"))
  write_csv(representative, file.path(directory, "GO_BP_representative_min_adjusted_p.csv"))
  top$term <- paste(top$group, top$Description, sep = ": ")
  top$minus_log10_p <- -log10(pmax(top$p.adjust, .Machine$double.xmin))
  p <- ggplot2::ggplot(top, ggplot2::aes(x = minus_log10_p, y = reorder(term, minus_log10_p), color = group)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = minus_log10_p, yend = term), color = "grey70") +
    ggplot2::geom_point(ggplot2::aes(size = Count)) + ggplot2::theme_bw() + ggplot2::labs(x = "-log10(adjusted p)", y = NULL)
  save_plot(p, file.path(directory, "GO_BP_top3_lollipop.pdf"), 12, max(5, nrow(top) * 0.2))
  run_go_bubble_panel(tab, directory)
  tab
}
run_source_gene_panels <- function(obj, selections, directory) {
  need(c("Seurat", "dplyr", "tidyr", "ggplot2"))
  for (panel in unique(selections$panel_id)) {
    genes <- unique(selections$gene[selections$panel_id == panel])
    present <- intersect(genes, rownames(obj))
    write_csv(data.frame(gene = genes, present = genes %in% present), file.path(directory, paste0(panel, "_gene_availability.csv")))
    if (!length(present)) next
    tab <- Seurat::FetchData(obj, vars = c("Day", "Lib", "level3", present)) |>
      tidyr::pivot_longer(dplyr::all_of(present), names_to = "gene", values_to = "expression")
    pooled <- tab |> dplyr::group_by(Day, gene) |>
      dplyr::summarise(mean = mean(expression), cell_SE = stats::sd(expression) / sqrt(dplyr::n()), .groups = "drop")
    by_subtype <- tab |> dplyr::group_by(Day, level3, gene) |> dplyr::summarise(mean = mean(expression), .groups = "drop")
    write_csv(pooled, file.path(directory, paste0(panel, "_expression_pooled.csv")))
    write_csv(by_subtype, file.path(directory, paste0(panel, "_expression_by_subtype.csv")))
    run_manual_gene_subtype_panel(by_subtype, unique(as.character(selections$key)), panel, directory)
    p <- ggplot2::ggplot(pooled, ggplot2::aes(factor(Day, levels = day_order), mean, group = 1)) +
      ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - cell_SE, ymax = mean + cell_SE), width = 0.2) +
      ggplot2::facet_wrap(~gene, scales = "free_y", ncol = 4) + ggplot2::theme_bw() +
      ggplot2::labs(x = NULL, y = "Mean normalized expression", subtitle = "Source-selected genes; descriptive cell-level SE, not independent-library inference")
    save_plot(p, file.path(directory, paste0(panel, "_expression_trends.pdf")), 13, max(5, ceiling(length(present) / 4) * 2.5))
    p <- ggplot2::ggplot(by_subtype, ggplot2::aes(factor(Day, levels = day_order), mean, color = level3, group = level3)) +
      ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~gene, scales = "free_y", ncol = 4) +
      ggplot2::scale_color_manual(values = palette_for(sort(unique(obj$level3)))) + ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Mean normalized expression")
    save_plot(p, file.path(directory, paste0(panel, "_subtype_expression_trends.pdf")), 14, max(5, ceiling(length(present) / 4) * 2.5))
  }
}
proportion_analysis <- function(obj, label, directory, allowed = NULL) {
  need(c("dplyr", "tidyr", "ggplot2"))
  meta <- obj@meta.data; check_columns(meta, c("Lib", "Day", "Group", label), "Proportion metadata")
  meta$cell_type <- as.character(meta[[label]])
  if (is.null(allowed)) allowed <- sort(unique(meta$cell_type))
  reference <- obj@misc$tbi_library_manifest
  if (!is.null(reference)) {
    tissues <- obj@misc$tbi_scope_tissues
    if (is.null(tissues)) tissues <- unique(meta$Group)
    reference <- reference[reference$Group %in% tissues, , drop = FALSE]
  }
  grid <- build_library_proportions(meta, label, allowed, reference)
  libraries <- unique(grid[, c("Lib", "Day", "Group", "Source", "Mice_per_library")])
  write_csv(unique(grid[, c("Lib", "Day", "Group", "denominator", "status")]), file.path(directory, "library_denominator_coverage.csv"))
  write_csv(summarize_library_coverage(grid), file.path(directory, "figure_caption_library_n.csv"))
  pooled_counts <- meta |> dplyr::count(Group, Day, cell_type, name = "n")
  pooled_totals <- meta |> dplyr::count(Group, Day, name = "denominator")
  pooled <- merge(unique(libraries[, c("Group", "Day")]), data.frame(cell_type = allowed), all = TRUE) |>
    dplyr::left_join(pooled_totals, by = c("Group", "Day")) |>
    dplyr::left_join(pooled_counts, by = c("Group", "Day", "cell_type")) |>
    dplyr::mutate(n = tidyr::replace_na(n, 0L), denominator = tidyr::replace_na(denominator, 0L), percentage = ifelse(denominator > 0, 100 * n / denominator, NA_real_))
  write_csv(grid, file.path(directory, "proportions_per_library_with_zeros.csv"))
  write_csv(pooled, file.path(directory, "proportions_pooled.csv"))
  summary <- grid |> dplyr::group_by(Group, Day, cell_type) |>
    dplyr::summarise(n_libraries = sum(is.finite(percentage)), n_unavailable = sum(!is.finite(percentage)),
                     mean = if (any(is.finite(percentage))) mean(percentage, na.rm = TRUE) else NA_real_,
                     sem = if (sum(is.finite(percentage)) > 1) stats::sd(percentage, na.rm = TRUE) / sqrt(sum(is.finite(percentage))) else NA_real_, .groups = "drop")
  tests <- list()
  for (tissue in intersect("Brain", unique(grid$Group))) for (ct in allowed) {
    tab <- grid[grid$Group == tissue & grid$cell_type == ct, ]
    days <- intersect(day_order, unique(tab$Day))
    if (length(days) < 2L || !"D0" %in% days) next
    for (pair in lapply(setdiff(days, "D0"), function(day) c("D0", day))) {
      x <- tab$percentage[tab$Day == pair[1]]; y <- tab$percentage[tab$Day == pair[2]]
      x <- x[is.finite(x)]; y <- y[is.finite(y)]
      status <- if (length(x) < 2 || length(y) < 2) "insufficient_independent_libraries" else "tested"
      fit <- if (status == "tested") tryCatch(stats::t.test(x, y, alternative = "two.sided", var.equal = TRUE), error = function(e) e) else NULL
      if (inherits(fit, "error")) status <- paste0("not_testable: ", conditionMessage(fit))
      tests[[length(tests) + 1L]] <- data.frame(Group = tissue, cell_type = ct, day1 = pair[1], day2 = pair[2],
        n1 = length(x), n2 = length(y), p_raw = if (status == "tested") fit$p.value else NA_real_, status = status)
    }
  }
  tests <- if (length(tests)) do.call(rbind, tests) else data.frame(Group = character(), cell_type = character(), day1 = character(), day2 = character(), n1 = integer(), n2 = integer(), p_raw = numeric(), status = character(), p_BH = numeric())
  if (nrow(tests)) {
    tests <- tests |> dplyr::group_by(Group) |> dplyr::mutate(p_BH = stats::p.adjust(p_raw, method = "BH")) |> dplyr::ungroup()
    write_csv(tests, file.path(directory, "proportion_Student_tests_raw_and_BH.csv"))
  }
  if (!nrow(tests)) write_csv(tests, file.path(directory, "proportion_Student_tests_raw_and_BH.csv"))
  summary$Day <- factor(summary$Day, levels = day_order)
  write_csv(summary, file.path(directory, "proportion_summary_mean_sem_n.csv"))
  summary$cell_type <- factor(summary$cell_type, levels = allowed)
  pooled$cell_type <- factor(pooled$cell_type, levels = allowed)
  p <- ggplot2::ggplot(summary[is.finite(summary$mean), ], ggplot2::aes(Day, mean, fill = cell_type)) +
    ggplot2::geom_col(width = 0.8) + ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, mean - sem), ymax = mean + sem), width = 0.2, na.rm = TRUE) +
    ggplot2::facet_grid(Group ~ cell_type, scales = "free_y") + ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = palette_for(allowed)) + ggplot2::labs(x = NULL, y = "Mean per-library proportion (%)") +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (nrow(tests)) {
    annotations <- tests[tests$day1 == "D0" & is.finite(tests[[cfg$proportion_label_p]]), ]
    if (nrow(annotations)) {
      annotations$p_label <- annotations[[cfg$proportion_label_p]]
      annotations$star <- ifelse(annotations$p_label < 0.001, "***", ifelse(annotations$p_label < 0.01, "**", ifelse(annotations$p_label < 0.05, "*", "ns")))
      heights <- summary |> dplyr::group_by(Group, cell_type) |>
        dplyr::summarise(height = if (any(is.finite(mean))) max(mean + ifelse(is.na(sem), 0, sem), na.rm = TRUE) * 1.12 + 0.1 else 0.1, .groups = "drop")
      annotations <- dplyr::left_join(annotations, heights, by = c("Group", "cell_type"))
      p <- p + ggplot2::geom_text(data = annotations, ggplot2::aes(x = factor(day2, levels = day_order), y = height, label = star), inherit.aes = FALSE, size = 3) +
        ggplot2::labs(subtitle = "Stars: two-sided Student t-test with D0 (equal variances); library/pool is the independent unit")
    }
  }
  save_plot(p, file.path(directory, "proportions_library_mean_sem.pdf"), max(10, length(allowed) * 1.4), 7)
  pooled$Day <- factor(pooled$Day, levels = day_order)
  p2 <- ggplot2::ggplot(pooled, ggplot2::aes(Day, percentage, fill = cell_type)) + ggplot2::geom_col() +
    ggplot2::facet_wrap(~Group, scales = "free_x") + ggplot2::scale_fill_manual(values = palette_for(allowed)) +
    ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Pooled-cell proportion (%)")
  save_plot(p2, file.path(directory, "proportions_pooled_stacked.pdf"), 9, 7)
  p3 <- ggplot2::ggplot(pooled, ggplot2::aes(Day, percentage, color = cell_type, group = cell_type)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~Group, scales = "free_x") +
    ggplot2::scale_color_manual(values = palette_for(allowed)) + ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Pooled-cell proportion (%)")
  save_plot(p3, file.path(directory, "proportions_pooled_temporal_lines.pdf"), 11, 7)
  tissue <- meta |> dplyr::count(Group, cell_type, name = "n") |> dplyr::group_by(Group) |>
    dplyr::mutate(percentage = 100 * n / sum(n)) |> dplyr::ungroup()
  p4 <- ggplot2::ggplot(tissue, ggplot2::aes(Group, percentage, fill = cell_type)) + ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(values = palette_for(allowed)) + ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = "Pooled-cell proportion (%)")
  save_plot(p4, file.path(directory, "Brain_Skull_overall_subtype_composition.pdf"), 7, 7)
  invisible(list(per_library = grid, pooled = pooled, tests = tests))
}
run_atlas <- function() {
  need(c("Seurat", "ggplot2", "dplyr", "openxlsx"))
  obj <- get_main(); d <- outdir("Figure1_S1")
  source(project_path("functions", "figure1_panels.R"), local = TRUE, encoding = "UTF-8")
  run_paper_atlas_panels(obj)
  tab <- obj@meta.data |> dplyr::count(level1, level2, Group, Day, name = "cells")
  write_csv(tab, file.path(d, "TableS1_counts_by_tissue_day.csv"))
  overall <- obj@meta.data |> dplyr::count(level1, level2, name = "cells") |> dplyr::mutate(percentage = 100 * cells / sum(cells))
  write_csv(overall, file.path(d, "TableS1_overall.csv"))
  openxlsx::write.xlsx(list(Overall = overall, Tissue_Day = tab), file.path(outdir("Supplementary_tables"), "TableS1.xlsx"), overwrite = TRUE)
  # Generic library-mean plots remain separate from the manuscript's pooled stacks.
  proportion_analysis(obj, "level2", d, level2_order)
  proportion_analysis(obj, "level1", ensure_dir(file.path(d, "broad_class_proportions")), c("Myeloid", "Lymphoid", "Non-immune"))
}
run_lineage_panels <- function(lineage, figure, tissue = NULL) {
  need(c("Seurat", "ggplot2")); obj <- get_lineage(lineage, tissue); d <- outdir(figure, lineage)
  order <- intersect(subtype_orders[[lineage]], unique(obj$level3))
  obj$level3 <- factor(obj$level3, levels = order)
  save_plot(Seurat::DimPlot(obj, group.by = "level3", label = TRUE, cols = palette_for(order), raster = TRUE), file.path(d, "subtype_UMAP.pdf"), 10, 8)
  save_plot(Seurat::DimPlot(obj, group.by = "level3", split.by = "Day", cols = palette_for(order), raster = TRUE), file.path(d, "subtype_UMAP_by_Day.pdf"), 16, 6)
  proportion_analysis(obj, "level3", d, order)
  run_markers_go(obj, d)
  run_source_lineage_panels(obj, lineage, d)
  if (lineage == "Macrophage") run_s100_cross_compartment_panel()
  if (lineage == "Microglia" && "Igf1" %in% rownames(obj)) {
    save_plot(Seurat::FeaturePlot(obj, features = "Igf1", split.by = "Day", raster = TRUE), file.path(d, "Igf1_expression_UMAP_by_day.pdf"), 15, 5)
    save_plot(Seurat::VlnPlot(obj, features = "Igf1", group.by = "level3", pt.size = 0), file.path(d, "Igf1_expression_by_microglia_state.pdf"), 10, 5)
  }
  if (lineage == "Erythroid") {
    genes <- intersect(c("Aqp1", "Hemgn", "Alas2", "Rsad2", "Hba-a1", "Hbb-bs"), rownames(obj))
    if (length(genes)) save_plot(Seurat::FeaturePlot(obj, features = genes, ncol = 3, raster = TRUE), file.path(d, "erythroid_selected_marker_featureplots.pdf"), 13, 8)
  }
}

# This module defines functions only; no model or real input is read at load time.
source(project_path("functions", "paper_panels.R"), local = TRUE, encoding = "UTF-8")
run_lymphoid_panels <- function() {
  run_lineage_panels("Bcell", "Figure5_S7"); run_lineage_panels("Tcell", "Figure5_S7")
  obj <- get_main(); cells <- colnames(obj)[obj$level2 %in% c("B Cell", "T Cell")]
  obj <- subset(obj, cells = cells)
  genes <- intersect(c("Cd79a", "Cd19", "Ube2c", "Ung", "Vpreb1", "Vpreb2", "Vpreb3", "Rag1", "Cd74", "Ms4a1", "Fcrl5", "Ly6a", "Jchain", "Mzb1", "Cd3d", "Cd3e", "Cd4", "Cd8a", "Nkg7", "Trdc", "Isg15"), rownames(obj))
  save_plot(Seurat::DotPlot(obj, features = genes, group.by = "level3") + Seurat::RotatedAxis(),
              file.path(outdir("Figure5_S7"), "B_T_combined_marker_dotplot.pdf"), 13, 9)
}
run_igf1_panels <- function() {
  need(c("Seurat", "dplyr", "ggplot2", "pheatmap"))
  obj <- get_main()
  if (!"Igf1" %in% rownames(obj)) stop("Igf1 is missing; S5E/F and S6E cannot be generated.")
  cells <- colnames(obj)[obj$level2 %in% c("Microglia", "Macrophage", "Neutrophil")]
  obj <- subset(obj, cells = cells)
  values <- Seurat::FetchData(obj, vars = c("Group", "Day", "Lib", "level2", "level3", "Igf1"))
  means <- values |> dplyr::group_by(Group, Day, level2, level3) |>
    dplyr::summarise(mean_expression = mean(Igf1), n_cells = dplyr::n(), .groups = "drop")
  libs <- values |> dplyr::group_by(Group, Day, Lib, level2, level3) |>
    dplyr::summarise(mean_expression = mean(Igf1), n_cells = dplyr::n(), .groups = "drop")
  d <- outdir("Figure4_S5", "Igf1")
  write_csv(means, file.path(d, "Igf1_pooled_mean_expression.csv"))
  write_csv(libs, file.path(d, "Igf1_expression_per_library.csv"))
  for (lineage in c("Microglia", "Macrophage")) {
    tab <- means[means$level2 == lineage & (lineage != "Microglia" | means$Group == "Brain"), ]
    tab$Day <- factor(tab$Day, levels = day_order)
    tab$level3 <- factor(tab$level3, levels = subtype_orders[[lineage]])
    p <- ggplot2::ggplot(tab, ggplot2::aes(Day, mean_expression, color = level3, group = level3)) +
      ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~Group, scales = "free_x") +
      ggplot2::scale_color_manual(values = palette_for(subtype_orders[[lineage]])) + ggplot2::theme_bw() +
      ggplot2::labs(x = NULL, y = "Mean normalized Igf1 expression")
    save_plot(p, file.path(d, paste0(if (lineage == "Microglia") "S5E" else "S5F", "_", lineage, "_Igf1_temporal.pdf")), 11, 6)
  }
  selected <- c("Microglia3", "Microglia7", "Macrophage2", "Macrophage5", "Neutrophil1")
  tab <- means[means$Group == "Brain" & means$level3 %in% selected, ]
  mat <- matrix(NA_real_, length(selected), length(day_order), dimnames = list(selected, day_order))
  mat[cbind(match(tab$level3, selected), match(tab$Day, day_order))] <- tab$mean_expression
  write_csv(data.frame(level3 = rownames(mat), mat, check.names = FALSE), file.path(d, "S6E_Igf1_raw_mean_matrix.csv"))
  if (!any(is.finite(mat))) stop("No selected S6E Igf1 profiles are available.")
  pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", na_col = "grey85",
    filename = file.path(d, "S6E_Igf1_raw_mean_heatmap.pdf"), width = 7, height = 5,
    main = "Igf1 mean normalized expression (no row scaling)")
}
run_reference <- function() {
  need(c("Seurat", "dplyr", "pheatmap")); obj <- get_main(); obj <- subset(obj, subset = Group == "Skull")
  ref <- read_input(cfg$input_reference)
  Seurat::DefaultAssay(ref) <- "RNA"
  check_columns(ref@meta.data, cfg$reference_label, "GSE192616 reference")
  # Use the merged 16-cell-type labels established in the source reference object.
  types <- c("granulocyte-monocyte progenitor", "monocyte progenitor", "monocyte", "immature neutrophil", "mature neutrophil",
    "basophil", "CDP", "pDC", "cDC", "pro B", "pre B", "B cell", "T cell", "NK cell", "erythroblast", "erythrocyte")
  labels <- as.character(ref@meta.data[[cfg$reference_label]])
  aliases <- c("common DC progenitor (CDP)" = "CDP", "plasmacytoid DC" = "pDC", "pro B cell" = "pro B", "pre B cell" = "pre B")
  hit <- labels %in% names(aliases); labels[hit] <- aliases[labels[hit]]
  ref$reference_type <- labels
  missing <- setdiff(types, labels); if (length(missing)) stop("Reference must have the final 16 merged types. Missing: ", paste(missing, collapse = ", "))
  ref_display <- ref  # S1C shows the complete annotated reference; S1D uses the 16 matched types.
  ref <- subset(ref, cells = colnames(ref)[labels %in% types])
  d <- outdir("FigureS1_reference")
  Seurat::Idents(obj) <- obj$level2
  marker_path <- input_path(cfg$reference_markers)
  markers <- if (file.exists(marker_path)) read_csv(marker_path) else Seurat::FindAllMarkers(obj, only.pos = TRUE, min.pct = 0, logfc.threshold = 0, test.use = "wilcox", random.seed = cfg$seed)
  check_columns(markers, c("gene", "cluster", "avg_log2FC"), "Reference marker input")
  top <- markers |> dplyr::group_by(cluster) |> dplyr::slice_max(avg_log2FC, n = 200, with_ties = FALSE) |> dplyr::ungroup()
  a <- Seurat::AverageExpression(obj, group.by = "level2", assays = "RNA", slot = "data")$RNA
  b <- Seurat::AverageExpression(ref, group.by = "reference_type", assays = "RNA", slot = "data")$RNA
  genes <- Reduce(intersect, list(unique(top$gene), rownames(a), rownames(b)))
  if (length(genes) < 10) stop("Insufficient shared reference genes.")
  a <- a[genes, intersect(c("Promyelocyte", "Macrophage", "Neutrophil", "Basophil", "DC", "B Cell", "T Cell", "Erythrocyte"), colnames(a)), drop = FALSE]
  b <- b[genes, types, drop = FALSE]
  correlation <- stats::cor(as.matrix(b), as.matrix(a), method = "spearman")
  if (anyNA(correlation)) stop("Undefined correlations; inspect constant profiles.")
  if (any(apply(correlation, 2, stats::sd) == 0)) stop("A reference correlation column is constant and cannot be column-scaled.")
  utils::write.csv(correlation, file.path(d, "Spearman_raw_coefficients.csv"))
  write_csv(data.frame(gene = genes), file.path(d, "correlation_genes.csv"))
  colors <- grDevices::colorRampPalette(c("#14448C", "#5076C1", "#92A6DE", "#CED7F2", "#F6F6F6", "#DE9494", "#B85A5B", "#7D2828"))(100)
  pheatmap::pheatmap(correlation, cluster_rows = FALSE, cluster_cols = FALSE, scale = "column", color = colors,
                      main = "Skull reference: Spearman correlation (column-scaled display)", filename = file.path(d, "S1D_Spearman_column_scaled.pdf"), width = 10, height = 10)
  save_plot(Seurat::DimPlot(ref_display, group.by = "reference_type", label = TRUE, raster = TRUE), file.path(d, "S1C_reference_UMAP.pdf"), 11, 9)
}
