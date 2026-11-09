trajectory_input_stamp <- function(day = NULL) {
  if (!is.null(day) && (length(day) != 1L || is.na(day) || !day %in% c("D0", "D1", "D7"))) stop("Invalid erythroid trajectory time point.")
  list(main_md5 = hash_input(main_final_path()),
       erythroid_md5 = hash_input(project_path("results", "objects", "Erythroid_final.rds")),
       trajectory_code_md5 = hash_input(project_path("functions", "trajectory.R")),
       day = day, seed = cfg$seed, reduction = "UMAP", cluster = "level3", start = "Erythroblast1")
}
read_slingshot_pseudotime <- function(obj, day = NULL) {
  assert_no_incomplete_stages("3")
  if (is.null(day)) assert_no_incomplete_stages("27")
  suffix <- if (is.null(day)) "All" else day
  d <- project_path("results", "Figure6_S12", "Erythroid_trajectory", suffix)
  path <- file.path(d, "slingshot_pseudotime.csv")
  stamp_path <- file.path(d, "slingshot_input_stamp.rds")
  if (!file.exists(stamp_path)) stop("Slingshot input provenance is missing; rerun script ", if (is.null(day)) "27" else "29", ".")
  expected <- list(inputs = trajectory_input_stamp(day), csv_md5 = hash_input(path))
  if (!identical(readRDS(stamp_path), expected)) stop("Slingshot input or pseudotime file changed; rerun script ", if (is.null(day)) "27" else "29", ".")
  ptime <- read_csv(path)
  check_columns(ptime, c("cell_id", "slingPseudotime_1"), "Slingshot pseudotime")
  if (anyNA(ptime$cell_id) || anyDuplicated(ptime$cell_id) || !setequal(ptime$cell_id, colnames(obj))) stop("Slingshot and erythroid object cells differ or cell IDs are duplicated.")
  if (!is.numeric(ptime$slingPseudotime_1)) stop("Slingshot pseudotime must be numeric.")
  ptime
}
run_slingshot <- function(day = NULL) {
  need(c("Seurat", "slingshot", "SingleCellExperiment", "S4Vectors", "ggplot2"))
  obj <- get_lineage("Erythroid", "Skull")
  if (!is.null(day)) obj <- subset(obj, cells = colnames(obj)[obj$Day == day])
  states <- unique(as.character(obj$level3))
  if (!"Erythroblast1" %in% states || length(states) < 2L) stop("Slingshot requires Erythroblast1 and at least one other state.")
  if (!"umap" %in% Seurat::Reductions(obj)) stop("The original erythroid UMAP is missing.")
  sce <- Seurat::as.SingleCellExperiment(obj)
  sce$subtype <- as.character(obj$level3)
  SingleCellExperiment::reducedDim(sce, "UMAP") <- Seurat::Embeddings(obj, "umap")
  set.seed(cfg$seed)
  sce <- slingshot::slingshot(sce, clusterLabels = "subtype", reducedDim = "UMAP", start.clus = "Erythroblast1", end.clus = NULL)
  pseudotime <- slingshot::slingPseudotime(sce)
  if (!ncol(pseudotime)) stop("No Slingshot lineages were inferred.")
  tab <- data.frame(cell_id = rownames(pseudotime), pseudotime)
  names(tab)[-1] <- paste0("slingPseudotime_", seq_len(ncol(pseudotime)))
  suffix <- if (is.null(day)) "All" else day
  d <- outdir("Figure6_S12", "Erythroid_trajectory", suffix)
  ptime_path <- file.path(d, "slingshot_pseudotime.csv")
  write_csv(tab, ptime_path); saveRDS(sce, file.path(d, "slingshot_sce.rds"))
  saveRDS(list(inputs = trajectory_input_stamp(day), csv_md5 = hash_input(ptime_path)), file.path(d, "slingshot_input_stamp.rds"))
  if (!any(is.finite(pseudotime[, 1]))) stop("The first Slingshot lineage has no finite pseudotime.")
  with_pdf_device(file.path(d, "Slingshot_UMAP.pdf"), 7, 6, {
  colors <- grDevices::colorRampPalette(c("#313695", "#FFFFBF", "#A50026"))(100)
  time <- pseudotime[, 1]; bins <- cut(time, breaks = 100, labels = FALSE); cols <- colors[bins]; cols[is.na(cols)] <- "grey85"
  graphics::plot(SingleCellExperiment::reducedDim(sce, "UMAP"), col = cols, pch = 16, asp = 1, xlab = "UMAP1", ylab = "UMAP2")
  graphics::lines(slingshot::SlingshotDataSet(sce), lwd = 2, col = "black")
  })
  audit_run(paste0("Slingshot_", suffix), project_path("results", "objects", "Erythroid_final.rds"), list(start = "Erythroblast1", reduction = "UMAP", cluster = "level3", seed = cfg$seed))
}
run_monocle_genes <- function(day = NULL) {
  set.seed(cfg$seed)
  need(c("Seurat", "monocle", "Biobase", "VGAM", "Matrix", "ggplot2"))
  suppressPackageStartupMessages(library(monocle))
  obj <- get_lineage("Erythroid", "Skull"); if (!is.null(day)) obj <- subset(obj, cells = colnames(obj)[obj$Day == day])
  suffix <- if (is.null(day)) "All" else day; d <- outdir("Figure6_S12", "Erythroid_trajectory", suffix)
  ptime <- read_slingshot_pseudotime(obj, day)
  obj$slingPseudotime_1 <- ptime$slingPseudotime_1[match(colnames(obj), ptime$cell_id)]
  # Monocle2 gene analysis uses Slingshot pseudotime, as in the manuscript source.
  valid <- is.finite(obj$slingPseudotime_1)
  excluded_ids <- colnames(obj)[!valid]
  write_csv(data.frame(cell_id = excluded_ids, reason = rep("Not assigned to Slingshot lineage1", length(excluded_ids))), file.path(d, "cells_excluded_from_lineage1_gene_models.csv"))
  obj <- subset(obj, cells = colnames(obj)[valid])
  if (ncol(obj) < 10L) stop("Too few cells with finite lineage1 pseudotime.")
  counts <- methods::as(Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"), "sparseMatrix")
  genes <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  cds <- monocle::newCellDataSet(counts, phenoData = methods::new("AnnotatedDataFrame", data = obj@meta.data),
    featureData = methods::new("AnnotatedDataFrame", data = genes), expressionFamily = VGAM::negbinomial.size(), lowerDetectionLimit = 0.5)
  cds <- estimateSizeFactors(cds); cds <- estimateDispersions(cds)
  cds <- monocle::detectGenes(cds, min_expr = 0.1)
  expressed <- rownames(Biobase::fData(cds))[Biobase::fData(cds)$num_cells_expressed >= 10]
  initial <- monocle::differentialGeneTest(cds[expressed, ], fullModelFormulaStr = "~slingPseudotime_1", cores = 1)
  significant <- initial[!is.na(initial$qval) & initial$qval < 0.01, , drop = FALSE]
  significant <- significant[order(significant$qval), , drop = FALSE]
  if (nrow(significant) < 4) stop("Fewer than four pseudotime-associated genes.")
  write_csv(cbind(gene = rownames(initial), initial), file.path(d, "Monocle2_initial_gene_tests.csv"))
  ordering <- utils::head(rownames(significant), 100)
  cds <- monocle::setOrderingFilter(cds, ordering)
  cds <- monocle::reduceDimension(cds, max_components = 2, method = "DDRTree")
  cds <- monocle::orderCells(cds)
  # Save the independent DDRTree coordinate before replacing display Pseudotime.
  cds$Monocle2_Pseudotime_before_replacement <- cds$Pseudotime
  cds$Pseudotime <- cds$slingPseudotime_1
  genes_heatmap <- utils::head(rownames(significant), 1500)
  smooth <- monocle::differentialGeneTest(cds[genes_heatmap, ], fullModelFormulaStr = "~sm.ns(slingPseudotime_1)", cores = 1)
  write_csv(cbind(gene = rownames(smooth), smooth), file.path(d, "Monocle2_smooth_pseudotime_gene_tests.csv"))
  heatmap <- with_pdf_device(file.path(d, "Monocle2_Slingshot_pseudotime_four_modules.pdf"), 8, 11, {
    monocle::plot_pseudotime_heatmap(cds[rownames(smooth), ], num_clusters = 4, show_rownames = FALSE, return_heatmap = TRUE)
  })
  if (is.null(heatmap$tree_row)) stop("Monocle2 did not return a gene clustering tree; four-module analysis is incomplete.")
  if (!is.null(heatmap$tree_row)) {
    cl <- stats::cutree(heatmap$tree_row, k = 4)
    modules <- data.frame(gene = names(cl), cluster = paste0("C", cl)); write_csv(modules, file.path(d, "pseudotime_four_module_genes.csv"))
    enrich_gene_groups(split(modules$gene, modules$cluster), d)
  }
  save_plot(monocle::plot_cell_trajectory(cds, color_by = "slingPseudotime_1"), file.path(d, "Monocle2_coordinates_Slingshot_pseudotime.pdf"), 7, 6)
  # Original Monocle2 notebook cells66-69 and the Figure6G manual display set.
  selected <- intersect(c("Npm1", "Atp5a1", "Gata1", "Klf1", "Alas2", "Hba-a1", "Bnip3l", "Ulk1"), rownames(cds))
  if (length(selected)) save_plot(monocle::plot_genes_in_pseudotime(cds[selected, ], color_by = "slingPseudotime_1"), file.path(d, "erythroid_representative_gene_trends.pdf"), 10, 10)
  saveRDS(cds, file.path(d, "Monocle2_Slingshot_gene_analysis.rds"))
}
run_timepoint_trajectories <- function() {
  for (day in c("D0", "D1", "D7")) { run_slingshot(day); run_monocle_genes(day) }
}
run_nvu <- function() {
  run_lineage_panels("NVU", "Figure6_S12", "Brain")
  need(c("Seurat", "dplyr", "tidyr", "ggplot2")); obj <- get_lineage("NVU", "Brain"); d <- outdir("Figure6_S12", "NVU")
  run_nvu_all_cell_panels(obj, get_main(), d)
  markers <- c("Plvap", "Kdr", "Igfbp3", "Cldn5", "Slco1a4", "Slc2a1", "Vtn", "Pdgfrb", "Kcnj8", "Acta2", "Tagln", "Myh11", "Pdgfra", "Dcn", "Col1a1", "Gfap", "Cd44", "Serpina3n", "Aldoc", "Aqp4", "Gja1")
  markers <- intersect(markers, rownames(obj))
  save_plot(Seurat::DotPlot(obj, features = markers, group.by = "level3") + Seurat::RotatedAxis(), file.path(d, "NVU_seven_state_marker_dotplot.pdf"), 13, 6)
  expression <- Seurat::FetchData(obj, vars = c("Day", "level3", markers)) |>
    tidyr::pivot_longer(dplyr::all_of(markers), names_to = "gene", values_to = "expression") |>
    dplyr::group_by(Day, level3, gene) |> dplyr::summarise(mean_expression = mean(expression), .groups = "drop")
  write_csv(expression, file.path(d, "NVU_marker_expression_by_day.csv"))
  p <- ggplot2::ggplot(expression, ggplot2::aes(factor(Day, levels = day_order), mean_expression, color = level3, group = level3)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~gene, scales = "free_y", ncol = 4) +
    ggplot2::scale_color_manual(values = palette_for(subtype_orders$NVU)) + ggplot2::theme_bw() + ggplot2::labs(x = NULL)
  save_plot(p, file.path(d, "NVU_marker_temporal_expression.pdf"), 15, 13)
}
