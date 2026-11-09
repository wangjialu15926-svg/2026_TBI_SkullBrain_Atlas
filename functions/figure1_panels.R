# Figure 1 / S1 panels adapted from merged_clustering_annotation.ipynb cells 15 and 17 (2026-11-06).
# All plotting inputs come from the audited main_final object. No old RDS is reloaded.
atlas_level1_colors <- function() {
  tab <- read_csv(project_path("config", "figure1_level1_colors.csv"))
  check_columns(tab, c("level1", "hex"), "Figure S1 level1 colors")
  if (nrow(tab) != 10L || anyNA(tab) || anyDuplicated(tab$level1) || any(!grepl("^#[0-9A-Fa-f]{6}$", tab$hex))) stop("Invalid Figure S1 level1 palette.")
  setNames(tab$hex, tab$level1)
}
prepare_atlas_metadata <- function(meta) {
  check_columns(meta, c("Source", "Group", "Day", "level2"), "Figure1 metadata")
  if (!nrow(meta) || anyNA(meta[, c("Source", "Group", "Day", "level2")])) stop("Missing atlas metadata.")
  if (any(!meta$Group %in% c("Brain", "Skull")) || any(!meta$Day %in% day_order) || any(!meta$level2 %in% level2_order)) stop("Unexpected atlas tissue/day/level2 label.")
  aliases <- c("ref" = "ref", "our" = "our", "Public" = "ref", "In-house" = "our")
  if (any(!as.character(meta$Source) %in% names(aliases))) stop("Unrecognized Source: use ref/our or Public/In-house.")
  meta$Source <- unname(aliases[as.character(meta$Source)])
  if (any(meta$Source == "ref" & meta$Group != "Brain")) stop("Published atlas source must be Brain.")
  # Ten display classes are local to these paper panels. The main object's three
  # analytical classes, consumed by CellChat, are preserved on disk.
  display <- as.character(meta$level2)
  display[display %in% c("Promyelocyte", "Macrophage", "Microglia", "Neutrophil", "Basophil", "DC")] <- "Myeloid Cell"
  display[display %in% c("B Cell", "T Cell")] <- "Lymphoid Cell"
  meta$level1 <- display
  meta
}
atlas_pooled_proportions <- function(meta, keys, label) {
  check_columns(meta, c(keys, label), "Atlas pooled proportions")
  if (!nrow(meta) || anyNA(meta[, c(keys, label), drop = FALSE])) stop("Cannot pool empty or missing atlas metadata.")
  values <- as.data.frame(lapply(meta[, c(keys, label), drop = FALSE], as.character), stringsAsFactors = FALSE)
  counts <- stats::aggregate(rep(1L, nrow(values)), by = values, FUN = sum)
  names(counts)[ncol(counts)] <- "n"
  totals <- stats::aggregate(counts$n, by = counts[, keys, drop = FALSE], FUN = sum)
  names(totals)[ncol(totals)] <- "total"
  counts <- merge(counts, totals, by = keys, all.x = TRUE, sort = FALSE)
  counts$prop <- counts$n / counts$total
  counts
}
atlas_library_level2_proportions <- function(meta, types) {
  brain <- meta[as.character(meta$Group) == "Brain", , drop = FALSE]
  if (!nrow(brain)) stop("No Brain cells for S1H.")
  # Denominator is every Brain cell in each library, before selecting the 13 panels.
  grid <- build_library_proportions(brain, "level2", level2_order)
  grid <- grid[grid$cell_type %in% types, , drop = FALSE]
  grid$level2 <- factor(grid$cell_type, levels = types); grid$prop <- grid$percentage
  grid$Day <- factor(grid$Day, levels = day_order)
  grid
}
atlas_adjacent_tests <- function(lib_prop, types, days = day_order) {
  check_columns(lib_prop, c("Lib", "Day", "level2", "prop"), "S1H library proportions")
  if (anyDuplicated(lib_prop[, c("Lib", "Day", "level2")])) stop("Duplicate library/day/cell-type observations.")
  if (length(days) < 2L || anyDuplicated(days)) stop("Invalid adjacent days.")
  rows <- list()
  for (ct in types) for (i in seq_len(length(days) - 1L)) {
    a <- lib_prop$prop[lib_prop$level2 == ct & lib_prop$Day == days[i]]
    b <- lib_prop$prop[lib_prop$level2 == ct & lib_prop$Day == days[i+1L]]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    status <- if (length(a) < 2L || length(b) < 2L) "insufficient_libraries" else "tested"
    fit <- if (status == "tested") tryCatch(stats::t.test(b, a, alternative = "two.sided", var.equal = TRUE), error = identity) else NULL
    if (inherits(fit, "error")) status <- paste0("not_testable: ", conditionMessage(fit))
    p <- if (status == "tested") fit$p.value else NA_real_
    star <- if (!is.finite(p)) "" else if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else ""
    rows[[length(rows)+1L]] <- data.frame(level2 = ct, from = days[i], to = days[i+1L], n1 = length(a), n2 = length(b), p = p, status = status, star = star)
  }
  result <- do.call(rbind, rows); result$level2 <- factor(result$level2, levels = types)
  result
}
run_paper_atlas_panels <- function(seuobj) {
  need(c("Seurat", "ggplot2", "dplyr", "tidyr", "patchwork", "scales"))
  original_meta <- seuobj@meta.data
  seuobj@meta.data <- prepare_atlas_metadata(original_meta)
  Seurat::DefaultAssay(seuobj) <- "RNA"
  if (!"umap" %in% Seurat::Reductions(seuobj)) stop("Main object has no UMAP.")
  exports <- list()
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)

  out_dir <- outdir("Figure1_S1", "Fig1_plots")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ================= 0. Colors, factor order and plotting utilities =================
  level2.cols <- level2_colors
  seuobj$level2 <- factor(as.character(seuobj$level2), levels = names(level2.cols))
  seuobj$Day    <- factor(seuobj$Day,   levels = c("D0", "D1", "D7", "M6"))
  seuobj$Group  <- factor(seuobj$Group, levels = c("Brain", "Skull"))

  # Save both PDF and PNG; raise export errors so incomplete output is not marked complete.
  save_pair <- function(plot, name, width, height, dir = out_dir) {
    ensure_dir(dir)
    pdf_file <- file.path(dir, paste0(name, ".pdf"))
    ok <- if (capabilities("cairo")) tryCatch({
      ggplot2::ggsave(pdf_file, plot = plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf, limitsize = FALSE)
      TRUE
    }, error = function(e) FALSE) else FALSE
    if (!ok) ggplot2::ggsave(pdf_file, plot = plot, width = width, height = height, units = "in", device = grDevices::pdf, limitsize = FALSE)
    ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot = plot, width = width, height = height, units = "in", dpi = 300, bg = "white", limitsize = FALSE)
    exports[[length(exports)+1L]] <<- data.frame(name = name, directory = basename(dir), width_in = width, height_in = height, png_dpi = 300L)
    message("Saved: ", name)
  }
  get_leg <- function(p) {
    g <- ggplotGrob(p)
    for (i in which(grepl("^guide-box", g$layout$name))) {
      if (!inherits(g$grobs[[i]], "zeroGrob")) return(g$grobs[[i]])
    }
    NULL
  }

  # ================= 1A. Cell counts for the schematic, printed to the console =================
  # Source values: "our" for in-house data and "ref" for published data.
  tab1A <- seuobj@meta.data %>%
    mutate(Dataset = case_when(Source == "ref" ~ "Published brain",
                               Source == "our" & Group == "Brain" ~ "In-house brain",
                               Source == "our" & Group == "Skull" ~ "In-house skull")) %>%
    count(Dataset, Day) %>%
    tidyr::pivot_wider(names_from = Day, values_from = n, values_fill = 0) %>%
    mutate(Total = rowSums(across(where(is.numeric))))
  print(as.data.frame(tab1A), row.names = FALSE)
  write_csv(tab1A, file.path(out_dir, "Fig1A_dataset_day_cell_counts.csv"))
  cat("Total cells:", ncol(seuobj), "\n")  # computed from current final object

  # ================= 1B. UMAP without cell labels =================
  emb      <- Embeddings(seuobj, reduction = "umap")
  xlim_all <- range(emb[, 1]); ylim_all <- range(emb[, 2])      # Use the same coordinate limits for all UMAP panels to aid comparison.

  plot_umap <- function(cells = NULL, title = NULL, label = FALSE) {
    DimPlot(seuobj, reduction = "umap", cells = cells, group.by = "level2", cols = level2.cols,
            label = label, label.size = 4, repel = TRUE, raster = FALSE) +
      ggtitle(title) +
      coord_fixed(xlim = xlim_all, ylim = ylim_all) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  W <- 10; H <- 8

  save_pair(plot_umap(title = "UMAP by Cell Type (level2)"), "Fig1B_UMAP_level2_all", W, H)

  cells_brain <- colnames(seuobj)[seuobj$Group == "Brain"]
  cells_skull <- colnames(seuobj)[seuobj$Group == "Skull"]
  save_pair(plot_umap(cells_brain, "Brain parenchyma"), "Fig1B_UMAP_level2_Brain", W, H)
  save_pair(plot_umap(cells_skull, "Skull"),            "Fig1B_UMAP_level2_Skull", W, H)

  for (d in levels(seuobj$Day)) {
    cells_d <- colnames(seuobj)[seuobj$Day == d]
    save_pair(plot_umap(cells_d, d), paste0("Fig1B_UMAP_level2_", d), W, H)
  }

  # ---- 1B. Cell-type count bars, descending counts with numeric labels ----
  cnt <- as.data.frame(table(seuobj$level2)); colnames(cnt) <- c("level2", "n")
  cnt <- cnt %>% arrange(desc(n)) %>% mutate(level2 = factor(level2, levels = rev(level2)))

  write_csv(cnt, file.path(out_dir, "Fig1B_level2_cell_counts.csv"))
  p_cnt <- ggplot(cnt, aes(x = level2, y = n, fill = level2)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.1, size = 3.5) +
    coord_flip() +
    scale_fill_manual(values = level2.cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "Number of cells") +
    theme_classic(base_size = 14)
  save_pair(p_cnt, "Fig1B_cell_number_barplot", 5, 6)

  # ================= 1C. Marker-gene dot plot =================
  # Use the cell-type marker genes listed in Results section 1.
  marker_list <- list(
    "Promyelocyte" = c("Elane", "Mpo"),
    "Macrophage"   = c("Fn1", "Lyz2", "Arg1"),
    "Microglia"    = c("Cx3cr1", "P2ry12", "Tmem119"),
    "Neutrophil"   = c("S100a9", "S100a8", "Ngp"),
    "Basophil"     = c("Ms4a2", "Fcer1a", "Itga2b"),
    "DC"           = c("Bst2", "Cd74"),
    "B Cell"       = c("Cd79a", "Cd19"),
    "T Cell"       = c("Cd3d", "Trbc2"),
    "Erythrocyte"  = c("Hba-a1", "Hbb-bs", "Aqp1"),
    "Astrocyte"    = c("Slc1a2", "Slc1a3"),
    "Endothelial Cell" = c("Cdh5", "Kdr", "Flt1"),
    "Pericyte"     = c("Rgs5", "Pdgfrb", "Vtn"),
    "Oligodendrocyte" = c("Mbp", "Plp1", "Mal"),
    "Neuron"       = c("Stmn2", "Ccnd2", "Sox4"),
    "CPEC"         = c("Ttr", "Otx2"),
    "Ependymal Cell" = c("Foxj1", "Tmem212")
  )
  genes_dot <- unique(unlist(marker_list))
  miss <- setdiff(genes_dot, rownames(seuobj))
  if (length(miss) > 0) stop("Genes absent from the object: ", paste(miss, collapse = ", "))

  DefaultAssay(seuobj) <- "RNA"        # Change this if the expression matrix uses a different assay.
  p_dot <- DotPlot(seuobj, features = genes_dot, group.by = "level2", dot.scale = 6,
                   col.min = -1, col.max = 2) +
    scale_color_viridis_c(name = "Average Exp.", limits = c(-1, 2), oob = scales::squish) +
    scale_radius(range = c(0, 6), name = "Percent Exp.") +
    scale_y_discrete(limits = rev(names(level2.cols))) +          # Display cell types from top to bottom in level2 order.
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_line(colour = "grey92"),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
  save_pair(p_dot, "Fig1C_marker_dotplot", 14, 6.5)

  # ================= 1D. Stacked cell-type proportions by tissue and time point =================
  bar_df <- atlas_pooled_proportions(seuobj@meta.data, c("Group", "Day"), "level2")
  bar_df$Group <- factor(bar_df$Group, levels = c("Brain", "Skull"))
  bar_df$Day <- factor(bar_df$Day, levels = day_order)
  bar_df$level2 <- factor(bar_df$level2, levels = names(level2.cols))
  write_csv(bar_df, file.path(out_dir, "Fig1D_pooled_proportions.csv"))

  p_D <- ggplot(bar_df, aes(x = Day, y = prop, fill = level2)) +
    geom_col(width = 0.9) +
    facet_grid(. ~ Group, scales = "free_x") +             # Brain and Skull facets have equal widths; individual Skull bars are wider.
    scale_fill_manual(values = level2.cols, breaks = names(level2.cols), name = "Cell type") +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(x) sprintf("%.2f", x),
                       expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Time", y = "Cell proportion") +
    theme_bw(base_size = 14) +
    theme(panel.grid       = element_blank(),
          panel.border     = element_rect(colour = "grey60", fill = NA, linewidth = 0.8),
          strip.background = element_rect(fill = "grey92", colour = "grey60"),
          strip.text       = element_text(size = 14),
          axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position  = "right",
          legend.title     = element_text(size = 14),
          legend.text      = element_text(size = 12))

  save_pair(p_D + theme(legend.position = "none"), "Fig1D_panel", 6.0, 6.1)          # No legend; retain the y axis.
  save_pair(p_D + theme(legend.position = "none",
                        axis.text.y = element_blank(), axis.title.y = element_blank()),
            "Fig1D_panel_noYlabel", 5.15, 6.1)                                       # Match the supplied reference panel dimensions.
  leg <- get_leg(p_D)
  if (is.null(leg)) stop("Cannot extract required legend: Fig1D")
  save_pair(patchwork::wrap_elements(full = leg), "Fig1D_legend", 2.6, 5.5)
  save_pair(p_D, "Fig1D_full", 8.6, 6.1)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)

  out_dir <- outdir("Figure1_S1", "FigS1_plots")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ================= 0. Add level1 if it is not already present =================
  level1_levels <- c("Myeloid Cell", "Lymphoid Cell", "Erythrocyte", "Astrocyte", "Endothelial Cell",
                     "Pericyte", "Oligodendrocyte", "Neuron", "CPEC", "Ependymal Cell")
  if (!"level1" %in% colnames(seuobj@meta.data)) {
    level1_map <- c(
      "Promyelocyte" = "Myeloid Cell", "Macrophage" = "Myeloid Cell", "Microglia" = "Myeloid Cell",
      "Neutrophil"   = "Myeloid Cell", "Basophil"   = "Myeloid Cell", "DC"        = "Myeloid Cell",
      "B Cell" = "Lymphoid Cell", "T Cell" = "Lymphoid Cell",
      "Erythrocyte" = "Erythrocyte", "Astrocyte" = "Astrocyte", "Endothelial Cell" = "Endothelial Cell",
      "Pericyte" = "Pericyte", "Oligodendrocyte" = "Oligodendrocyte", "Neuron" = "Neuron",
      "CPEC" = "CPEC", "Ependymal Cell" = "Ependymal Cell"
    )
    seuobj$level1 <- unname(level1_map[as.character(seuobj$level2)])
  }

  # ================= 1. Colors, factor order and plotting utilities =================
  level1.cols <- atlas_level1_colors()
  level2.cols <- level2_colors
  chk <- function(x, cols, nm) {
    bad <- setdiff(unique(as.character(x)), names(cols))
    if (length(bad) > 0) stop(nm, " has labels without defined colors: ", paste(bad, collapse = ", "))
  }
  chk(seuobj$level1, level1.cols, "level1"); chk(seuobj$level2, level2.cols, "level2")

  seuobj$level1 <- factor(as.character(seuobj$level1), levels = level1_levels)
  seuobj$level2 <- factor(as.character(seuobj$level2), levels = names(level2.cols))
  seuobj$Day    <- factor(seuobj$Day,   levels = c("D0", "D1", "D7", "M6"))
  seuobj$Group  <- factor(seuobj$Group, levels = c("Brain", "Skull"))

  # Save PDF using cairo_pdf when available; propagate export failures.
  save_pair <- function(plot, name, width, height, dir = out_dir) {
    ensure_dir(dir)
    pdf_file <- file.path(dir, paste0(name, ".pdf"))
    ok <- if (capabilities("cairo")) tryCatch({
      ggplot2::ggsave(pdf_file, plot = plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf, limitsize = FALSE)
      TRUE
    }, error = function(e) FALSE) else FALSE
    if (!ok) ggplot2::ggsave(pdf_file, plot = plot, width = width, height = height, units = "in", device = grDevices::pdf, limitsize = FALSE)
    ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot = plot, width = width, height = height, units = "in", dpi = 300, bg = "white", limitsize = FALSE)
    exports[[length(exports)+1L]] <<- data.frame(name = name, directory = basename(dir), width_in = width, height_in = height, png_dpi = 300L)
    message("Saved: ", name)
  }
  get_leg <- function(p) {
    g <- ggplotGrob(p)
    for (i in which(grepl("^guide-box", g$layout$name))) {
      if (!inherits(g$grobs[[i]], "zeroGrob")) return(g$grobs[[i]])
    }
    NULL
  }
  # Save stacked plots as _panel (no legend), _legend (legend only), and _full (panel plus legend).
  save_stack <- function(p, name, panel_w, panel_h, leg_w, leg_h) {
    save_pair(p + theme(legend.position = "none"), paste0(name, "_panel"), panel_w, panel_h)
    leg <- get_leg(p)
    if (!is.null(leg)) {
      save_pair(patchwork::wrap_elements(full = leg), paste0(name, "_legend"), leg_w, leg_h)
      save_pair(p, paste0(name, "_full"), panel_w, panel_h + leg_h)
    } else stop("Cannot extract required legend: ", name)
  }

  emb      <- Embeddings(seuobj, reduction = "umap")
  xlim_all <- range(emb[, 1]); ylim_all <- range(emb[, 2])      # Use identical coordinate limits for all UMAP panels.

  plot_umap <- function(group.by, cols, cells = NULL, title = NULL, label = FALSE) {
    DimPlot(seuobj, reduction = "umap", cells = cells, group.by = group.by, cols = cols,
            label = label, label.size = 4, repel = TRUE, raster = FALSE) +
      ggtitle(title) +
      coord_fixed(xlim = xlim_all, ylim = ylim_all) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  W <- 10; H <- 8

  # ================= S1A. UMAP at level1 =================
  save_pair(plot_umap("level1", level1.cols, title = NULL) + labs(colour = "Cell class"),
            "S1A_UMAP_level1", W, H)

  # ================= S1B. UMAP by data source, colored by level2 =================
  # Source values: "our" for in-house data and "ref" for published data.
  cells_ref       <- colnames(seuobj)[seuobj$Source == "ref"]
  cells_our_brain <- colnames(seuobj)[seuobj$Source == "our" & seuobj$Group == "Brain"]
  cells_our_skull <- colnames(seuobj)[seuobj$Source == "our" & seuobj$Group == "Skull"]

  save_pair(plot_umap("level2", level2.cols, cells_ref,       "Published brain dataset"), "S1B_UMAP_ref_brain", W, H)
  save_pair(plot_umap("level2", level2.cols, cells_our_brain, "In-house brain"),          "S1B_UMAP_our_brain", W, H)
  save_pair(plot_umap("level2", level2.cols, cells_our_skull, "In-house skull"),          "S1B_UMAP_our_skull", W, H)

  # ================= S1E. FeaturePlot for level2 marker genes =================
  marker_df <- data.frame(
    celltype = c("Promyelocyte", "Macrophage", "Microglia", "Neutrophil", "Basophil", "DC",
                 "B Cell", "T Cell", "Erythrocyte", "Astrocyte", "Endothelial Cell", "Pericyte",
                 "Oligodendrocyte", "Neuron", "CPEC", "Ependymal Cell"),
    gene     = c("Mpo", "Lyz2", "Cx3cr1", "S100a9", "Itga2b", "Cd74",
                 "Cd79a", "Cd3d", "Hba-a1", "Slc1a2", "Flt1", "Pdgfrb",
                 "Mbp", "Stmn2", "Ttr", "Tmem212"),
    stringsAsFactors = FALSE)
  miss <- setdiff(marker_df$gene, rownames(seuobj))
  if (length(miss) > 0) stop("Genes absent from the object: ", paste(miss, collapse = ", "))

  DefaultAssay(seuobj) <- "RNA"        # Change this if the expression matrix uses a different assay.
  umap_df <- as.data.frame(emb); colnames(umap_df) <- c("UMAP_1", "UMAP_2")

  # Use rasterized points when scattermore is available to reduce PDF size.
  point_layer <- if (requireNamespace("scattermore", quietly = TRUE)) {
    scattermore::geom_scattermore(pointsize = 2, pixels = c(1500, 1500))
  } else {
    geom_point(size = 0.05, shape = 16)
  }

  plot_feature <- function(ct, gene) {
    df <- umap_df
    expr <- FetchData(seuobj, vars = gene)
    if (!all(rownames(df) %in% rownames(expr))) stop("Feature expression cell IDs do not match UMAP.")
    df$expr <- expr[rownames(df), gene]
    df       <- df[order(df$expr), ]                         # Draw high-expression cells last so they appear on top.
    df$label <- sprintf('"%s"~"-"~italic("%s")', ct, gene)   # Title: cell type - gene, with the gene in italics.
    ggplot(df, aes(UMAP_1, UMAP_2, colour = expr)) +
      point_layer +
      scale_colour_viridis_c(option = "plasma") +
      facet_wrap(~label, labeller = label_parsed) +
      coord_fixed(xlim = xlim_all, ylim = ylim_all) +
      theme_void() +
      theme(strip.background = element_rect(fill = "grey92", colour = "grey60"),
            strip.text   = element_text(size = 16, margin = margin(6, 0, 6, 0)),
            legend.position = "none")
  }

  feat_plots <- Map(plot_feature, marker_df$celltype, marker_df$gene)
  for (i in seq_along(feat_plots)) {
    save_pair(feat_plots[[i]],
              paste0("S1E_feature_", sprintf("%02d", i), "_", gsub(" ", "", marker_df$celltype[i]), "_", marker_df$gene[i]),
              5, 5)
  }
  save_pair(wrap_plots(feat_plots, ncol = 4), "S1E_feature_all_4x4", 20, 20)

  # ================= Shared theme for stacked plots =================
  theme_stack <- theme_bw(base_size = 14) +
    theme(panel.grid       = element_blank(),
          panel.border     = element_rect(colour = "grey60", fill = NA, linewidth = 0.8),
          strip.background = element_rect(fill = "grey92", colour = "grey60"),
          strip.text       = element_text(size = 14),
          axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position  = "bottom")

  # ================= S1F. Level1 proportions in Brain and Skull over time =================
  f_df <- atlas_pooled_proportions(seuobj@meta.data, c("Group", "Day"), "level1")
  f_df$Group <- factor(f_df$Group, levels = c("Brain", "Skull"))
  f_df$Day <- factor(f_df$Day, levels = day_order)
  f_df$level1 <- factor(f_df$level1, levels = level1_levels)
  write_csv(f_df, file.path(out_dir, "S1F_pooled_proportions.csv"))

  p_F <- ggplot(f_df, aes(x = Day, y = prop, fill = level1)) +
    geom_col(width = 0.9) +
    facet_grid(. ~ Group, scales = "free_x") +              # Use equal-width Brain and Skull facets, matching the source figure.
    scale_fill_manual(values = level1.cols, name = "Cell type") +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(x) sprintf("%.2f", x),
                       expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = NULL, y = "Cell proportion") + theme_stack +
    guides(fill = guide_legend(nrow = 4))
  # Use the source panel aspect ratio (about 4.3 x 4.9) plus space for y-axis labels.
  save_stack(p_F, "S1F_level1_proportion", panel_w = 5.2, panel_h = 4.9, leg_w = 6, leg_h = 1.6)

  # ================= S1G. Level2 in in-house versus published Brain data at D0 and D1 =================
  g_df <- seuobj@meta.data %>%
    filter(Group == "Brain", Day %in% c("D0", "D1")) %>%
    mutate(Dataset = factor(ifelse(Source == "our", "This study", "Jha Ruchira M et al., 2024"),
                            levels = c("This study", "Jha Ruchira M et al., 2024")),
           Group = "Brain", Day = droplevels(Day)) %>%
    count(Group, Dataset, Day, level2, name = "n") %>%
    group_by(Dataset, Day) %>% mutate(prop = n / sum(n)) %>% ungroup()

  write_csv(g_df, file.path(out_dir, "S1G_pooled_proportions.csv"))
  p_G <- ggplot(g_df, aes(x = Day, y = prop, fill = level2)) +
    geom_col(width = 0.9) +
    facet_wrap(vars(Group, Dataset), nrow = 1) +            # Two-line facet labels: Brain / data source.
    scale_fill_manual(values = level2.cols, name = "Cell class") +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = function(x) sprintf("%.2f", x),
                       expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = NULL, y = "Cell proportion") + theme_stack +
    theme(strip.text = element_text(size = 10)) +
    guides(fill = guide_legend(nrow = 6))
  # Use the source panel aspect ratio (about 3.6 x 5.4) plus space for y-axis labels.
  save_stack(p_G, "S1G_level2_our_vs_ref_brain", panel_w = 4.5, panel_h = 5.4, leg_w = 6, leg_h = 2.2)

  # ================= S1H. Brain level2 proportions over time, with adjacent-time tests =================
  ct_H   <- c("Macrophage", "Microglia", "Neutrophil", "DC", "B Cell", "T Cell", "Astrocyte",
              "Endothelial Cell", "Pericyte", "Oligodendrocyte", "Neuron", "CPEC", "Ependymal Cell")
  day_lv <- c("D0", "D1", "D7", "M6")

  lib_prop <- atlas_library_level2_proportions(seuobj@meta.data, ct_H)
  write_csv(lib_prop, file.path(out_dir, "S1H_proportions_per_library.csv"))

  sum_df <- lib_prop %>%
    group_by(level2, Day) %>%
    summarise(mean = mean(prop), sem = sd(prop) / sqrt(n()), .groups = "drop") %>%
    mutate(ci = match(as.character(level2), ct_H),
           di = match(as.character(Day), day_lv),
           x  = ci + (di - 2.5) * 0.2)

  # Adjacent comparisons: D0 vs D1, D1 vs D7, D7 vs M6; two-sided equal-variance Student tests on library proportions, unadjusted p.
  pv <- atlas_adjacent_tests(lib_prop, ct_H, day_lv)
  print(as.data.frame(pv), row.names = FALSE)
  write_csv(pv, file.path(out_dir, "S1H_adjacent_Student_raw_p.csv"))

  top_df   <- sum_df %>% transmute(level2, Day, top = mean + ifelse(is.na(sem), 0, sem))
  step     <- 0.045 * max(c(0, top_df$top), na.rm = TRUE)
  pair_idx <- c("D0 D1" = 1, "D1 D7" = 2, "D7 M6" = 3)

  br <- pv %>%
    filter(!is.na(p), p < 0.05) %>%
    mutate(idx = unname(pair_idx[paste(from, to)])) %>%
    left_join(top_df, by = c("level2", "from" = "Day")) %>% rename(top_from = top) %>%
    left_join(top_df, by = c("level2", "to"   = "Day")) %>% rename(top_to   = top) %>%
    group_by(level2) %>%
    # Place brackets 1 and 3 lower because their endpoints differ; raise bracket 2 (D1-D7) when a neighboring bracket is significant.
    mutate(lev = ifelse(idx == 2 & any(idx %in% c(1, 3)), 2, 1)) %>%
    ungroup() %>%
    mutate(y  = pmax(top_from, top_to) + step * ifelse(lev == 2, 2.2, 1),
           ci = match(as.character(level2), ct_H),
           x0 = ci + (match(from, day_lv) - 2.5) * 0.2 + 0.02,     # Inset both endpoints by 0.02.
           x1 = ci + (match(to,   day_lv) - 2.5) * 0.2 - 0.02)

  p_H <- ggplot(sum_df, aes(x = x, y = mean, fill = level2)) +
    geom_col(width = 0.18) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.07, linewidth = 0.4) +
    geom_text(aes(y = 0, label = Day), vjust = 1.4, size = 2.5) +
    geom_segment(data = br, aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
    geom_text(data = br, aes(x = (x0 + x1) / 2, y = y, label = star), inherit.aes = FALSE,
              vjust = -0.1, size = 4) +
    scale_fill_manual(values = level2.cols, guide = "none") +
    scale_x_continuous(breaks = seq_along(ct_H), labels = ct_H, expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    labs(x = NULL, y = "Proportion (%)") +
    theme_classic(base_size = 14) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
          axis.text.x  = element_text(angle = 45, hjust = 1))
  save_pair(p_H, "S1H_level2_proportion_adjacent_test", 16, 5.5)
  export_manifest <- do.call(rbind, exports)
  write_csv(export_manifest, file.path(outdir("Figure1_S1"), "Figure1_S1_export_manifest.csv"))
  invisible(export_manifest)
}
