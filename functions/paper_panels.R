# Additional manuscript panels recovered from original notebooks.
# Source cells and author selections are recorded in provenance/config CSVs.
paper_save_pair <- function(p, directory, name, width, height, dpi = 300L) {
  ensure_dir(directory)
  save_plot(p, file.path(directory, paste0(name, ".pdf")), width, height)
  ggplot2::ggsave(file.path(directory, paste0(name, ".png")), plot = p,
    width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
}

# Original NVU notebook cells13/15/16: within-NVU composition times the
# three-parent fraction in all cells. The two factors are exported separately.
calculate_nvu_all_cell_proportions <- function(main_meta, nvu_meta) {
  check_columns(main_meta, c("Day", "level2"), "All-cell NVU denominator")
  check_columns(nvu_meta, c("Day", "level3"), "Local NVU composition")
  if (!nrow(main_meta) || !nrow(nvu_meta) || anyNA(main_meta[, c("Day", "level2")]) || anyNA(nvu_meta[, c("Day", "level3")])) stop("Missing NVU proportion metadata.")
  if (any(!main_meta$Day %in% day_order) || any(!nvu_meta$Day %in% main_meta$Day) || any(!nvu_meta$level3 %in% subtype_orders$NVU)) stop("Invalid NVU days or subtype labels.")
  days <- intersect(day_order, unique(as.character(main_meta$Day)))
  totals <- do.call(rbind, lapply(days, function(day) {
    tab <- main_meta[main_meta$Day == day, , drop = FALSE]
    data.frame(Day = day, n_all = nrow(tab), n_three = sum(tab$level2 %in% lineage_parent$NVU))
  }))
  totals$three_type_prop <- totals$n_three / totals$n_all
  grid <- expand.grid(Day = days, level3 = subtype_orders$NVU, stringsAsFactors = FALSE)
  counts <- as.data.frame(table(factor(nvu_meta$Day, levels = days), factor(nvu_meta$level3, levels = subtype_orders$NVU)))
  names(counts) <- c("Day", "level3", "n_subtype")
  grid <- merge(grid, counts, by = c("Day", "level3"), sort = FALSE)
  denominators <- setNames(vapply(days, function(day) sum(nvu_meta$Day == day), integer(1)), days)
  grid$n_local_nvu <- unname(denominators[grid$Day])
  if (any(grid$n_local_nvu == 0L)) stop("A required day lacks local NVU subtype cells; cannot distribute the three-parent total.")
  grid$subtype_prop_inner <- grid$n_subtype / grid$n_local_nvu
  grid <- merge(grid, totals, by = "Day", sort = FALSE)
  grid$plot_height <- grid$subtype_prop_inner * grid$three_type_prop
  grid$percentage <- 100 * grid$plot_height
  list(totals = totals, composition = grid)
}
run_nvu_all_cell_panels <- function(obj, main, directory) {
  need(c("ggplot2", "scales", "Seurat"))
  r <- calculate_nvu_all_cell_proportions(main@meta.data, obj@meta.data)
  write_csv(r$totals, file.path(directory, "Fig6K_all_cell_denominator_by_Day.csv"))
  write_csv(r$composition, file.path(directory, "Fig6K_NVU_subtype_weighted_proportions.csv"))
  tab <- r$composition
  order <- c("Angiogenic-like_endo", "homeostatic_endo", "homeostatic_pericyte", "SMC", "VLMCs", "DAA", "homeostatic_Astro")
  tab$Day <- factor(tab$Day, levels = day_order); tab$level3 <- factor(tab$level3, levels = order)
  p <- ggplot2::ggplot(tab, ggplot2::aes(Day, plot_height, fill = level3)) +
    ggplot2::geom_col(width = 0.75, colour = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = palette_for(order), breaks = rev(order), drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = ggplot2::expansion(mult = c(0, .05))) +
    ggplot2::theme_classic() + ggplot2::labs(x = "Day", y = "Proportion in all cells", fill = "Subtype")
  paper_save_pair(p, directory, "Fig6K_NVU_proportion_in_all_cells", 7, 6)
  p <- ggplot2::ggplot(tab, ggplot2::aes(Day, percentage, color = level3, group = level3)) +
    ggplot2::geom_line(linewidth = .8) + ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = palette_for(order)) + ggplot2::theme_classic() +
    ggplot2::labs(x = "Day", y = "Proportion in all cells (%)", color = "Subtype")
  paper_save_pair(p, directory, "S12E_NVU_proportion_in_all_cells_lines", 10, 6)
  p <- Seurat::DimPlot(obj, group.by = "level2", cols = level2_colors[lineage_parent$NVU], raster = FALSE) + ggplot2::coord_fixed()
  paper_save_pair(p, directory, "S12D_NVU_three_level2_UMAP", 8, 7)
  run_nvu_source_expression_panels(obj, directory)
}
run_nvu_source_expression_panels <- function(obj, directory) {
  groups <- list(Astro = list(genes = c("Slc1a3","Gfap","Mt1","Mt2"), cells = c("homeostatic_Astro","DAA")),
    Endo = list(genes = c("Cldn5","Plvap"), cells = c("homeostatic_endo","Angiogenic-like_endo")),
    Pericyte = list(genes = c("Kcnj8","Acta2","Dcn","Pdgfrb"), cells = c("homeostatic_pericyte","SMC","VLMCs")))
  for (name in names(groups)) {
    spec <- groups[[name]]; genes <- intersect(spec$genes, rownames(obj))
    write_csv(data.frame(gene = spec$genes, present = spec$genes %in% genes), file.path(directory, paste0("source_NVU_",name,"_genes.csv")))
    if (!length(genes)) next
    tab <- Seurat::FetchData(obj, vars = c("Day","level3",genes)) |>
      dplyr::filter(level3 %in% spec$cells) |>
      tidyr::pivot_longer(dplyr::all_of(genes), names_to="gene", values_to="expression") |>
      dplyr::group_by(Day,level3,gene) |>
      dplyr::summarise(mean=mean(expression), n_cells=dplyr::n(), cell_SE=stats::sd(expression)/sqrt(dplyr::n()), .groups="drop")
    write_csv(tab, file.path(directory,paste0("source_NVU_",name,"_mean_cell_SE.csv")))
    p <- ggplot2::ggplot(tab,ggplot2::aes(factor(Day,levels=day_order),mean,color=level3,group=level3)) +
      ggplot2::geom_line(linewidth=.8) + ggplot2::geom_point(size=2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-cell_SE,ymax=mean+cell_SE),width=.15,na.rm=TRUE) +
      ggplot2::facet_wrap(~gene,scales="free_y") + ggplot2::scale_color_manual(values=palette_for(spec$cells)) +
      ggplot2::theme_classic() + ggplot2::labs(x="Day",y="Average expression",subtitle="Descriptive cell-level SE")
    paper_save_pair(p,directory,paste0("source_NVU_",name,"_mean_SEM_lines"),10,8)
  }
}

select_go_display_terms <- function(tab, top_n = 3L, selected_descriptions = NULL) {
  check_columns(tab, c("group", "Description", "p.adjust", "Count"), "GO display input")
  if (length(top_n) != 1L || !is.finite(top_n) || top_n < 1 || top_n != as.integer(top_n)) stop("Invalid GO display top_n.")
  tab <- tab[is.finite(tab$p.adjust) & tab$p.adjust < cfg$go_p, , drop = FALSE]
  if ("ONTOLOGY" %in% names(tab)) tab <- tab[tab$ONTOLOGY == "BP", , drop = FALSE]
  if (!is.null(selected_descriptions)) return(tab[tab$Description %in% selected_descriptions, , drop = FALSE])
  tab <- tab[order(tab$group, tab$p.adjust, -tab$Count, tab$Description), , drop = FALSE]
  tab[ave(seq_len(nrow(tab)), tab$group, FUN = seq_along) <= top_n, , drop = FALSE]
}
run_go_bubble_panel <- function(tab, directory) {
  # Edit these two values to select the GO terms to display; no historical CSV is required.
  top_n <- 3L
  selected_descriptions <- NULL
  selected <- select_go_display_terms(tab, top_n, selected_descriptions)
  write_csv(selected, file.path(directory, "GO_BP_display_selected.csv"))
  if (!nrow(selected)) return(invisible(NULL))
  p <- ggplot2::ggplot(selected, ggplot2::aes(group, Description, size = Count, color = -log10(pmax(p.adjust, .Machine$double.xmin)))) +
    ggplot2::geom_point() + ggplot2::scale_color_viridis_c() + ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = NULL, color = "-log10(adjusted p)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  paper_save_pair(p, directory, "GO_BP_selected_bubble", max(8, length(unique(selected$group))*.7), max(5, length(unique(selected$Description))*.22))
}

manual_gene_subset <- function(summary, choices) {
  check_columns(summary, c("gene", "level3", "Day", "mean"), "Manual gene expression summary")
  check_columns(choices, c("gene", "level3"), "Author gene/subtype choices")
  merge(summary, unique(choices[, c("gene", "level3"), drop = FALSE]), by = c("gene", "level3"), sort = FALSE)
}
run_manual_gene_subtype_panel <- function(summary, key, panel, directory) {
  choices <- read_csv(project_path("config", "manual_gene_subtypes.csv"))
  choices <- choices[choices$key == key & choices$panel_id == panel, , drop = FALSE]
  if (!nrow(choices)) return(invisible(NULL))
  choices$status <- ifelse(choices$level3 %in% c("CD8-ODG", paste0("Neutrophil", 10:13), "Microglia9", "Macrophage13"), "excluded_from_confirmed_final_annotation", "selected_by_author")
  active <- choices[choices$status == "selected_by_author", , drop = FALSE]
  selected <- manual_gene_subset(summary, active)
  availability <- paste(selected$gene, selected$level3, sep = "\t")
  choices$observed_in_expression_data <- paste(choices$gene, choices$level3, sep = "\t") %in% availability
  write_csv(choices, file.path(directory, paste0(panel, "_author_selection_audit.csv")))
  write_csv(selected, file.path(directory, paste0(panel, "_author_selected_expression.csv")))
  if (!nrow(selected)) return(invisible(NULL))
  selected$time_display <- unname(c(D0 = 0, D1 = 1, D7 = 3, M6 = 11)[as.character(selected$Day)])
  days <- intersect(day_order, unique(as.character(selected$Day)))
  p <- ggplot2::ggplot(selected, ggplot2::aes(time_display, mean, color = level3, group = level3)) +
    ggplot2::geom_point(size = 2) + ggplot2::geom_line(linewidth = .8) +
    ggplot2::facet_wrap(~gene, scales = "free_y", ncol = 3) + ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = unname(c(D0=0,D1=1,D7=3,M6=11)[days]), labels = days) +
    ggplot2::scale_color_manual(values = palette_for(unique(selected$level3))) +
    ggplot2::labs(x = "Time (source display spacing)", y = "Average expression", color = "Author-selected subtype")
  paper_save_pair(p, directory, paste0(panel, "_author_selected_subtype_trends"), 12, max(5, ceiling(length(unique(selected$gene))/3)*2.5))
}

run_source_lineage_panels <- function(obj, lineage, directory) {
  need(c("Seurat", "ggplot2"))
  order <- intersect(subtype_orders[[lineage]], as.character(unique(obj$level3)))
  for (tissue in intersect(c("Brain", "Skull"), unique(as.character(obj$Group)))) {
    local <- subset(obj, cells = colnames(obj)[obj$Group == tissue])
    p <- Seurat::DimPlot(local, group.by = "level3", label = FALSE, cols = palette_for(order), raster = FALSE) + ggplot2::coord_fixed()
    paper_save_pair(p, directory, paste0("subtype_UMAP_", tissue), 10, 8)
    p <- Seurat::DimPlot(local, group.by = "level3", split.by = "Day", label = FALSE, cols = palette_for(order), raster = FALSE) + ggplot2::coord_fixed()
    paper_save_pair(p, directory, paste0("subtype_UMAP_", tissue, "_by_Day"), if (tissue == "Brain") 16 else 12, 6)
  }
  genes <- read_csv(project_path("config", "paper_marker_genes.csv"))
  genes <- genes[genes$lineage == lineage, , drop = FALSE]; genes <- genes[order(genes$order), , drop = FALSE]
  if (nrow(genes)) {
    genes$present <- genes$gene %in% rownames(obj)
    write_csv(genes, file.path(directory, "source_marker_gene_availability.csv"))
    present <- unique(genes$gene[genes$present])
    if (length(present)) {
      p <- Seurat::DotPlot(obj, features = present, group.by = "level3") + Seurat::RotatedAxis()
      paper_save_pair(p, directory, "source_selected_marker_dotplot", max(8, length(present)*.35), 6)
    }
  }
  # Source subtype proportion panels were manually selected from these curves.
  # Export every subtype as a separate panel so no six-subtype selection is guessed.
  path <- file.path(directory, "proportions_pooled.csv")
  if (file.exists(path)) {
    tab <- read_csv(path); tab$Day <- factor(tab$Day, levels = day_order)
    p <- ggplot2::ggplot(tab, ggplot2::aes(Day, percentage, color = Group, group = Group)) +
      ggplot2::geom_line(linewidth = 1.1) + ggplot2::geom_point(shape = 21, fill = "white", size = 2.5) +
      ggplot2::facet_wrap(~cell_type, scales = "free_y", ncol = 3) + ggplot2::theme_bw() +
      ggplot2::scale_color_manual(values = c(Brain = "#2B83BA", Skull = "#A6316B")) + ggplot2::labs(x = NULL, y = "Pooled proportion (%)")
    paper_save_pair(p, directory, "source_subtype_temporal_Brain_vs_Skull", 12, max(5, ceiling(length(order)/3)*3))
  }
  run_library_proportion_line_panels(directory, order)
}
run_library_proportion_line_panels <- function(directory, order) {
  path <- file.path(directory, "proportion_summary_mean_sem_n.csv")
  if (!file.exists(path)) return(invisible(NULL))
  tab <- read_csv(path); tab <- tab[tab$Group == "Brain", , drop = FALSE]
  if (!nrow(tab)) return(invisible(NULL))
  tab$Day <- factor(tab$Day, levels = day_order)
  tests_path <- file.path(directory, "proportion_Student_tests_raw_and_BH.csv")
  tests <- if (file.exists(tests_path)) read_csv(tests_path) else data.frame()
  for (ct in intersect(order, unique(tab$cell_type))) {
    local <- tab[tab$cell_type == ct, , drop = FALSE]
    p <- ggplot2::ggplot(local, ggplot2::aes(Day, mean, group = 1)) +
      ggplot2::geom_line(color = palette_for(ct), linewidth = 1) + ggplot2::geom_point(color = palette_for(ct), size = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean-sem, ymax = mean+sem), width = .15, na.rm = TRUE) +
      ggplot2::theme_classic() + ggplot2::labs(title = ct, x = NULL, y = "Mean per-library proportion (%)")
    if (nrow(tests)) {
      sig <- tests[tests$Group == "Brain" & tests$cell_type == ct & tests$day1 == "D0" & is.finite(tests$p_raw), , drop = FALSE]
      if (nrow(sig)) {
        sig$star <- ifelse(sig$p_raw < .001, "***", ifelse(sig$p_raw < .01, "**", ifelse(sig$p_raw < .05, "*", "ns")))
        height <- max(c(0, local$mean + ifelse(is.na(local$sem),0,local$sem)), na.rm = TRUE)*1.1 + .1
        p <- p + ggplot2::geom_text(data = sig, ggplot2::aes(x = factor(day2, levels = day_order), y = height, label = star), inherit.aes = FALSE)
      }
    }
    name <- paste0("Brain_library_mean_line_", gsub("[^[:alnum:]_-]", "_", ct))
    paper_save_pair(p, directory, name, 5, 4)
  }
}

communication_direction_metrics <- function(tab, days, classes) {
  check_columns(tab, c("Day", "Source_Class", "Target_Class", "Strength"), "Direction metrics")
  directions <- unique(c(paste(classes[1], classes[2], sep = "->"), paste(classes[2], classes[1], sep = "->")))
  rows <- lapply(c(directions, "Total"), function(direction) {
    local <- if (direction == "Total") tab else tab[paste(tab$Source_Class, tab$Target_Class, sep = "->") == direction, , drop = FALSE]
    do.call(rbind, lapply(days, function(day) {
      values <- local$Strength[as.character(local$Day) == day]
      data.frame(Direction = direction, Day = day, Count = length(values), Weight = sum(values))
    }))
  })
  do.call(rbind, rows)
}
communication_role_frequencies <- function(tab, by_day = FALSE) {
  check_columns(tab, c("Source", "Target", "Day"), "LR role frequencies")
  if (!nrow(tab)) {
    value <- data.frame(Cell = character(), Role = character(), frequency = integer())
    if (by_day) value$Day <- character()
    return(value)
  }
  sender <- data.frame(Cell = tab$Source, Role = "Sender")
  receiver <- data.frame(Cell = tab$Target, Role = "Receiver")
  if (by_day) { sender$Day <- as.character(tab$Day); receiver$Day <- as.character(tab$Day) }
  rows <- rbind(sender, receiver)
  counts <- stats::aggregate(rep(1L, nrow(rows)), rows, sum); names(counts)[ncol(counts)] <- "frequency"; counts
}
run_source_cellchat_panels <- function(tab, tissue, category, directory) {
  need(c("ggplot2", "dplyr"))
  days <- if (tissue == "Brain") day_order else day_order[1:3]
  classes <- switch(category, myeloid_lymphoid = c("Myeloid", "Lymphoid"), myeloid_neuroglia = c("Myeloid", "Neuroglia"),
    lymphoid_neuroglia = c("Lymphoid", "Neuroglia"), intra_myeloid = rep("Myeloid",2), intra_lymphoid = rep("Lymphoid",2), intra_neuroglia = rep("Neuroglia",2))
  metrics <- communication_direction_metrics(tab, days, classes)
  write_csv(metrics, file.path(directory, "source_directional_Count_Weight.csv"))
  metrics$Day <- factor(metrics$Day, levels = days)
  for (metric in c("Count", "Weight")) {
    p <- ggplot2::ggplot(metrics, ggplot2::aes(Day, .data[[metric]], color = Direction, group = Direction)) +
      ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(size = 2) + ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = metric)
    paper_save_pair(p, directory, paste0("source_directional_", metric), 8, 5)
    p <- ggplot2::ggplot(metrics[metrics$Direction == "Total", , drop = FALSE], ggplot2::aes(Day, .data[[metric]])) +
      ggplot2::geom_col(fill = if (metric == "Count") "#E64B35" else "#4DBBD5") + ggplot2::theme_classic() + ggplot2::labs(x = NULL, y = metric)
    paper_save_pair(p, directory, paste0("source_total_", metric, "_bars"), 5, 5)
  }
  focus <- read_csv(project_path("config", "cellchat_focus_lr.csv"))
  selected <- focus$LR_pair[focus$tissue == tissue & focus$category == category]
  hits <- tolower(tab$LR_pair) %in% tolower(selected) | tolower(paste(tab$Ligand,tab$Receptor,sep="_")) %in% tolower(selected)
  for (scope in c("all_significant", "curated_LR")) {
    local <- if (scope == "curated_LR") tab[hits, , drop = FALSE] else tab
    frequencies <- communication_role_frequencies(local)
    write_csv(frequencies, file.path(directory, paste0("source_sender_receiver_", scope, ".csv")))
    by_day <- communication_role_frequencies(local, by_day = TRUE)
    write_csv(by_day, file.path(directory, paste0("source_sender_receiver_", scope, "_by_Day.csv")))
    if (nrow(by_day)) {
      by_day$Day <- factor(by_day$Day, levels = days)
      p_day <- ggplot2::ggplot(by_day, ggplot2::aes(frequency, Cell, fill = Role)) +
        ggplot2::geom_col(position = "dodge", width = .7) + ggplot2::facet_wrap(~Day, nrow = 1) +
        ggplot2::scale_fill_manual(values = c(Sender = "#E64B35", Receiver = "#3C5488")) +
        ggplot2::theme_bw() + ggplot2::labs(x = "Number of L-R interactions", y = NULL)
      paper_save_pair(p_day, directory, paste0("source_sender_receiver_", scope, "_by_Day_bars"), 16, max(6, length(unique(by_day$Cell))*.25+2))
    }
    if (!nrow(frequencies)) next
    p <- ggplot2::ggplot(frequencies, ggplot2::aes(frequency, reorder(Cell, frequency), fill = Role)) +
      ggplot2::geom_col(position = "dodge", width = .7) +
      ggplot2::scale_fill_manual(values = c(Sender = "#185FA5", Receiver = "#0F6E56")) +
      ggplot2::theme_bw() + ggplot2::labs(x = "Number of L-R interactions", y = NULL)
    paper_save_pair(p, directory, paste0("source_sender_receiver_", scope, "_bars"), 8, max(6, length(unique(frequencies$Cell))*.3+2))
  }
  run_author_cellchat_core_panels(tab[hits, , drop = FALSE], tissue, category, directory, days)
}

filter_author_cellchat_core <- function(tab, choice) {
  check_columns(choice, c("rule", "group", "level3"), "Author CellChat core cells")
  a <- choice$level3[choice$group == "A"]; b <- choice$level3[choice$group == "B"]
  rule <- unique(choice$rule)
  if (length(rule) != 1L || !rule %in% c("either", "cross")) stop("Invalid author CellChat core rule.")
  keep <- if (rule == "either") tab$Source %in% a | tab$Target %in% a else
    (tab$Source %in% a & tab$Target %in% b) | (tab$Source %in% b & tab$Target %in% a)
  tab[keep, , drop = FALSE]
}
run_author_cellchat_core_panels <- function(curated, tissue, category, directory, days) {
  choices <- read_csv(project_path("config", "paper_cellchat_core_cells.csv"))
  choices <- choices[choices$tissue == tissue & choices$category == category, , drop = FALSE]
  exclusions <- read_csv(project_path("config", "paper_cellchat_display_exclusions.csv"))
  for (panel in unique(choices$panel_id)) {
    choice <- choices[choices$panel_id == panel, , drop = FALSE]
    drop <- exclusions[exclusions$tissue == tissue & exclusions$category == category & exclusions$panel_id == panel, , drop = FALSE]
    pair_keys <- paste(curated$Source, curated$Target, sep = "|")
    all_drop <- paste(drop$Source[drop$scope == "all"], drop$Target[drop$scope == "all"], sep = "|")
    all <- curated[!pair_keys %in% all_drop, , drop = FALSE]
    main <- filter_author_cellchat_core(all, choice)
    main_drop <- paste(drop$Source[drop$scope == "main"], drop$Target[drop$scope == "main"], sep = "|")
    main <- main[!paste(main$Source, main$Target, sep = "|") %in% main_drop, , drop = FALSE]
    write_csv(choice, file.path(directory, paste0(panel, "_core_cells_used.csv")))
    write_csv(drop, file.path(directory, paste0(panel, "_display_exclusions_used.csv")))
    for (scope in c("ALL", "MAIN")) {
      local <- if (scope == "ALL") all else main
      name <- paste0(panel, "_source_LR_", scope)
      write_csv(local, file.path(directory, paste0(name, ".csv")))
      if (!nrow(local)) next
      local$Day <- factor(local$Day, levels = days)
      local$Pair <- paste(abbreviate_cell(local$Source), abbreviate_cell(local$Target), sep = " -> ")
      p <- ggplot2::ggplot(local, ggplot2::aes(LR_pair, Pair, size = Strength, color = Day, group = Day)) +
        ggplot2::geom_point(position = ggplot2::position_dodge(width = .7)) + ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + ggplot2::labs(x = NULL, y = NULL)
      paper_save_pair(p, directory, name, if (scope == "MAIN") 16 else 12, max(7, length(unique(local$Pair))*.2+2))
    }
  }
}

run_s100_cross_compartment_panel <- function() {
  need(c("Seurat", "dplyr", "tidyr", "ggplot2")); obj <- get_lineage("Macrophage")
  genes <- intersect(c("S100a8", "S100a9"), rownames(obj))
  if (!length(genes)) stop("Neither S100a8 nor S100a9 is available for Figure2K.")
  tab <- Seurat::FetchData(obj, vars = c("Group", "Day", genes)) |>
    tidyr::pivot_longer(dplyr::all_of(genes), names_to = "gene", values_to = "expression") |>
    dplyr::group_by(Group, Day, gene) |> dplyr::summarise(mean = mean(expression), .groups = "drop")
  d <- outdir("Figure2_S2", "Macrophage"); write_csv(tab, file.path(d, "Fig2K_S100a8_S100a9_Brain_vs_Skull.csv"))
  p <- ggplot2::ggplot(tab, ggplot2::aes(factor(Day, levels = day_order), mean, color = Group, group = Group)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~gene, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(Brain = "#2B83BA", Skull = "#A6316B")) + ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = "Average normalized expression")
  paper_save_pair(p, d, "Fig2K_S100a8_S100a9_Brain_vs_Skull", 10, 5)
  fc <- tab |> dplyr::group_by(Group, gene) |>
    dplyr::mutate(D0_mean = mean[as.character(Day) == "D0"][1],
      fold_change = ifelse(is.finite(D0_mean) & D0_mean > 0, mean / D0_mean, NA_real_),
      status = ifelse(is.finite(D0_mean) & D0_mean > 0, "available", "D0_zero_or_missing")) |> dplyr::ungroup()
  write_csv(fc, file.path(d, "Fig2K_S100a8_S100a9_foldchange_vs_D0.csv"))
  p <- ggplot2::ggplot(fc, ggplot2::aes(factor(Day, levels = day_order), fold_change, color = Group, group = Group)) +
    ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(na.rm = TRUE) + ggplot2::facet_wrap(~gene, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(Brain = "#d62828", Skull = "#2166ac")) + ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = "Fold change (relative to D0)")
  paper_save_pair(p, d, "Fig2K_S100a8_S100a9_foldchange_vs_D0", 7, 4.5)
}

all_module_pair_gene_overlaps <- function(brain, skull, corresponding) {
  rows <- list()
  for (b in unique(brain$cluster)) for (s in unique(skull$cluster)) {
    genes <- intersect(brain$gene[brain$cluster == b], skull$gene[skull$cluster == s])
    if (!length(genes)) next
    rows[[length(rows)+1L]] <- data.frame(brain = b, skull = s, is_corresponding = any(corresponding$brain == b & corresponding$skull == s), n_overlap_genes = length(genes), gene = genes)
  }
  if (!length(rows)) return(data.frame(brain=character(), skull=character(), is_corresponding=logical(), n_overlap_genes=integer(), gene=character()))
  do.call(rbind, rows)
}
