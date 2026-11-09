cellchat_path <- function() project_path("results", "cellchat", "cellchat_all_full_min50.rds")
run_cellchat_build <- function() {
  need(c("CellChat", "Seurat", "future"))
  previous_plan <- future::plan(); on.exit(future::plan(previous_plan), add = TRUE)
  future::plan(future::sequential)
  obj <- get_main(); result <- list(); inventory <- list()
  for (tissue in c("Brain", "Skull")) for (day in if (tissue == "Brain") day_order else day_order[1:3]) {
    cells <- colnames(obj)[obj$Group == tissue & obj$Day == day]
    meta <- obj@meta.data[cells, , drop = FALSE]
    counts <- table(meta$level3); valid <- names(counts)[counts >= cfg$cellchat_min_cells]
    kept <- cells[meta$level3 %in% valid]
    if (length(valid) < 2) stop("Insufficient CellChat groups: ", tissue, " ", day)
    meta <- obj@meta.data[kept, , drop = FALSE]; meta$level3 <- factor(meta$level3, levels = sort(valid))
    data <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")[, kept, drop = FALSE]
    cc <- CellChat::createCellChat(data, meta = meta, group.by = "level3")
    cc@DB <- CellChat::CellChatDB.mouse
    cc <- CellChat::subsetData(cc); cc <- CellChat::identifyOverExpressedGenes(cc)
    cc <- CellChat::identifyOverExpressedInteractions(cc)
    cc <- CellChat::computeCommunProb(cc, type = "triMean", trim = 0.1, population.size = TRUE,
                                      nboot = cfg$cellchat_nboot, seed.use = cfg$cellchat_seed)
    cc <- CellChat::filterCommunication(cc, min.cells = cfg$cellchat_min_cells)
    cc <- CellChat::computeCommunProbPathway(cc, thresh = cfg$cellchat_p)
    cc <- CellChat::aggregateNet(cc, thresh = cfg$cellchat_p)
    name <- paste(tissue, day, sep = "_"); result[[name]] <- cc
    saveRDS(cc, file.path(outdir("cellchat", "objects"), paste0(name, ".rds")))
    inventory[[name]] <- data.frame(Object = name, cells_before_min50 = length(cells), cells_after_min50 = length(kept),
                                    subtypes = length(valid), downsampled = FALSE,
                                    excluded_subtypes = paste(names(counts)[counts < cfg$cellchat_min_cells], collapse = ";"))
  }
  attr(result, "main_md5") <- hash_input(main_final_path()); attr(result, "downsampled") <- FALSE
  attr(result, "settings") <- cfg[c("cellchat_min_cells", "cellchat_nboot", "cellchat_p", "cellchat_seed")]
  saveRDS(result, cellchat_path())
  write_csv(do.call(rbind, inventory), file.path(outdir("cellchat"), "full_cellchat_input_inventory.csv"))
}
get_cellchat <- function() {
  assert_no_incomplete_stages(c("3", "4"))
  need("CellChat")
  cc <- read_input(cellchat_path())
  if (!identical(attr(cc, "downsampled"), FALSE) || !identical(attr(cc, "main_md5"), hash_input(main_final_path()))) stop("CellChat list is sampled or stale; rerun script 4.")
  expected <- cfg[c("cellchat_min_cells", "cellchat_nboot", "cellchat_p", "cellchat_seed")]
  if (!identical(attr(cc, "settings"), expected)) stop("CellChat settings changed; rerun script 4.")
  required <- c(paste0("Brain_", day_order), paste0("Skull_", day_order[1:3]))
  if (anyDuplicated(names(cc)) || !setequal(names(cc), required)) stop("The full CellChat list must contain exactly the seven study objects.")
  cc
}
cellchat_category <- function(a, b) {
  output <- rep(NA_character_, length(a))
  output[a == "Myeloid" & b == "Myeloid"] <- "intra_myeloid"
  output[a == "Lymphoid" & b == "Lymphoid"] <- "intra_lymphoid"
  output[a == "Neuroglia" & b == "Neuroglia"] <- "intra_neuroglia"
  output[(a == "Myeloid" & b == "Lymphoid") | (a == "Lymphoid" & b == "Myeloid")] <- "myeloid_lymphoid"
  output[(a == "Myeloid" & b == "Neuroglia") | (a == "Neuroglia" & b == "Myeloid")] <- "myeloid_neuroglia"
  output[(a == "Lymphoid" & b == "Neuroglia") | (a == "Neuroglia" & b == "Lymphoid")] <- "lymphoid_neuroglia"
  output
}
extract_sig_communications <- function(cc, tissue, day, p_threshold = cfg$cellchat_p) {
  # Apply one explicit p mask to every plot, table and interaction count.
  probabilities <- cc@net$prob; pvalues <- cc@net$pval
  if (length(dim(probabilities)) != 3L || !identical(dim(probabilities), dim(pvalues))) stop("CellChat probability/p-value arrays must have matching 3D dimensions.")
  if (!is.numeric(p_threshold) || length(p_threshold) != 1L || !is.finite(p_threshold) || p_threshold <= 0 || p_threshold > 1) stop("Invalid communication p threshold.")
  names_lr <- dimnames(probabilities)[[3]]
  if (any(vapply(dimnames(probabilities)[1:2], function(x) is.null(x) || anyNA(x) || anyDuplicated(x) || any(!nzchar(x)), logical(1)))) stop("CellChat source/target names are missing or duplicated.")
  if (!is.null(dimnames(pvalues)) && !identical(dimnames(probabilities), dimnames(pvalues))) stop("CellChat probability/p-value dimension names differ.")
  if (anyDuplicated(names_lr)) stop("CellChat interaction_name IDs must be unique; do not collapse ligand/receptor display names.")
  lr <- cc@LR$LRsig
  if (is.null(names_lr)) stop("CellChat probabilities lack interaction IDs.")
  lr_ids <- rownames(lr)
  if (is.null(lr_ids) || !all(names_lr %in% lr_ids)) lr_ids <- if ("interaction_name" %in% names(lr)) as.character(lr$interaction_name) else character()
  if (anyDuplicated(lr_ids) || !all(names_lr %in% lr_ids)) stop("Missing or duplicated CellChat interaction identifiers.")
  chunks <- list()
  for (i in seq_along(names_lr)) {
    prob <- matrix(probabilities[, , i], nrow = dim(probabilities)[1], ncol = dim(probabilities)[2], dimnames = dimnames(probabilities)[1:2])
    pv <- matrix(pvalues[, , i], nrow = nrow(prob), ncol = ncol(prob))
    index <- which(is.finite(prob) & prob > 0 & is.finite(pv) & pv < p_threshold, arr.ind = TRUE)
    if (!nrow(index)) next
    row <- match(names_lr[i], lr_ids); if (is.na(row)) stop("Unmatched CellChat interaction ID: ", names_lr[i])
    chunks[[length(chunks) + 1L]] <- data.frame(Source = rownames(prob)[index[, 1]], Target = colnames(prob)[index[, 2]],
      LR_pair = names_lr[i], Ligand = as.character(lr$ligand[row]), Receptor = as.character(lr$receptor[row]),
      Pathway = as.character(lr$pathway_name[row]), Strength = prob[index], P_value = pv[index], Tissue = tissue, Day = day)
  }
  if (!length(chunks)) return(data.frame(Source = character(), Target = character(), LR_pair = character(), Ligand = character(), Receptor = character(), Pathway = character(), Strength = numeric(), P_value = numeric(), Tissue = character(), Day = character()))
  do.call(rbind, chunks)
}
get_communication_table <- function() {
  need("dplyr")
  # Validate the actual object settings before accepting any derived cache.
  lists <- get_cellchat()
  cache_path <- file.path(outdir("cellchat"), "significant_communications_nine_categories.rds")
  stamp <- list(main_md5 = hash_input(main_final_path()), cellchat_md5 = hash_input(cellchat_path()), p = cfg$cellchat_p,
    classifier_code_md5 = hash_input(project_path("functions", "cellchat.R")), label_code_md5 = hash_input(project_path("functions", "common.R")))
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(attr(cached, "input_stamp"), stamp)) return(cached)
  }
  main <- get_main()
  labels <- unique(main@meta.data[, c("level3", "level1")])
  if (anyDuplicated(labels$level3)) stop("A level3 label has multiple level1 classes.")
  groups <- setNames(as.character(labels$level1), labels$level3)
  groups[unique(c(subtype_orders$NVU, "Oligodendrocyte", "Neuron", "CPEC", "Ependymal Cell"))] <- "Neuroglia"
  chunks <- list(); network_checks <- list()
  for (nm in names(lists)) {
    tissue <- sub("_.*$", "", nm); day <- sub("^[^_]+_", "", nm)
    tab <- extract_sig_communications(lists[[nm]], tissue, day)
    network_checks[[nm]] <- data.frame(Object = nm, Strict_p_interaction_count = nrow(tab),
      Strict_p_weight = if (nrow(tab)) sum(tab$Strength) else 0,
      AggregateNet_count = sum(lists[[nm]]@net$count), AggregateNet_weight = sum(lists[[nm]]@net$weight),
      p_exactly_on_threshold = sum(lists[[nm]]@net$pval == cfg$cellchat_p, na.rm = TRUE))
    if (!nrow(tab)) next
    tab$Source_Class <- unname(groups[tab$Source]); tab$Target_Class <- unname(groups[tab$Target])
    tab$Source_Class[is.na(tab$Source_Class)] <- "Unclassified"; tab$Target_Class[is.na(tab$Target_Class)] <- "Unclassified"
    tab$Category <- cellchat_category(tab$Source_Class, tab$Target_Class)
    # Skull categories are deliberately limited to immune interactions, not a claim of no non-immune cells.
    if (tissue == "Skull") tab$Category[!tab$Category %in% c("intra_myeloid", "myeloid_lymphoid", "intra_lymphoid")] <- NA_character_
    chunks[[nm]] <- tab
  }
  result <- dplyr::bind_rows(chunks)
  write_csv(dplyr::bind_rows(network_checks), file.path(outdir("cellchat"), "strict_p_vs_aggregateNet_diagnostics.csv"))
  if (!nrow(result)) stop("No significant interactions under the configured threshold.")
  write_csv(result[is.na(result$Category), ], file.path(outdir("cellchat"), "interactions_outside_nine_categories.csv"))
  result <- result[!is.na(result$Category), ]
  attr(result, "input_stamp") <- stamp
  saveRDS(result, cache_path)
  result
}
score_time_table <- function(long, keys, days) {
  if (!is.character(days) || !identical(days[1], "D0") || length(days) < 2L || anyNA(days) || anyDuplicated(days)) stop("Provide unique chronological days starting with D0.")
  if (!is.character(keys) || !length(keys) || anyDuplicated(keys)) stop("Provide unique grouping keys.")
  if (!nrow(long)) {
    empty <- as.data.frame(setNames(rep(list(character()), length(keys)), keys))
    for (day in days) empty[[day]] <- numeric()
    empty$Baseline <- logical(); empty$Change_Score <- numeric(); empty$Induced_Score <- numeric(); empty$Peak <- character()
    for (day in setdiff(days, "D0")) { empty[[paste0(day, "_vs_D0_log2FC")]] <- numeric(); empty[[paste0("Lost_", day)]] <- logical() }
    empty$Lost_all_later <- logical(); return(empty)
  }
  check_columns(long, c(keys, "Day", "Strength"), "Communication scoring input")
  long <- as.data.frame(long)
  if (anyNA(long[, c(keys, "Day"), drop = FALSE]) || any(!long$Day %in% days)) stop("Communication keys/days are missing or outside the configured time course.")
  if (!is.numeric(long$Strength) || any(!is.finite(long$Strength)) || any(long$Strength < 0)) stop("Communication strengths must be finite and nonnegative.")
  tab <- unique(long[, keys, drop = FALSE])
  for (day in days) {
    slice <- long[long$Day == day, c(keys, "Strength"), drop = FALSE]
    if (nrow(slice)) {
      total <- stats::aggregate(slice$Strength, by = slice[, keys, drop = FALSE], FUN = sum)
      names(total)[ncol(total)] <- day
      tab <- merge(tab, total, by = keys, all.x = TRUE, sort = FALSE)
      tab[[day]][is.na(tab[[day]])] <- 0
    } else tab[[day]] <- 0
  }
  later <- setdiff(days, "D0"); values <- as.matrix(tab[, days, drop = FALSE]); after <- as.matrix(tab[, later, drop = FALSE])
  tab$Baseline <- tab$D0 > 0
  ratios <- log2(sweep(after + cfg$cellchat_epsilon, 1, tab$D0 + cfg$cellchat_epsilon, "/"))
  for (day in later) tab[[paste0(day, "_vs_D0_log2FC")]] <- ratios[, day]
  tab$Change_Score <- ifelse(tab$Baseline, rowSums(abs(ratios)), 0)
  tab$Induced_Score <- ifelse(!tab$Baseline, apply(after, 1, max), 0)
  tab$Peak <- days[max.col(values, ties.method = "first")]
  for (day in later) tab[[paste0("Lost_", day)]] <- tab$Baseline & tab[[day]] == 0
  tab$Lost_all_later <- tab$Baseline & rowSums(after) == 0
  tab
}
table_category_spec <- data.frame(
  table = paste0("S", 3:11), tissue = c(rep("Brain", 6), rep("Skull", 3)),
  category = c("intra_myeloid", "myeloid_neuroglia", "myeloid_lymphoid", "intra_neuroglia", "lymphoid_neuroglia", "intra_lymphoid",
                 "intra_myeloid", "myeloid_lymphoid", "intra_lymphoid"),
  plot = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE))
abbreviate_cell <- function(x) {
  x <- sub("^Macrophage", "Mφ", x); x <- sub("^Neutrophil", "Neu", x); x <- sub("^Microglia", "MG", x)
  map <- c("Monocyte" = "Mono", "CD4 T Cell" = "CD4T", "CD8 T Cell" = "CD8T", "NK Cell" = "NK", "γδ T Cell" = "γδT",
           "B Cell" = "B", "B1 Cell" = "B1", "Plasma Cell" = "Plasma", "T progenitors" = "Tprog")
  hit <- x %in% names(map); x[hit] <- map[x[hit]]; x
}
cellchat_heatmap <- function(tab, key, days, path, type) {
  need("pheatmap"); if (!nrow(tab)) return(invisible(NULL))
  if (type == "baseline") {
    tab <- tab[tab$Baseline, , drop = FALSE]; tab <- tab[order(-tab$Change_Score), , drop = FALSE]
  } else {
    tab <- tab[!tab$Baseline, , drop = FALSE]; tab <- tab[order(-tab$Induced_Score), , drop = FALSE]
  }
  if (!nrow(tab)) return(invisible(NULL))
  tab <- utils::head(tab, if (type == "baseline") 50L else 30L)
  later <- setdiff(days, "D0")
  if (type == "baseline") {
    mat <- as.matrix(tab[, paste0(later, "_vs_D0_log2FC"), drop = FALSE]); colnames(mat) <- later
    loss <- as.matrix(tab[, paste0("Lost_", later), drop = FALSE]); labels <- matrix("", nrow(mat), ncol(mat)); labels[loss] <- "•"
    finite_values <- abs(mat[!loss]); limit <- if (length(finite_values)) max(1, finite_values) else 1
    mat[mat < -limit] <- -limit; mat[mat > limit] <- limit; mat[loss] <- -limit * 1.02
    breaks <- c(-limit * 1.03, -limit * 1.001, seq(-limit, limit, length.out = 100))
    colors <- c("#082567", grDevices::colorRampPalette(c("#313695", "white", "#A50026"))(100))
  } else {
    mat <- as.matrix(tab[, later, drop = FALSE]); labels <- FALSE; breaks <- NA
    colors <- grDevices::colorRampPalette(c("white", "#FB6A4A", "#99000D"))(100)
  }
  rownames(mat) <- make.unique(as.character(tab[[key]]))
  annotation <- if ("Function_Group" %in% names(tab)) data.frame(Function_Group = tab$Function_Group, row.names = rownames(mat)) else NA
  pheatmap::pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, color = colors, breaks = breaks,
    annotation_row = annotation, display_numbers = labels, number_color = "white", cellwidth = 35, cellheight = 10, border_color = "grey90",
    filename = path, width = max(7, ncol(mat) * 35 / 72 + 5), height = max(3, nrow(mat) * 10 / 72 + 2))
}
run_cellchat_category <- function(tissue, category, make_plots = TRUE) {
  need(c("dplyr", "tidyr", "openxlsx", "ggplot2"))
  all <- get_communication_table(); tab <- all[all$Tissue == tissue & all$Category == category, ]
  spec <- table_category_spec[table_category_spec$tissue == tissue & table_category_spec$category == category, ]
  if (nrow(spec) != 1L) stop("Invalid category specification.")
  d <- outdir("cellchat", paste(tissue, category, sep = "_")); days <- if (tissue == "Brain") day_order else day_order[1:3]
  interactions <- score_time_table(tab, c("Source", "Target", "LR_pair", "Ligand", "Receptor", "Pathway"), days)
  lr <- score_time_table(tab, c("LR_pair", "Ligand", "Receptor", "Pathway"), days)
  pathways <- score_time_table(tab, "Pathway", days)
  # Curated functional groups are imported from the final source notebook definitions.
  map_path <- project_path("config", "cellchat_pathway_groups.csv")
  map <- read_csv(map_path)
  for (nm in c("interactions", "lr", "pathways")) {
    value <- get(nm)
    if (nrow(value)) {
      lookup <- map[map$tissue == tissue & map$category == category, ]
      value$Function_Group <- lookup$function_group[match(value$Pathway, lookup$pathway)]
      value$Function_Group[is.na(value$Function_Group)] <- "Other"
    }
    assign(nm, value)
  }
  write_csv(interactions, file.path(d, "Interaction.csv")); write_csv(lr, file.path(d, "LR.csv")); write_csv(pathways, file.path(d, "Pathway.csv"))
  legend <- data.frame(Field = c("Strength", "Interaction", "Baseline", "Change_Score", "Induced_Score", "Peak", "Lost_*", "p threshold", "Self-interactions"),
    Definition = c("Sum of significant communication probabilities; no downsampling", "Unique sender x receiver x CellChat interaction_name",
      "D0 strength > 0", "Baseline: sum abs(log2((later+1e-8)/(D0+1e-8)))", "Induced: row-wise maximum later strength",
      "First time point with maximum strength; chronological tie break", "D0>0 and corresponding later strength=0",
      paste0("p < ", cfg$cellchat_p), "Included; source and target may be the same subtype"))
  openxlsx::write.xlsx(list(Legend = legend, Cell_interactions = interactions, LR_pairs = lr, Pathways = pathways),
                        file.path(outdir("Supplementary_tables"), paste0("Table", spec$table, "_", tissue, "_", category, ".xlsx")), overwrite = TRUE)
  metrics <- tab |> dplyr::group_by(Day) |> dplyr::summarise(interaction_number = dplyr::n(), interaction_weight = sum(Strength), .groups = "drop")
  metrics <- dplyr::left_join(data.frame(Day = days), metrics, by = "Day")
  metrics$interaction_number[is.na(metrics$interaction_number)] <- 0L; metrics$interaction_weight[is.na(metrics$interaction_weight)] <- 0
  write_csv(metrics, file.path(d, "interaction_metrics.csv"))
  if (make_plots) run_source_cellchat_panels(tab, tissue, category, d)
  # All plots consume exactly the same p-filtered tables as the workbook.
  if (make_plots && nrow(tab)) {
    for (type in c("baseline", "induced")) {
      cellchat_heatmap(pathways, "Pathway", days, file.path(d, paste0("Pathway_", type, ".pdf")), type)
      cellchat_heatmap(lr, "LR_pair", days, file.path(d, paste0("LR_", type, ".pdf")), type)
    }
    roles <- dplyr::bind_rows(tab |> dplyr::transmute(Day, Cell = Source, Role = "Sender"), tab |> dplyr::transmute(Day, Cell = Target, Role = "Receiver")) |>
      dplyr::count(Day, Cell, Role, name = "frequency")
    write_csv(roles, file.path(d, "sender_receiver_frequencies.csv"))
    roles$Abbreviation <- abbreviate_cell(roles$Cell)
    p <- ggplot2::ggplot(roles, ggplot2::aes(factor(Day, levels = days), Abbreviation, size = frequency, color = Role)) +
      ggplot2::geom_point() + ggplot2::facet_wrap(~Role) + ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = NULL)
    save_plot(p, file.path(d, "sender_receiver_frequency_bubbles.pdf"), 10, max(6, length(unique(roles$Cell)) * 0.2))
    focus <- read_csv(project_path("config", "cellchat_focus_lr.csv"))
    selected <- unique(focus$LR_pair[focus$tissue == tissue & focus$category == category])
    display_names <- paste(tab$Ligand, tab$Receptor, sep = "_")
    hits <- tolower(tab$LR_pair) %in% tolower(selected) | tolower(display_names) %in% tolower(selected)
    detail <- if (length(selected)) tab[hits, ] else tab
    write_csv(data.frame(selected_LR = selected, detected = tolower(selected) %in% tolower(c(tab$LR_pair, display_names))), file.path(d, "curated_LR_detection_audit.csv"))
    if (nrow(detail)) {
      if (!length(selected)) detail <- detail |> dplyr::group_by(LR_pair) |> dplyr::filter(dplyr::cur_group_id() <= 30) |> dplyr::ungroup()
      detail$Pair <- paste(abbreviate_cell(detail$Source), abbreviate_cell(detail$Target), sep = " → ")
      lookup <- map[map$tissue == tissue & map$category == category, ]
      detail$Function_Group <- lookup$function_group[match(detail$Pathway, lookup$pathway)]
      detail$Function_Group[is.na(detail$Function_Group)] <- "Other"
      detail$Day <- factor(detail$Day, levels = days)
      write_csv(detail, file.path(d, "curated_LR_bubble_data.csv"))
      p <- ggplot2::ggplot(detail, ggplot2::aes(LR_pair, Pair, size = Strength, color = Day, group = Day)) +
        ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.7)) +
        ggplot2::facet_grid(. ~ Function_Group, scales = "free_x", space = "free_x") + ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
      save_plot(p, file.path(d, "curated_LR_cell_pair_bubbles.pdf"), max(12, length(unique(detail$LR_pair)) * 0.35), max(7, length(unique(detail$Pair)) * 0.2))
    }
  }
}
run_cellchat_summary <- function() {
  need(c("dplyr", "tidyr", "openxlsx", "ggplot2")); tab <- get_communication_table(); d <- outdir("cellchat", "global")
  counts <- tab |> dplyr::group_by(Tissue, Category) |>
    dplyr::summarise(Interactions = dplyr::n_distinct(paste(Source, Target, LR_pair, sep = "|")),
      LR_pairs = dplyr::n_distinct(LR_pair), Pathways = dplyr::n_distinct(Pathway), .groups = "drop")
  skeleton <- data.frame(Tissue = table_category_spec$tissue, Category = table_category_spec$category)
  counts <- dplyr::left_join(skeleton, counts, by = c("Tissue", "Category")) |>
    dplyr::mutate(dplyr::across(c(Interactions, LR_pairs, Pathways), ~ tidyr::replace_na(., 0L))) |>
    dplyr::group_by(Tissue) |> dplyr::mutate(Interaction_percentage = if (sum(Interactions) > 0) 100 * Interactions / sum(Interactions) else 0) |> dplyr::ungroup()
  write_csv(counts, file.path(d, "TableS2_category_union_counts.csv"))
  metrics <- tab |> dplyr::group_by(Tissue, Category, Day) |> dplyr::summarise(Number = dplyr::n(), Weight = sum(Strength), .groups = "drop")
  skeleton_day <- do.call(rbind, lapply(seq_len(nrow(skeleton)), function(i) {
    data.frame(Tissue = skeleton$Tissue[i], Category = skeleton$Category[i], Day = if (skeleton$Tissue[i] == "Brain") day_order else day_order[1:3])
  }))
  metrics <- dplyr::left_join(skeleton_day, metrics, by = c("Tissue", "Category", "Day")) |>
    dplyr::mutate(Number = tidyr::replace_na(Number, 0L), Weight = tidyr::replace_na(Weight, 0))
  write_csv(metrics, file.path(d, "category_metrics_by_day.csv"))
  immune_axes <- metrics |> dplyr::filter(Category %in% c("intra_myeloid", "myeloid_lymphoid", "intra_lymphoid")) |>
    dplyr::group_by(Tissue, Day) |> dplyr::mutate(Weight_percentage_within_immune_axes = if (sum(Weight) > 0) 100 * Weight / sum(Weight) else 0) |> dplyr::ungroup()
  write_csv(immune_axes, file.path(d, "immune_axis_strength_percentages.csv"))
  openxlsx::write.xlsx(list(Category_union = counts, Per_timepoint = metrics), file.path(outdir("Supplementary_tables"), "TableS2.xlsx"), overwrite = TRUE)
  p <- ggplot2::ggplot(metrics, ggplot2::aes(factor(Day, levels = day_order), Weight, color = Category, group = Category)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~Tissue, scales = "free_y") + ggplot2::theme_bw() + ggplot2::labs(x = NULL)
  save_plot(p, file.path(d, "global_interaction_weights.pdf"), 11, 6)
  p2 <- ggplot2::ggplot(metrics, ggplot2::aes(factor(Day, levels = day_order), Number, color = Category, group = Category)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~Tissue, scales = "free_y") + ggplot2::theme_bw() + ggplot2::labs(x = NULL)
  save_plot(p2, file.path(d, "global_interaction_numbers.pdf"), 11, 6)
}
run_s9e <- function() {
  need(c("dplyr", "tidyr", "ggplot2")); obj <- get_main(); cc <- get_cellchat()
  # Reuse the tested full-cell/label checks and exact directional extraction from the S9E patch.
  source(project_path("functions", "s9e_core.R"), local = TRUE, encoding = "UTF-8")
  r <- calculate_s9e(obj, cc, p_threshold = cfg$cellchat_p, min_cells = cfg$cellchat_min_cells)
  d <- outdir("FigureS9E"); write_csv(r$data, file.path(d, "raw_strengths_and_proportions.csv"))
  write_csv(r$details, file.path(d, "directional_strength_details.csv")); write_csv(r$input_checks, file.path(d, "input_cell_checks.csv"))
  scaled <- r$data; for (nm in setdiff(names(scaled), "Timepoint")) scaled[[nm]] <- s9e_scale_01(scaled[[nm]])
  write_csv(scaled, file.path(d, "scaled_curves.csv"))
  long <- tidyr::pivot_longer(scaled, -Timepoint, names_to = "Metric", values_to = "Value")
  p <- ggplot2::ggplot(long, ggplot2::aes(factor(Timepoint, levels = day_order), Value, color = Metric, group = Metric)) +
    ggplot2::geom_line(ggplot2::aes(linetype = Metric), linewidth = 1.1) + ggplot2::geom_point(ggplot2::aes(shape = Metric), size = 3) +
    ggplot2::scale_color_manual(values = c(CD8T_proportion = "#E41A1C", Ccl5_Ccr5 = "#4DAF4A", Cd86_Cd28 = "#FF7F00", Microglia1_proportion = "#377EB8", Microglia2_proportion = "#984EA3")) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + ggplot2::labs(x = NULL, y = "Scaled value (0–1)")
  save_plot(p, file.path(d, "CD8T_Microglia_CCL5CCR5_CD86CD28_convergence.pdf"), 8, 6)
}
