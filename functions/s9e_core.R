# Figure S9E: compute communication strengths and proportions from input objects.
# Run with the final full Brain Seurat object and the non-downsampled CellChat list.

calculate_s9e <- function(seurat_obj, cc_list, p_threshold = 0.05, min_cells = 50L) {
  timepoints <- c("D0", "D1", "D7", "M6")
  groups <- c("CD8 T Cell", "Microglia1", "Microglia2")
  meta <- seurat_obj@meta.data
  required <- c("Group", "Day", "level3")
  if (!all(required %in% colnames(meta))) {
    stop("S9E requires Group, Day and level3 in the final Seurat metadata.")
  }
  if (!is.list(cc_list) || !all(paste0("Brain_", timepoints) %in% names(cc_list))) {
    stop("Load the non-downsampled CellChat list with Brain_D0/D1/D7/M6 first.")
  }
  brain <- meta[!is.na(meta$Group) & as.character(meta$Group) == "Brain", , drop = FALSE]
  brain$Day <- as.character(brain$Day)
  brain$level3 <- as.character(brain$level3)
  if (anyNA(brain[, c("Day", "level3")]) || anyDuplicated(rownames(brain))) {
    stop("Brain metadata has missing Day/level3 values or duplicated cell IDs.")
  }
  if (!all(timepoints %in% brain$Day)) stop("Brain metadata lacks a required time point.")
  if (!all(groups %in% brain$level3)) stop("Check final CD8 T/Microglia1/Microglia2 labels.")

  # Denominator: all final Brain cells at each Day, before the CellChat min50 filter.
  totals <- vapply(timepoints, function(day) sum(brain$Day == day), integer(1))
  proportions <- vapply(groups, function(group) {
    vapply(timepoints, function(day) {
      100 * sum(brain$Day == day & brain$level3 == group) / totals[[day]]
    }, numeric(1))
  }, numeric(length(timepoints)))
  cell_prop <- data.frame(Day = timepoints, proportions)
  names(cell_prop)[-1] <- groups

  directions <- data.frame(
    Metric = c("Ccl5_Ccr5", "Ccl5_Ccr5", "Cd86_Cd28", "Cd86_Cd28"),
    Source = c("CD8 T Cell", "CD8 T Cell", "Microglia1", "Microglia2"),
    Target = c("Microglia1", "Microglia2", "CD8 T Cell", "CD8 T Cell"),
    Pair = c("CCL5_CCR5", "CCL5_CCR5", "CD86_CD28", "CD86_CD28"),
    stringsAsFactors = FALSE
  )
  details <- list()
  input_checks <- list()
  for (day in timepoints) {
    cc <- cc_list[[paste0("Brain_", day)]]
    cc_meta <- cc@meta
    cc_cells <- rownames(cc_meta)
    day_meta <- brain[brain$Day == day, , drop = FALSE]
    if (length(cc_cells) == 0L || anyDuplicated(cc_cells)) {
      stop("Invalid CellChat cell IDs for Brain_", day)
    }
    if (!all(cc_cells %in% rownames(day_meta))) {
      stop("CellChat cells do not match final Brain metadata at ", day,
           ". Load objects built from the same final annotations.")
    }
    if (!"level3" %in% colnames(cc_meta)) stop("CellChat metadata lacks level3 at ", day)
    cc_labels <- as.character(cc_meta$level3)
    final_labels <- day_meta[cc_cells, "level3"]
    if (anyNA(cc_labels) || !identical(cc_labels, unname(final_labels))) {
      stop("CellChat and final Seurat level3 labels disagree at ", day)
    }
    # The min50 construction may exclude entire subtypes; it must retain every
    # cell within each included subtype. A sampled list therefore fails here.
    type_counts <- table(day_meta$level3)
    expected_types <- names(type_counts)[type_counts >= min_cells]
    if (!setequal(unique(cc_labels), expected_types)) stop("CellChat omitted an eligible subtype or retained an undersized subtype at ", day)
    expected_cells <- rownames(day_meta)[day_meta$level3 %in% expected_types]
    if (!setequal(cc_cells, expected_cells)) {
      stop("CellChat omits cells within included subtypes at ", day,
           ". Use the non-downsampled list from the full CellChat build (numbered script 4).")
    }
    idents <- as.character(cc@idents)
    if (length(idents) != length(cc_cells)) stop("Invalid CellChat idents at ", day)
    if (!is.null(names(cc@idents))) {
      if (!setequal(names(cc@idents), cc_cells)) stop("CellChat idents IDs disagree at ", day)
      idents <- as.character(cc@idents[cc_cells])
    }
    if (!identical(idents, cc_labels)) stop("CellChat idents differ from level3 at ", day)
    input_checks[[day]] <- data.frame(
      Timepoint = day, Brain_cells = nrow(day_meta),
      CellChat_cells = length(cc_cells), Included_subtypes = length(unique(cc_labels)),
      Full_cells_within_included_subtypes = TRUE
    )
    if (!exists("extract_sig_communications", mode = "function")) stop("Load functions/cellchat.R before the S9E core.")
    communication <- extract_sig_communications(cc, "Brain", day, p_threshold = p_threshold)
    for (i in seq_len(nrow(directions))) {
      spec <- directions[i, ]
      present <- all(c(spec$Source, spec$Target) %in% cc_labels)
      if (!present || nrow(communication) == 0L) {
        selected <- communication[FALSE, , drop = FALSE]
      } else {
        selected <- communication[
          as.character(communication$Source) == spec$Source &
          as.character(communication$Target) == spec$Target &
          tolower(as.character(communication$LR_pair)) == tolower(spec$Pair), , drop = FALSE]
      }
      strength <- if (nrow(selected) == 0L) 0 else sum(selected$Strength)
      if (!is.finite(strength) || strength < 0) stop("Invalid strength at ", day)
      details[[length(details) + 1L]] <- data.frame(
        Timepoint = day, Metric = spec$Metric, Source = spec$Source,
        Target = spec$Target, Pair = spec$Pair, Strength = strength,
        N_significant_entries = nrow(selected),
        Status = if (!present) "group_not_in_CellChat" else if (nrow(selected) == 0L)
          "no_significant_entry" else "detected", stringsAsFactors = FALSE
      )
    }
  }
  details <- do.call(rbind, details)
  result <- data.frame(Timepoint = timepoints, stringsAsFactors = FALSE)
  for (metric in c("Ccl5_Ccr5", "Cd86_Cd28")) {
    result[[metric]] <- vapply(timepoints, function(day) {
      sum(details$Strength[details$Timepoint == day & details$Metric == metric])
    }, numeric(1))
  }
  result$CD8T_proportion <- cell_prop[["CD8 T Cell"]]
  result$Microglia1_proportion <- cell_prop[["Microglia1"]]
  result$Microglia2_proportion <- cell_prop[["Microglia2"]]
  list(data = result, cell_prop = cell_prop, details = details,
       input_checks = do.call(rbind, input_checks))
}

s9e_scale_01 <- function(x) {
  if (any(!is.finite(x))) stop("S9E scaling input must be finite.")
  limits <- range(x)
  if (limits[1] == limits[2]) return(rep(0, length(x)))
  (x - limits[1]) / (limits[2] - limits[1])
}

