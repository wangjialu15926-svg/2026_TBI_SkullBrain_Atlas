tdeg_key <- function(lineage, tissue) paste(lineage, tissue, sep = "_")
fit_tdeg_gene <- function(gene, values, time_normalized, gam_k, n_points = 100L) {
  need("mgcv")
  if (length(values) != length(time_normalized)) stop("Expression and time vectors differ in length.")
  if (length(gam_k) != 1L || !is.finite(gam_k) || gam_k < 3 || gam_k != as.integer(gam_k)) stop("Invalid GAM basis dimension.")
  if (length(n_points) != 1L || !is.finite(n_points) || n_points < 2 || n_points != as.integer(n_points)) stop("Invalid GAM prediction-grid length.")
  dat <- data.frame(expression = as.numeric(values), time_normalized = as.numeric(time_normalized))
  dat <- dat[is.finite(dat$expression) & is.finite(dat$time_normalized), , drop = FALSE]
  if (nrow(dat) < 4L) return(list(gene = gene, success = FALSE, reason = "insufficient_data", warnings = character()))
  notes <- character(); s <- mgcv::s
  tryCatch({
    fit <- withCallingHandlers(mgcv::gam(expression ~ s(time_normalized, bs = "cr", k = gam_k), data = dat, method = "REML"),
      warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") })
    if (isFALSE(fit$converged)) stop("GAM did not converge")
    pred <- data.frame(time_normalized = seq(0, 1, length.out = n_points))
    est <- stats::predict(fit, newdata = pred, se.fit = TRUE)
    pred$fitted <- as.numeric(est$fit); pred$se <- as.numeric(est$se.fit)
    if (any(!is.finite(pred$fitted)) || any(!is.finite(pred$se))) stop("Non-finite GAM prediction or uncertainty")
    pred$lower95 <- pred$fitted - 1.96 * pred$se; pred$upper95 <- pred$fitted + 1.96 * pred$se
    list(gene = gene, success = TRUE, range = diff(range(pred$fitted)), predictions = pred, n_samples = nrow(dat), reason = "", warnings = notes)
  }, error = function(e) list(gene = gene, success = FALSE, reason = conditionMessage(e), warnings = notes))
}
read_tdeg_result <- function(lineage, tissue) {
  key <- tdeg_key(lineage, tissue)
  stages <- c(Macrophage_Brain = "8", Macrophage_Skull = "9", Neutrophil_Brain = "11", Neutrophil_Skull = "12",
              Microglia_Brain = "15", Bcell_Brain = "18", Bcell_Skull = "19", Tcell_Brain = "20", Tcell_Skull = "21")
  if (!key %in% names(stages)) stop("Unknown tDEG lineage/compartment: ", key)
  assert_no_incomplete_stages(c("3", unname(stages[key])))
  path <- project_path("results", "tDEG", tdeg_key(lineage, tissue), "tDEG_final.rds")
  value <- read_input(path)
  if (!identical(value$key, key)) stop("tDEG file identity does not match requested lineage/compartment: ", key)
  if (!identical(value$main_md5, hash_input(main_final_path()))) stop("Stale tDEG main input: ", value$key)
  source <- project_path("results", "objects", paste0(lineage, "_final.rds"))
  if (!identical(value$lineage_md5, hash_input(source))) stop("Stale or missing tDEG subtype fingerprint; rerun ", tdeg_key(lineage, tissue))
  expected <- cfg[c("tdeg_de_p", "tdeg_de_abs_fc", "tdeg_min_pct", "tdeg_logfc_threshold", "range_quantile", "prediction_points", "seed", "module_score_seed", "module_score_nbin", "module_score_ctrl")]
  if (!identical(value$analysis_settings, expected)) stop("tDEG settings changed; rerun ", tdeg_key(lineage, tissue))
  value
}
run_tdeg <- function(lineage, tissue, n_modules, gam_k, figure) {
  need(c("Seurat", "mgcv", "dplyr", "tidyr", "ggplot2", "pheatmap"))
  s <- mgcv::s
  obj <- get_lineage(lineage, tissue); key <- tdeg_key(lineage, tissue); d <- outdir("tDEG", key)
  module_order <- cfg$tdeg_module_order[[key]]
  if (!setequal(module_order, paste0("C", seq_len(n_modules))) || anyDuplicated(module_order)) stop("Module display order does not cover the requested modules in ", key)
  days <- intersect(day_order, unique(as.character(obj$Day)))
  expected <- if (tissue == "Brain") day_order else day_order[1:3]
  if (!identical(days, expected)) stop("Unexpected time points for ", key)
  Seurat::Idents(obj) <- factor(obj$Day, levels = days)
  pair_results <- list()
  pairs <- utils::combn(days, 2, simplify = FALSE)
  for (pair in pairs) {
    tab <- Seurat::FindMarkers(obj, ident.1 = pair[1], ident.2 = pair[2], test.use = "wilcox",
                                logfc.threshold = cfg$tdeg_logfc_threshold, min.pct = cfg$tdeg_min_pct, only.pos = FALSE, random.seed = cfg$seed)
    if (!nrow(tab)) next
    tab$gene <- rownames(tab); tab$comparison <- paste(pair, collapse = "_vs_")
    tab$up.down <- ifelse(tab$avg_log2FC > 0, "up", "down"); pair_results[[paste(pair, collapse = "_vs_")]] <- tab
  }
  if (!length(pair_results)) stop("No pairwise DEGs in ", key)
  deg_all <- dplyr::bind_rows(pair_results)
  up_all <- deg_all |> dplyr::filter(!is.na(p_val_adj), p_val_adj < cfg$tdeg_de_p, abs(avg_log2FC) > cfg$tdeg_de_abs_fc)
  candidates <- unique(up_all$gene)
  if (!length(candidates)) stop("No pairwise-DEG candidates in ", key)
  write_csv(deg_all, file.path(d, "pairwise_DEG_all.csv")); write_csv(up_all, file.path(d, "pairwise_DEG_significant.csv"))
  # Stable sample keys prevent Seurat separator changes from corrupting Day/Lib parsing.
  sample_info <- unique(obj@meta.data[, c("Day", "Lib")])
  sample_info <- sample_info[order(match(sample_info$Day, days), sample_info$Lib), ]
  sample_info$sample_id <- sprintf("S%04d", seq_len(nrow(sample_info)))
  lookup <- paste(sample_info$Day, sample_info$Lib, sep = "\t")
  obj$.tbi_sample <- sample_info$sample_id[match(paste(obj$Day, obj$Lib, sep = "\t"), lookup)]
  means <- Seurat::AverageExpression(obj, group.by = ".tbi_sample", assays = "RNA", slot = "data", return.seurat = FALSE)$RNA
  sample_info <- sample_info[match(colnames(means), sample_info$sample_id), ]
  if (anyNA(sample_info$sample_id)) stop("Pseudo-bulk columns do not match library metadata.")
  days_numeric <- c(D0 = 0, D1 = 1, D7 = 7, M6 = 180)[days]
  time <- log10(days_numeric + 1); time <- (time - min(time)) / (max(time) - min(time))
  sample_info$time_normalized <- unname(time[as.character(sample_info$Day)])
  write_csv(sample_info, file.path(d, "pseudobulk_sample_info.csv"))
  means <- means[intersect(candidates, rownames(means)), , drop = FALSE]
  fits <- lapply(rownames(means), function(gene) fit_tdeg_gene(gene, means[gene, ], sample_info$time_normalized, gam_k, cfg$prediction_points))
  diagnostics <- data.frame(gene = vapply(fits, `[[`, character(1), "gene"), success = vapply(fits, `[[`, logical(1), "success"),
    reason = vapply(fits, `[[`, character(1), "reason"), warnings = vapply(fits, function(x) paste(unique(x$warnings), collapse = "; "), character(1)))
  write_csv(diagnostics, file.path(d, "GAM_fit_diagnostics.csv"))
  successful <- Filter(function(x) isTRUE(x$success), fits)
  failed <- Filter(function(x) !isTRUE(x$success), fits)
  if (length(failed)) write_csv(data.frame(gene = vapply(failed, `[[`, character(1), "gene"), reason = vapply(failed, `[[`, character(1), "reason")), file.path(d, "GAM_failed_genes.csv"))
  if (!length(successful)) stop("No successful GAM fits in ", key, "; check the number of independent library profiles.")
  ranges <- data.frame(gene = vapply(successful, `[[`, character(1), "gene"), range = vapply(successful, `[[`, numeric(1), "range"))
  threshold <- unname(stats::quantile(ranges$range, cfg$range_quantile))
  retained <- ranges$gene[ranges$range > threshold]
  successful <- Filter(function(x) x$gene %in% retained, successful)
  if (length(successful) <= n_modules) stop("Too few dynamic genes for requested modules in ", key)
  write_csv(ranges, file.path(d, "gene_ranges_statistics.csv"))
  # Candidate restriction plus successful-fit/range filter implements the two criteria.
  trend <- vapply(successful, function(x) x$predictions$fitted, numeric(cfg$prediction_points))
  colnames(trend) <- vapply(successful, `[[`, character(1), "gene")
  scaled <- apply(trend, 2, function(x) (x - min(x)) / diff(range(x)))
  hc <- stats::hclust(stats::dist(t(scaled), method = "euclidean"), method = "ward.D2")
  clusters <- stats::cutree(hc, k = n_modules)
  modules <- data.frame(gene = names(clusters), cluster = paste0("C", clusters))
  modules <- dplyr::left_join(modules, ranges, by = "gene")
  write_csv(modules, file.path(d, "module_genes.csv"))
  curves <- dplyr::bind_rows(lapply(successful, function(x) {
    p <- x$predictions; r <- x$range
    data.frame(gene = x$gene, time = p$time_normalized, fitted = (p$fitted - min(p$fitted)) / r,
                 lower95 = pmax(0, pmin(1, (p$lower95 - min(p$fitted)) / r)),
                 upper95 = pmax(0, pmin(1, (p$upper95 - min(p$fitted)) / r)))
  })) |> dplyr::left_join(modules[, c("gene", "cluster")], by = "gene")
  average <- curves |> dplyr::group_by(cluster, time) |>
    dplyr::summarise(mean = mean(fitted), lower95 = mean(lower95), upper95 = mean(upper95), .groups = "drop") |>
    dplyr::mutate(se_display = (upper95 - lower95) / (2 * 1.96),
                   lower50 = mean - 0.674 * se_display, upper50 = mean + 0.674 * se_display,
                   lower75 = mean - 1.150 * se_display, upper75 = mean + 1.150 * se_display)
  average$cluster <- factor(average$cluster, levels = module_order)
  write_csv(average, file.path(d, "module_mean_display_bands.csv"))
  p <- ggplot2::ggplot(average, ggplot2::aes(time, mean, color = cluster, fill = cluster)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower75, ymax = upper75), alpha = 0.18, color = NA) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower50, ymax = upper50), alpha = 0.35, color = NA) +
    ggplot2::geom_line(linewidth = 1) + ggplot2::facet_wrap(~cluster, ncol = 1) +
    ggplot2::scale_color_manual(values = module_colors) + ggplot2::scale_fill_manual(values = module_colors) +
    ggplot2::scale_x_continuous(breaks = unname(time), labels = names(time)) + ggplot2::theme_bw() +
    ggplot2::labs(x = "Time", y = "Scaled fitted expression", title = key,
                    subtitle = "Range > candidate fitted-range Q25; 50%/75% display bands") + ggplot2::theme(legend.position = "none")
  save_plot(p, file.path(d, paste0(figure, "_GAM_modules_gradient_bands.pdf")), 7, max(7, n_modules * 1.6))
  # Preserve the source display-band aggregation; these are not a new module-level inferential CI.
  Seurat::Idents(obj) <- obj$level3
  score_columns <- character()
  for (cl in paste0("C", seq_len(n_modules))) {
    genes <- modules$gene[modules$cluster == cl]
    obj <- Seurat::AddModuleScore(obj, features = list(genes), name = paste0(cl, "_score"),
                                   nbin = cfg$module_score_nbin, ctrl = cfg$module_score_ctrl, seed = cfg$module_score_seed)
    score_columns <- c(score_columns, paste0(cl, "_score1"))
  }
  score <- stats::aggregate(obj@meta.data[, score_columns, drop = FALSE], list(level3 = obj$level3), mean)
  score <- score[match(intersect(subtype_orders[[lineage]], as.character(score$level3)), as.character(score$level3)), , drop = FALSE]
  write_csv(score, file.path(d, "module_scores_by_subtype.csv"))
  score_matrix <- as.matrix(score[, paste0(module_order, "_score1"), drop = FALSE]); rownames(score_matrix) <- score$level3
  pheatmap::pheatmap(t(score_matrix), cluster_rows = FALSE, cluster_cols = FALSE,
                      filename = file.path(d, "module_scores_by_subtype.pdf"), width = 11, height = 6)
  enrich_gene_groups(split(modules$gene, modules$cluster), d, universe = candidates)
  # Representative-gene panels use reproducible top-range selection; manuscript overrides are configurable.
  selections_path <- project_path("config", "representative_genes.csv")
  selections <- read_csv(selections_path)
  chosen_all <- selections[selections$key == key & nzchar(selections$gene), ]
  chosen <- unique(chosen_all[, c("key", "cluster", "gene"), drop = FALSE])
  if (nrow(chosen)) {
    chosen$status <- ifelse(chosen$gene %in% modules$gene, "retained_tDEG", "not_retained_after_recomputation")
    chosen$actual_cluster <- modules$cluster[match(chosen$gene, modules$gene)]
    if (any(chosen$cluster != "AUTO" & chosen$cluster != chosen$actual_cluster, na.rm = TRUE)) stop("Explicit representative gene module labels changed in ", key)
    write_csv(chosen, file.path(d, "source_selected_gene_retention_audit.csv"))
    chosen <- chosen[chosen$status == "retained_tDEG", , drop = FALSE]
  } else {
    chosen <- modules |> dplyr::group_by(cluster) |> dplyr::slice_max(range, n = 2, with_ties = FALSE) |> dplyr::ungroup()
    writeLines("No manuscript gene override supplied: panels use the two largest-range genes per module. Review against the manuscript.", file.path(d, "representative_gene_selection_NOTE.txt"))
  }
  selected_curves <- curves[curves$gene %in% chosen$gene, ]
  if (nrow(selected_curves)) {
    p <- ggplot2::ggplot(selected_curves, ggplot2::aes(time, fitted, color = gene)) + ggplot2::geom_line() +
      ggplot2::facet_wrap(~cluster, scales = "free_y") + ggplot2::scale_x_continuous(breaks = unname(time), labels = names(time)) + ggplot2::theme_bw()
    save_plot(p, file.path(d, "representative_gene_fitted_curves.pdf"), 11, 8)
  }
  write_csv(chosen, file.path(d, "representative_gene_selection.csv"))
  if (nrow(chosen_all)) run_source_gene_panels(obj, chosen_all, d)
  saveRDS(list(modules = modules, successful_fits = successful, candidates = candidates, sample_info = sample_info,
               range_threshold = threshold, hclust = hc, key = key, n_modules = n_modules, gam_k = gam_k,
               main_md5 = hash_input(main_final_path()), lineage_md5 = hash_input(project_path("results", "objects", paste0(lineage, "_final.rds"))),
               analysis_settings = cfg[c("tdeg_de_p", "tdeg_de_abs_fc", "tdeg_min_pct", "tdeg_logfc_threshold", "range_quantile", "prediction_points", "seed", "module_score_seed", "module_score_nbin", "module_score_ctrl")], final_genes = modules$gene,
               candidate_definition = "union pairwise adjusted-p<0.05 and abs-log2FC>0.5; successful GAM and fitted-range>Q25"), file.path(d, "tDEG_final.rds"))
  audit_run(paste0("tDEG_", key), c(main_final_path(), project_path("results", "objects", paste0(lineage, "_final.rds"))), cfg)
}
build_module_overlap <- function(brain, skull, corresponding) {
  check_columns(brain, c("gene", "cluster"), "Brain modules")
  check_columns(skull, c("gene", "cluster"), "Skull modules")
  check_columns(corresponding, c("brain", "skull"), "Corresponding module pairs")
  for (tab in list(brain, skull)) {
    if (!nrow(tab) || anyNA(tab[, c("gene", "cluster")]) || any(!nzchar(tab$gene)) || anyDuplicated(tab$gene)) stop("Each gene must have one module within a compartment.")
  }
  if (anyNA(corresponding[, c("brain", "skull")]) || anyDuplicated(corresponding[, c("brain", "skull")]) ||
      any(!corresponding$brain %in% brain$cluster) || any(!corresponding$skull %in% skull$cluster)) stop("Invalid corresponding module pair.")
  bs <- split(as.character(brain$gene), as.character(brain$cluster))
  ss <- split(as.character(skull$gene), as.character(skull$cluster))
  sets <- c(setNames(bs, paste0("Brain_", names(bs))), setNames(ss, paste0("Skull_", names(ss))))
  genes <- unique(unlist(sets, use.names = FALSE)); membership <- data.frame(gene = genes)
  for (nm in names(sets)) membership[[nm]] <- genes %in% sets[[nm]]
  membership$intersect_type <- "Other"
  summary <- list(); shared <- list()
  for (i in seq_len(nrow(corresponding))) {
    pair <- corresponding[i, ]; bg <- bs[[pair$brain]]; sg <- ss[[pair$skull]]
    common <- intersect(bg, sg)
    membership$intersect_type[membership$gene %in% common] <- "Corresponding"
    summary[[i]] <- data.frame(brain = pair$brain, skull = pair$skull, n_brain = length(bg), n_skull = length(sg),
      n_shared = length(common), fraction_brain = length(common)/length(bg), fraction_skull = length(common)/length(sg),
      jaccard = length(common)/length(union(bg, sg)))
    shared[[i]] <- data.frame(brain = rep(pair$brain, length(common)), skull = rep(pair$skull, length(common)), gene = common)
  }
  list(membership = membership, sets = names(sets), summary = do.call(rbind, summary), shared = do.call(rbind, shared))
}
run_module_overlap <- function() {
  need(c("ComplexUpset", "ggplot2")); pairs <- c("Macrophage", "Neutrophil", "Bcell", "Tcell")
  for (lineage in pairs) {
    a <- read_tdeg_result(lineage, "Brain")
    b <- read_tdeg_result(lineage, "Skull")
    if (a$main_md5 != hash_input(main_final_path()) || b$main_md5 != a$main_md5) stop("Stale tDEG inputs in overlap: ", lineage)
    config <- read_csv(project_path("config", "module_overlap_pairs.csv"))
    check_columns(config, c("lineage", "brain", "skull"), "module_overlap_pairs")
    corresponding <- config[config$lineage == lineage, , drop = FALSE]
    if (!nrow(corresponding)) stop("No trend correspondence configured for ", lineage)
    overlap <- build_module_overlap(a$modules, b$modules, corresponding)
    d <- outdir("module_overlap", lineage)
    write_csv(overlap$membership, file.path(d, "gene_module_membership.csv"))
    write_csv(corresponding, file.path(d, "corresponding_module_pairs_used.csv"))
    write_csv(overlap$summary, file.path(d, "corresponding_module_overlap_summary.csv"))
    write_csv(overlap$shared, file.path(d, "corresponding_module_shared_genes.csv"))
    write_csv(all_module_pair_gene_overlaps(a$modules, b$modules, corresponding), file.path(d, "source_all_module_pairs_overlap_genes.csv"))
    p <- ComplexUpset::upset(overlap$membership, intersect = overlap$sets, min_size = 1,
      base_annotations = list("Intersection size" = ComplexUpset::intersection_size(
        mapping = ggplot2::aes(fill = intersect_type)) + ggplot2::scale_fill_manual(
          values = c(Corresponding = "#e07b54", Other = "grey60"), name = "Trend correspondence")))
    save_plot(p, file.path(d, "Brain_Skull_module_UpSet.pdf"), 15, 8)
  }
}
