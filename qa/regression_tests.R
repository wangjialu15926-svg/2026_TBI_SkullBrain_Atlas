# Core regression checks. No Seurat/CellChat model is fitted by this harness.
suppressWarnings(Sys.setlocale("LC_CTYPE", if (.Platform$OS.type == "windows") "English_United States.utf8" else "C.UTF-8"))
args <- commandArgs(trailingOnly = TRUE)
TBI_ROOT <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
source(file.path(TBI_ROOT, "config", "settings.R"), encoding = "UTF-8")
for (module in c("common", "preprocessing", "plots", "figure1_panels", "tdeg", "cellchat", "s9e_core", "trajectory", "validation"))
  source(file.path(TBI_ROOT, "functions", paste0(module, ".R")), encoding = "UTF-8")
bundle_root <- TBI_ROOT
output_base <- ensure_dir(file.path(bundle_root, "qa", "test_results"))
output <- ensure_dir(file.path(output_base, paste0("run_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d%H%M%S"))))
results <- list()
test <- function(name, code) {
  error <- tryCatch({ force(code); NULL }, error = function(e) conditionMessage(e))
  results[[length(results) + 1L]] <<- data.frame(test = name, pass = is.null(error), detail = if (is.null(error)) "PASS" else error)
  cat(if (is.null(error)) "PASS " else "FAIL ", name, if (!is.null(error)) paste(":", error) else "", "\n", sep = "")
}
equals <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE)))
errors <- function(code, pattern = NULL) {
  message <- tryCatch({ force(code); NULL }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(message)); if (!is.null(pattern)) stopifnot(grepl(pattern, message, fixed = TRUE))
}
test("all active R sources parse", {
  files <- c(list.files(bundle_root, "\\.R$", full.names = TRUE), list.files(file.path(bundle_root, "functions"), "\\.R$", full.names = TRUE), file.path(bundle_root, "config", "settings.R"))
  invisible(lapply(files, parse, encoding = "UTF-8"))
})
test("Greek and short labels canonicalize", equals(canonical_labels(c("Mφ12", "Neu9", "MG1", "γδ T", "B1-like MBC")), c("Macrophage12", "Neutrophil9", "Microglia1", "γδ T Cell", "B1 Cell")))
test("missing input fingerprints fail", errors(hash_input(file.path(tempdir(), "does_not_exist.rds"))))
test("input paths preserve drive and UNC roots", {
  equals(input_path("F:/data/a.rds"), "F:/data/a.rds")
  equals(input_path("\\\\server\\data\\a.rds"), "\\\\server\\data\\a.rds")
  errors(input_path(NA_character_)); errors(input_path(character()))
})
test("UTF8 CSV read preserves labels", {
  path <- file.path(output, "unicode.csv"); write_csv(data.frame(label = c("γδ T Cell", "Mφ12")), path)
  equals(read_csv(path)$label, c("γδ T Cell", "Mφ12"))
})

meta <- data.frame(Lib = c("TBI_A", "TBI_A", "TBI_B", "TBI_B"), Day = c("D0", "D0", "D1", "D1"), Group = "Brain", Source = "In-house", Mice_per_library = 4L, state = c("X", "X", "X", "Y"))
manifest <- rbind(unique(meta[, 1:5]), data.frame(Lib = "TBI_C", Day = "D7", Group = "Brain", Source = "In-house", Mice_per_library = 4L))
test("proportions distinguish zero state from absent parent", {
  p <- build_library_proportions(meta, "state", c("X", "Y"), manifest)
  equals(p$percentage[p$Lib == "TBI_A" & p$cell_type == "Y"], 0)
  stopifnot(all(is.na(p$percentage[p$Lib == "TBI_C"])), all(p$status[p$Lib == "TBI_C"] == "no_parent_lineage_cells"))
  equals(p$percentage[p$Lib == "TBI_B"], c(50, 50)); equals(unique(p$denominator[p$Lib == "TBI_A"]), 2L)
})
test("proportion library metadata mismatches fail", {
  bad <- meta; bad$Day[1] <- "D7"; errors(build_library_proportions(bad, "state", c("X", "Y"), manifest), "disagree")
  errors(build_library_proportions(meta, "state", "X", manifest), "cover")
  errors(build_library_proportions(meta, "state", c("X", "Y"), rbind(manifest, transform(manifest[1, ], Day = "D7"))), "exactly one")
})
test("caption n counts libraries and reports absent parent libraries", {
  grid <- build_library_proportions(meta, "state", c("X", "Y"), manifest)
  n <- summarize_library_coverage(grid)
  equals(n$n_libraries[match(c("D0", "D1", "D7"), n$Day)], c(1L, 1L, 0L))
  equals(n$n_unavailable[match(c("D0", "D1", "D7"), n$Day)], c(0L, 0L, 1L))
  equals(n$included_Lib[n$Day == "D0"], "TBI_A"); equals(n$unavailable_Lib[n$Day == "D7"], "TBI_C")
  bad <- grid; bad$denominator[1] <- 99; errors(summarize_library_coverage(bad), "Inconsistent")
})

long <- data.frame(ID = c("baseline", "baseline", "baseline", "induced", "induced", "lost"), Day = c("D0", "D1", "D1", "D1", "D7", "D0"), Strength = c(1, 1, 2, 3, 3, 2))
test("communication aggregation and chronological peak", {
  z <- score_time_table(long, "ID", c("D0", "D1", "D7")); a <- z[z$ID == "baseline", ]; b <- z[z$ID == "induced", ]; c <- z[z$ID == "lost", ]
  equals(a$D1, 3); equals(a$Change_Score, abs(log2((3 + cfg$cellchat_epsilon)/(1 + cfg$cellchat_epsilon))) + abs(log2(cfg$cellchat_epsilon/(1 + cfg$cellchat_epsilon))))
  equals(b$Induced_Score, 3); equals(b$Peak, "D1"); stopifnot(c$Lost_all_later, a$Lost_D7, !b$Baseline)
})
test("communication singleton and empty schemas", {
  equals(score_time_table(long[1, ], "ID", c("D0", "D1"))$Lost_D1, TRUE)
  z <- score_time_table(long[FALSE, ], "ID", day_order); stopifnot(nrow(z) == 0, all(c("Peak", "M6", "Lost_M6") %in% names(z)))
})
test("invalid communication strengths and days fail", {
  for (value in c(-1, NA_real_, Inf)) { bad <- long; bad$Strength[1] <- value; errors(score_time_table(bad, "ID", day_order)) }
  errors(score_time_table(long, "ID", c("D1", "D0"))); errors(score_time_table(long, "ID", c("D0", "D1")))
})
methods::setClass("TBIChatTest", slots = c(net = "list", LR = "list", meta = "data.frame", idents = "factor"))
methods::setClass("TBIObjectTest", slots = c(meta.data = "data.frame", misc = "list", mat = "matrix"))
methods::setMethod("dim", "TBIObjectTest", function(x) dim(x@mat))
methods::setMethod("dimnames", "TBIObjectTest", function(x) dimnames(x@mat))
methods::setMethod("$", "TBIObjectTest", function(x, name) x@meta.data[[name]])
methods::setReplaceMethod("$", "TBIObjectTest", function(x, name, value) { x@meta.data[[name]] <- value; x })
subset.TBIObjectTest <- function(x, cells, ...) { x@mat <- x@mat[, cells, drop = FALSE]; x@meta.data <- x@meta.data[cells, , drop = FALSE]; x }
object <- function(labels) {
  meta <- data.frame(level3 = labels, row.names = paste0("c", seq_along(labels)))
  methods::new("TBIObjectTest", meta.data = meta, misc = list(), mat = matrix(0, 1, nrow(meta), dimnames = list("g", rownames(meta))))
}
groups <- c("CD8 T Cell", "Microglia1", "Microglia2")
brain <- do.call(rbind, lapply(day_order, function(day) {
  data.frame(Group = "Brain", Day = day, level3 = rep(c(groups, "Tiny"), c(50, 50, 50, 3)), row.names = paste0(day, "_", seq_len(153)))
}))
main <- methods::new("TBIObjectTest", meta.data = brain, misc = list(), mat = matrix(0, 1, nrow(brain), dimnames = list("g", rownames(brain))))
make_chat <- function(day) {
  m <- brain[brain$Day == day & brain$level3 != "Tiny", , drop = FALSE]
  dn <- list(groups, groups, c("CCL5_CCR5", "CD86_CD28"))
  prob <- array(0, c(3, 3, 2), dimnames = dn); pval <- array(1, c(3, 3, 2), dimnames = dn)
  prob[1, 2, 1] <- 0.2; pval[1, 2, 1] <- .01
  prob[1, 3, 1] <- 0.4; pval[1, 3, 1] <- .01
  prob[2, 1, 2] <- 0.3; pval[2, 1, 2] <- .01
  prob[3, 1, 2] <- 0.9; pval[3, 1, 2] <- .05
  lr <- data.frame(ligand = c("Ccl5", "Cd86"), receptor = c("Ccr5", "Cd28"), pathway_name = c("CCL", "CD86"), row.names = dn[[3]])
  methods::new("TBIChatTest", net = list(prob = prob, pval = pval), LR = list(LRsig = lr), meta = m, idents = factor(m$level3, levels = groups))
}
chats <- setNames(lapply(day_order, make_chat), paste0("Brain_", day_order))
test("strict p mask, LR IDs and directional S9E strengths", {
  z <- calculate_s9e(main, chats); equals(z$data$Ccl5_Ccr5, rep(.6, 4)); equals(z$data$Cd86_Cd28, rep(.3, 4))
  equals(z$data$CD8T_proportion, rep(100 * 50/153, 4)); equals(nrow(z$details), 16L)
})
test("zero significant interactions are typed and safe", {
  empty <- chats; empty[[1]]@net$pval[] <- 1
  z <- extract_sig_communications(empty[[1]], "Brain", "D0"); stopifnot(nrow(z) == 0, is.numeric(z$Strength))
  equals(calculate_s9e(main, empty)$data$Ccl5_Ccr5[1], 0)
})
test("1x1 CellChat arrays retain dimensions", {
  cc <- chats[[1]]; cc@net$prob <- array(.2, c(1, 1, 1), dimnames = list("A", "B", "CCL5_CCR5")); cc@net$pval <- cc@net$prob * 0 + .01
  equals(extract_sig_communications(cc, "Brain", "D0")$Strength, .2)
})
test("probability and pvalue name mismatch fails", {
  cc <- chats[[1]]; dimnames(cc@net$pval)[[1]] <- rev(groups); errors(extract_sig_communications(cc, "Brain", "D0"), "dimension names")
})
test("within-state downsampling rejected", {
  bad <- chats; bad[[1]]@meta <- bad[[1]]@meta[-1, , drop = FALSE]; bad[[1]]@idents <- bad[[1]]@idents[-1]
  errors(calculate_s9e(main, bad), "omits cells")
})
test("omitting an entire eligible subtype rejected", {
  bad <- chats; keep <- bad[[1]]@meta$level3 != "Microglia2"; bad[[1]]@meta <- bad[[1]]@meta[keep, , drop = FALSE]; bad[[1]]@idents <- bad[[1]]@idents[keep]
  errors(calculate_s9e(main, bad), "eligible subtype")
})
test("annotation mismatch rejected before S9E calculation", {
  bad <- chats; bad[[1]]@meta$level3[1] <- "Microglia1"; errors(calculate_s9e(main, bad), "labels disagree")
})
test("constant and variable S9E display normalization", {
  equals(s9e_scale_01(c(5, 5)), c(0, 0)); equals(s9e_scale_01(c(2, 3, 4)), c(0, .5, 1)); errors(s9e_scale_01(c(1, NA)))
})
test("real Brain GAM k4 has finite predictions and bands", {
  time <- rep(c(0, .13, .4, 1), each = 3); fit <- fit_tdeg_gene("g", 1 + 2*time + sin(4*time) + rep(c(-.05, 0, .05), 4), time, 4)
  stopifnot(fit$success, nrow(fit$predictions) == 100, fit$range > 0, all(fit$predictions$lower95 <= fit$predictions$upper95))
})
test("real Skull GAM k3 works", {
  time <- rep(c(0, .33, 1), each = 3); fit <- fit_tdeg_gene("g", time^2 + rep(c(-.01, 0, .01), 3), time, 3)
  stopifnot(fit$success, all(is.finite(fit$predictions$se)))
})
test("GAM insufficient or degenerate data diagnosed", {
  stopifnot(!fit_tdeg_gene("g", 1:3, c(0, .5, 1), 3)$success, !fit_tdeg_gene("g", 1:6, rep(0, 6), 3)$success)
  errors(fit_tdeg_gene("g", 1:3, 1:4, 3)); errors(fit_tdeg_gene("g", 1:5, 1:5, 2))
})
test("zero excluded trajectory cells makes valid empty CSV", {
  ids <- character(); tab <- data.frame(cell_id = ids, reason = rep("Not assigned to Slingshot lineage1", length(ids)))
  write_csv(tab, file.path(output, "empty_exclusions.csv")); equals(nrow(read_csv(file.path(output, "empty_exclusions.csv"))), 0L)
})
test("PDF devices close on success and exception", {
  before <- grDevices::dev.cur(); with_pdf_device(file.path(output, "test.pdf"), 3, 3, graphics::plot(1:3)); equals(grDevices::dev.cur(), before)
  errors(with_pdf_device(file.path(output, "failure.pdf"), 3, 3, stop("intentional")), "intentional"); equals(grDevices::dev.cur(), before)
})
real_need <- need; need <- function(packages) invisible(TRUE) # Only mock object methods below.
test("legacy macrophage12 deleted and13 renamed idempotently", {
  z <- clean_subtype(object(c("Macrophage1", "Macrophage12", "Macrophage13")), "Macrophage", "level3", "legacy")
  equals(z$object$level3, c("Macrophage1", "Macrophage12")); equals(z$removed, "c2")
  again <- clean_subtype(z$object, "Macrophage", "level3", "final"); equals(again$object$level3, z$object$level3)
})
test("current final macrophage12 retained", equals(clean_subtype(object("Macrophage12"), "Macrophage", "level3", "final")$object$level3, "Macrophage12"))
test("missing and unapproved subtype labels rejected", {
  errors(clean_subtype(object(c("Neutrophil1", "Neutrophil13")), "Neutrophil", "level3", "final"), "unexpected")
  errors(clean_subtype(object(c("Neutrophil1", NA)), "Neutrophil", "level3", "final"), "missing labels")
  errors(clean_subtype(object(c("Neutrophil1", "")), "Neutrophil", "level3", "final"))
})
test("approved Neu10-12 exclusions remain explicit", {
  z <- clean_subtype(object(c("Neutrophil1", paste0("Neutrophil", 10:12))), "Neutrophil", "level3", "legacy")
  equals(z$object$level3, "Neutrophil1"); equals(length(z$removed), 3L)
})
test("final annotation selects explicit main input and Skull-only erythroid", {
  old_root <- TBI_ROOT; old_cfg <- cfg
  TBI_ROOT <- file.path(output, "annotation_fixture"); ensure_dir(project_path("config")); ensure_dir(project_path("inputs", "objects")); ensure_dir(project_path("results", "objects"))
  tryCatch({
    labels <- unlist(subtype_orders, use.names = FALSE)
    parents <- unlist(lapply(names(subtype_orders), function(name) rep(lineage_parent[[name]][1], length(subtype_orders[[name]]))), use.names = FALSE)
    fixture <- object(c(labels, "Erythrocyte"))
    fixture$level2 <- c(parents, "Erythrocyte"); fixture$Group <- c(ifelse(parents == "Erythrocyte", "Skull", "Brain"), "Brain")
    fixture$Day <- "D0"; fixture$Lib <- paste0("TBI_", fixture$Group)
    fixture@misc$tag <- "configured"; saveRDS(fixture, input_path(cfg$input_main))
    fixture@misc$tag <- "generated"; saveRDS(fixture, project_path("results", "objects", "integrated_level2.rds"))
    write_csv(data.frame(cell_id = character(), reason = character()), project_path("config", "manual_excluded_cells.csv"))
    cfg$input_mode <- "final_main"; cfg$main_source <- "configured"
    run_annotation(); equals(readRDS(main_final_path())@misc$tag, "configured")
    ery <- readRDS(project_path("results", "objects", "Erythroid_final.rds")); stopifnot(all(ery$Group == "Skull"), ncol(ery) == 3L)
    final <- readRDS(main_final_path())
    stopifnot(setequal(ery$level3, subtype_orders$Erythroid), all(final$level3[final$level2=="Erythrocyte"]=="Erythrocyte"),
      all(ery$global_level3=="Erythrocyte"), all(!c("Erythroblast1","Erythroblast2") %in% final$level3))
    equals(ery@misc$tbi_main_md5, hash_input(main_final_path()))
    counts <- read_csv(project_path("results", "audit", "S2C_S3C_S5B_library_n.csv"))
    stopifnot(setequal(counts$panel, c("S2C", "S3C", "S5B")), all(counts$Group[counts$panel %in% c("S2C", "S5B")] == "Brain"))
    equals(counts$n_libraries[counts$panel == "S3C" & counts$Group == "Skull"], 0L)
    equals(counts$n_libraries[counts$panel == "S3C" & counts$Group == "Brain"], 1L)
    cfg$main_source <- "generated"; run_annotation(); equals(readRDS(main_final_path())@misc$tag, "generated")
  }, finally = { TBI_ROOT <- old_root; cfg <- old_cfg })
})
test("stale subtype fingerprints and Day/Lib/Group mismatches rejected", {
  old_root <- TBI_ROOT; old_main <- get_main
  TBI_ROOT <- file.path(output, "fingerprint_fixture"); ensure_dir(project_path("results", "objects"))
  tryCatch({
    fixture <- object("Microglia1"); fixture$Lib <- "TBI_A"; fixture$Day <- "D0"; fixture$Group <- "Brain"; fixture$level2 <- "Microglia"
    saveRDS(fixture, main_final_path()); get_main <- function() fixture
    path <- project_path("results", "objects", "Microglia_final.rds")
    saveRDS(fixture, path); errors(get_lineage("Microglia"), "fingerprint")
    for (field in c("Day", "Lib", "Group")) {
      bad <- fixture; bad@misc$tbi_main_md5 <- hash_input(main_final_path()); bad@meta.data[[field]] <- "wrong"; saveRDS(bad, path)
      errors(get_lineage("Microglia"), paste0(field, " metadata disagree"))
    }
  }, finally = { TBI_ROOT <- old_root; get_main <- old_main })
})
test("metadata rows reorder by cell ID and library maps are checked", {
  x <- object(c("Microglia1", "Microglia2")); x$Lib <- c("TBI_A", "TBI_B"); x$Day <- c("D0", "D1"); x$Group <- "Brain"; x$level2 <- "Microglia"
  x@meta.data <- x@meta.data[2:1, , drop = FALSE]; y <- normalize_metadata(x)
  equals(y$Lib, c("TBI_A", "TBI_B")); equals(y$Mice_per_library, c(4L, 4L))
  y$Lib <- "TBI_A"; errors(normalize_metadata(y), "multiple")
})
test("CellChat cache rejects changed settings and missing study objects", {
  old_root <- TBI_ROOT
  TBI_ROOT <- file.path(output, "cellchat_cache_fixture"); ensure_dir(project_path("results", "objects")); ensure_dir(project_path("results", "cellchat"))
  tryCatch({
    saveRDS(1, main_final_path())
    cc <- setNames(rep(list(1), 7), c(paste0("Brain_", day_order), paste0("Skull_", day_order[1:3])))
    attr(cc, "downsampled") <- FALSE; attr(cc, "main_md5") <- hash_input(main_final_path())
    attr(cc, "settings") <- cfg[c("cellchat_min_cells", "cellchat_nboot", "cellchat_p", "cellchat_seed")]
    saveRDS(cc, cellchat_path()); equals(length(get_cellchat()), 7L)
    bad <- cc; attr(bad, "settings")$cellchat_min_cells <- 51L; saveRDS(bad, cellchat_path()); errors(get_cellchat(), "settings changed")
    bad <- cc; names(bad)[7] <- "Skull_M6"; saveRDS(bad, cellchat_path()); errors(get_cellchat(), "seven study objects")
    bad <- cc; attr(bad, "main_md5") <- NULL; saveRDS(bad, cellchat_path()); errors(get_cellchat(), "stale")
  }, finally = { TBI_ROOT <- old_root })
})

test("derived communication cache cannot bypass CellChat settings validation", {
  old_root <- TBI_ROOT; TBI_ROOT <- file.path(output, "derived_cache_fixture")
  for (folder in c("results/objects", "results/cellchat", "functions")) ensure_dir(project_path(folder))
  tryCatch({
    saveRDS(1, main_final_path())
    for (file in c("cellchat.R", "common.R")) writeLines("fixture", project_path("functions", file))
    cc <- setNames(rep(list(1), 7), c(paste0("Brain_", day_order), paste0("Skull_", day_order[1:3])))
    attr(cc, "downsampled") <- FALSE; attr(cc, "main_md5") <- hash_input(main_final_path())
    attr(cc, "settings") <- cfg[c("cellchat_min_cells", "cellchat_nboot", "cellchat_p", "cellchat_seed")]
    make_cache <- function() {
      cached <- data.frame(cache_sentinel = 17L)
      attr(cached, "input_stamp") <- list(main_md5 = hash_input(main_final_path()), cellchat_md5 = hash_input(cellchat_path()), p = cfg$cellchat_p,
        classifier_code_md5 = hash_input(project_path("functions", "cellchat.R")), label_code_md5 = hash_input(project_path("functions", "common.R")))
      saveRDS(cached, project_path("results", "cellchat", "significant_communications_nine_categories.rds"))
    }
    bad <- cc; attr(bad, "settings")$cellchat_seed <- 1L; saveRDS(bad, cellchat_path()); make_cache()
    errors(get_communication_table(), "settings changed")
    saveRDS(cc, cellchat_path()); make_cache(); equals(get_communication_table()$cache_sentinel, 17L)
  }, finally = { TBI_ROOT <- old_root })
})

test("main expression consumers use RNA even when imported assay is integrated", {
  old_root <- TBI_ROOT; TBI_ROOT <- file.path(output, "main_assay_fixture"); ensure_dir(project_path("results", "objects"))
  tryCatch({
    x <- object(c("Microglia1", "Microglia2")); x@misc$tbi_schema <- "TBI_final_level3_v1"; x@misc$active_assay <- "integrated"
    saveRDS(x, main_final_path())
    # Substitute only the unavailable Seurat assay setter; execute the real getter body.
    replace_setter <- function(expr) {
      if (missing(expr)) return(expr)
      if (identical(expr, quote(Seurat::DefaultAssay))) return(as.name("fixture_assay"))
      if (is.call(expr)) for (i in seq_along(expr)) expr[[i]] <- replace_setter(expr[[i]])
      expr
    }
    `fixture_assay<-` <- function(object, value) { object@misc$active_assay <- value; object }
    getter <- get_main; body(getter) <- replace_setter(body(getter)); environment(getter) <- environment()
    equals(getter()@misc$active_assay, "RNA")
  }, finally = { TBI_ROOT <- old_root })
})
test("tDEG cache rejects changed subtype input and settings", {
  old_root <- TBI_ROOT
  TBI_ROOT <- file.path(output, "tdeg_cache_fixture"); ensure_dir(project_path("results", "objects")); ensure_dir(project_path("results", "tDEG", "Microglia_Brain"))
  tryCatch({
    saveRDS(1, main_final_path()); path <- project_path("results", "objects", "Microglia_final.rds"); saveRDS(2, path)
    value <- list(key = "Microglia_Brain", main_md5 = hash_input(main_final_path()), lineage_md5 = hash_input(path),
      analysis_settings = cfg[c("tdeg_de_p", "tdeg_de_abs_fc", "tdeg_min_pct", "tdeg_logfc_threshold", "range_quantile", "prediction_points", "seed", "module_score_seed", "module_score_nbin", "module_score_ctrl")])
    result <- project_path("results", "tDEG", "Microglia_Brain", "tDEG_final.rds"); saveRDS(value, result)
    equals(read_tdeg_result("Microglia", "Brain")$key, "Microglia_Brain")
    bad <- value; bad$analysis_settings$range_quantile <- .5; saveRDS(bad, result); errors(read_tdeg_result("Microglia", "Brain"), "settings changed")
    saveRDS(value, result); saveRDS(3, path); errors(read_tdeg_result("Microglia", "Brain"), "subtype fingerprint")
  }, finally = { TBI_ROOT <- old_root })
})
need <- real_need
test("empty marker output has usable columns and malformed nonempty output fails", {
  empty <- normalize_marker_result(data.frame())
  equals(nrow(empty), 0L); stopifnot(all(c("gene", "cluster", "p_val_adj", "avg_log2FC") %in% names(empty)))
  equals(nrow(empty[!is.na(empty$p_val_adj) & empty$p_val_adj < .05 & empty$avg_log2FC > .25, ]), 0L)
  equals(length(split(empty$gene, empty$cluster)), 0L)
  errors(normalize_marker_result(data.frame(gene = "Spp1")), "lacks columns")
})
test("production t-tests use Student equal variance and pooled degrees of freedom", {
  calls <- list()
  walk <- function(expr) {
    if (missing(expr) || !is.call(expr)) return(invisible(NULL))
    if (identical(expr[[1]], quote(stats::t.test)) || identical(expr[[1]], quote(t.test))) calls[[length(calls) + 1L]] <<- expr
    for (i in seq_along(expr)) walk(expr[[i]])
  }
  walk(body(proportion_analysis)); stopifnot(length(calls) == 1L)
  x <- c(1, 2, 3, 4); y <- c(2, 6, 10, 14, 18)
  for (call in calls) {
    stopifnot(identical(call$var.equal, TRUE))
    fit <- eval(call, envir = list(x = x, y = y))
    df <- length(x) + length(y) - 2L
    pooled <- ((length(x)-1)*stats::var(x) + (length(y)-1)*stats::var(y))/df
    statistic <- (mean(x)-mean(y))/sqrt(pooled*(1/length(x)+1/length(y)))
    equals(unname(fit$parameter), df); equals(fit$p.value, 2*stats::pt(-abs(statistic), df))
    stopifnot(!isTRUE(all.equal(fit$p.value, stats::t.test(x, y, var.equal = FALSE)$p.value)))
  }
})
test("stage failure and success statuses recorded", {
  old_root <- TBI_ROOT; TBI_ROOT <- file.path(output, "fixture_project")
  tryCatch({
    errors(with_stage("unit_failure", stop("intentional")), "intentional")
    equals(readRDS(project_path("results", "audit", "unit_failure.status.rds"))$status, "failed")
    with_stage("unit_success", 1 + 1)
    equals(readRDS(project_path("results", "audit", "unit_success.status.rds"))$status, "success")
    with_stage("03", 1 + 1)
    equals(readRDS(project_path("results", "audit", "3.status.rds"))$status, "success")
    errors(with_stage("03", stop("numbered failure")), "numbered failure")
    equals(readRDS(project_path("results", "audit", "3.status.rds"))$status, "failed")
  }, finally = { TBI_ROOT <- old_root })
})
test("all configured random seeds equal 1234 and PCA overrides package seed", {
  stopifnot(identical(cfg$seed, 1234L), identical(cfg$module_score_seed, 1234L), identical(cfg$cellchat_seed, 1234L))
  calls <- list()
  walk <- function(expr) {
    if (missing(expr) || !is.call(expr)) return(invisible(NULL))
    if (identical(expr[[1]], quote(Seurat::RunPCA))) calls[[length(calls)+1L]] <<- expr
    for (i in seq_along(expr)) walk(expr[[i]])
  }
  for (fun in list(run_qc, run_integration, run_subtype_reclustering_review)) walk(body(fun))
  stopifnot(length(calls) == 3L)
  for (call in calls) stopifnot(identical(call$seed.use, quote(cfg$seed)))
})
test("module overlap highlights curated pairs and uses final module genes only", {
  brain <- data.frame(gene = c("a", "b", "c", "d"), cluster = c("C1", "C1", "C2", "C2"))
  skull <- data.frame(gene = c("a", "c", "e"), cluster = c("C1", "C1", "C2"))
  pairs <- data.frame(brain = c("C1", "C2"), skull = c("C1", "C2"))
  x <- build_module_overlap(brain, skull, pairs)
  equals(x$summary$n_shared, c(1L, 0L)); equals(x$summary$fraction_brain, c(.5, 0))
  equals(x$summary$jaccard, c(1/3, 0)); equals(x$shared$gene, "a")
  equals(x$membership$intersect_type[match(c("a", "c", "e"), x$membership$gene)], c("Corresponding", "Other", "Other"))
  equals(x$sets, c("Brain_C1", "Brain_C2", "Skull_C1", "Skull_C2"))
  bad <- pairs; bad$skull[1] <- "C9"; errors(build_module_overlap(brain, skull, bad), "Invalid corresponding")
  errors(build_module_overlap(rbind(brain, brain[1, ]), skull, pairs), "one module")
})
test("all four lineages have source-supported trend module correspondence", {
  pairs <- read_csv(project_path("config", "module_overlap_pairs.csv"))
  stopifnot(nrow(pairs) == 24L, !anyDuplicated(pairs[, c("lineage", "brain", "skull")]))
  equals(sort(unique(pairs$lineage)), sort(c("Macrophage", "Neutrophil", "Bcell", "Tcell")))
  for (i in seq_len(nrow(pairs))) {
    stopifnot(pairs$brain[i] %in% cfg$tdeg_module_order[[paste0(pairs$lineage[i], "_Brain")]],
              pairs$skull[i] %in% cfg$tdeg_module_order[[paste0(pairs$lineage[i], "_Skull")]])
  }
})
test("atlas aliases and ten display classes preserve analytical input metadata", {
  m <- data.frame(Source = c("Public", "In-house", "our", "ref"), Group = c("Brain", "Skull", "Brain", "Brain"),
    Day = c("D0", "D1", "D7", "M6"), level2 = c("Macrophage", "B Cell", "Erythrocyte", "Astrocyte"),
    level1 = c("Myeloid", "Lymphoid", "Non-immune", "Non-immune"))
  before <- m; p <- prepare_atlas_metadata(m)
  equals(p$Source, c("ref", "our", "our", "ref")); equals(p$level1, c("Myeloid Cell", "Lymphoid Cell", "Erythrocyte", "Astrocyte"))
  equals(m, before)
  bad <- m; bad$Source[1] <- "unknown"; errors(prepare_atlas_metadata(bad), "Unrecognized Source")
  bad <- m; bad$Group[1] <- "Skull"; errors(prepare_atlas_metadata(bad), "Published atlas source")
  stopifnot(identical(names(atlas_level1_colors()), c("Myeloid Cell", "Lymphoid Cell", "Erythrocyte", "Astrocyte", "Endothelial Cell", "Pericyte", "Oligodendrocyte", "Neuron", "CPEC", "Ependymal Cell")))
})
test("Figure1D pooling weights cell counts rather than averaging library fractions", {
  # Library A: X=9/10; B: X=0/2. Pooled X is 9/12, not the library mean 45%.
  m <- data.frame(Group = "Brain", Day = "D0", level2 = c(rep("X", 9), rep("Y", 3)))
  p <- atlas_pooled_proportions(m, c("Group", "Day"), "level2")
  equals(p$prop[p$level2 == "X"], .75); equals(sum(p$prop), 1)
  errors(atlas_pooled_proportions(m[FALSE, ], c("Group", "Day"), "level2"), "empty")
})
test("S1H denominator includes all Brain cells and fills absent types with zero", {
  m <- data.frame(Lib = c(rep("TBI_A", 10), rep("TBI_B", 2)), Day = "D0", Group = "Brain", Source = "our", Mice_per_library = 4L,
    level2 = c(rep("Macrophage", 9), "Erythrocyte", rep("Erythrocyte", 2)))
  p <- atlas_library_level2_proportions(m, "Macrophage")
  equals(p$prop[match(c("TBI_A", "TBI_B"), p$Lib)], c(90, 0)); equals(sort(p$denominator), c(2L, 10L))
})
test("S1H uses three adjacent raw Student tests with explicit library counts", {
  p <- data.frame(Lib = paste0("L", seq_len(12)), Day = rep(day_order, each = 3), level2 = "Macrophage", prop = c(1, 2, 4, 3, 7, 9, 2, 5, 8, 4, 10, 11))
  x <- atlas_adjacent_tests(p, "Macrophage")
  equals(paste(x$from, x$to), c("D0 D1", "D1 D7", "D7 M6")); equals(x$n1, rep(3L, 3)); equals(x$n2, rep(3L, 3))
  equals(x$p[1], stats::t.test(p$prop[4:6], p$prop[1:3], alternative = "two.sided", var.equal = TRUE)$p.value)
  stopifnot(!"p_BH" %in% names(x), all(x$status == "tested"))
  errors(atlas_adjacent_tests(rbind(p, p[1, ]), "Macrophage"), "Duplicate")
  p <- p[c(1, 4, 7, 10), ]; x <- atlas_adjacent_tests(p, "Macrophage")
  stopifnot(all(is.na(x$p)), all(x$status == "insufficient_libraries"))
})


test("failed or running upstream stages block previous CellChat and tDEG caches", {
  old_root <- TBI_ROOT; old_need <- need
  TBI_ROOT <- file.path(output, "failed_stage_fixture")
  for (folder in c("results/objects", "results/cellchat", "results/audit", "results/tDEG/Microglia_Brain")) ensure_dir(project_path(folder))
  need <- function(...) invisible(NULL)
  tryCatch({
    saveRDS(1, main_final_path()); saveRDS(2, project_path("results/objects/Microglia_final.rds"))
    cc <- setNames(rep(list(1), 7), c(paste0("Brain_", day_order), paste0("Skull_", day_order[1:3])))
    attr(cc, "downsampled") <- FALSE; attr(cc, "main_md5") <- hash_input(main_final_path())
    attr(cc, "settings") <- cfg[c("cellchat_min_cells", "cellchat_nboot", "cellchat_p", "cellchat_seed")]
    saveRDS(cc, cellchat_path())
    value <- list(key = "Microglia_Brain", main_md5 = hash_input(main_final_path()), lineage_md5 = hash_input(project_path("results/objects/Microglia_final.rds")),
      analysis_settings = cfg[c("tdeg_de_p", "tdeg_de_abs_fc", "tdeg_min_pct", "tdeg_logfc_threshold", "range_quantile", "prediction_points", "seed", "module_score_seed", "module_score_nbin", "module_score_ctrl")])
    path <- project_path("results/tDEG/Microglia_Brain/tDEG_final.rds"); saveRDS(value, path)
    status <- function(stage, value) saveRDS(list(status = value), project_path("results/audit", paste0(stage, ".status.rds")))
    for (state in c("failed", "running")) {
      status(3, "success"); status(4, state); status(15, state)
      errors(get_cellchat(), "Script 4"); errors(read_tdeg_result("Microglia", "Brain"), "Script 15")
      status(4, "success"); status(15, "success"); status(3, state)
      errors(get_cellchat(), "Script 3"); errors(read_tdeg_result("Microglia", "Brain"), "Script 3")
    }
    status(3, "success"); equals(length(get_cellchat()), 7L)
    equals(read_tdeg_result("Microglia", "Brain")$key, "Microglia_Brain")
    value$key <- "Neutrophil_Brain"; saveRDS(value, path)
    errors(read_tdeg_result("Microglia", "Brain"), "identity")
  }, finally = { TBI_ROOT <- old_root; need <- old_need })
})

test("final_main rejects unannotated Brain NVU cells while allowing Brain erythrocytes", {
  old_root <- TBI_ROOT; old_cfg <- cfg; old_need <- need
  TBI_ROOT <- file.path(output, "unannotated_nvu_fixture")
  for (folder in c("config", "inputs/objects", "results/objects")) ensure_dir(project_path(folder))
  need <- function(...) invisible(NULL)
  tryCatch({
    labels <- unlist(subtype_orders, use.names = FALSE)
    parents <- unlist(lapply(names(subtype_orders), function(name) rep(lineage_parent[[name]][1], length(subtype_orders[[name]]))), use.names = FALSE)
    fixture <- object(c(labels, "Erythrocyte", "Astrocyte"))
    fixture$level2 <- c(parents, "Erythrocyte", "Astrocyte")
    fixture$Group <- c(ifelse(parents == "Erythrocyte", "Skull", "Brain"), "Brain", "Brain")
    fixture$Day <- "D0"; fixture$Lib <- paste0("TBI_", fixture$Group)
    write_csv(data.frame(cell_id = character(), reason = character()), project_path("config/manual_excluded_cells.csv"))
    cfg$input_mode <- "final_main"; cfg$main_source <- "configured"
    saveRDS(fixture, input_path(cfg$input_main)); errors(run_annotation(), "Unannotated cells remain in NVU")
    fixture <- subset(fixture, cells = head(colnames(fixture), -1)); saveRDS(fixture, input_path(cfg$input_main))
    run_annotation(); stopifnot(any(readRDS(main_final_path())$level3 == "Erythrocyte"))
  }, finally = { TBI_ROOT <- old_root; cfg <- old_cfg; need <- old_need })
})

test("trajectory input fingerprints reject changed objects, seeds and CSV despite identical cell IDs", {
  old_root <- TBI_ROOT; old_cfg <- cfg
  TBI_ROOT <- file.path(output, "trajectory_cache_fixture")
  for (folder in c("functions", "results/objects", "results/audit")) ensure_dir(project_path(folder))
  tryCatch({
    writeLines("fixture", project_path("functions/trajectory.R"))
    saveRDS(1, main_final_path()); saveRDS(2, project_path("results/objects/Erythroid_final.rds"))
    obj <- object(c("Erythroblast1", "Erythroblast2", "Erythrocyte"))
    d <- outdir("Figure6_S12", "Erythroid_trajectory", "All")
    path <- file.path(d, "slingshot_pseudotime.csv"); stamp_path <- file.path(d, "slingshot_input_stamp.rds")
    tab <- data.frame(cell_id = colnames(obj), slingPseudotime_1 = c(0, 1, 2)); write_csv(tab, path)
    errors(read_slingshot_pseudotime(obj), "provenance is missing")
    save_stamp <- function() saveRDS(list(inputs = trajectory_input_stamp(), csv_md5 = hash_input(path)), stamp_path)
    save_stamp(); equals(read_slingshot_pseudotime(obj)$slingPseudotime_1, c(0, 1, 2))
    saveRDS(3, project_path("results/objects/Erythroid_final.rds")); errors(read_slingshot_pseudotime(obj), "changed")
    save_stamp(); cfg$seed <- 1L; errors(read_slingshot_pseudotime(obj), "changed"); cfg <- old_cfg
    save_stamp(); write_csv(transform(tab, slingPseudotime_1 = c(5, 4, 3)), path); errors(read_slingshot_pseudotime(obj), "changed")
    write_csv(rbind(tab, tab[1, ]), path); save_stamp(); errors(read_slingshot_pseudotime(obj), "duplicated")
    write_csv(tab, path); save_stamp(); saveRDS(list(status = "failed"), project_path("results/audit/27.status.rds"))
    errors(read_slingshot_pseudotime(obj), "Script 27")
    errors(trajectory_input_stamp("M6"), "Invalid")
    # Per-day analysis is allowed inside running script 29 and still validates its own input stamp.
    saveRDS(list(status = "running"), project_path("results/audit/29.status.rds"))
    d1 <- outdir("Figure6_S12", "Erythroid_trajectory", "D1")
    write_csv(tab, file.path(d1, "slingshot_pseudotime.csv"))
    saveRDS(list(inputs = trajectory_input_stamp("D1"), csv_md5 = hash_input(file.path(d1, "slingshot_pseudotime.csv"))), file.path(d1, "slingshot_input_stamp.rds"))
    equals(read_slingshot_pseudotime(obj, "D1")$cell_id, colnames(obj))
  }, finally = { TBI_ROOT <- old_root; cfg <- old_cfg })
})

test("confirmed global 58 labels and supplied counts agree without forcing observed totals", {
  confirmed <- read_csv(project_path("config", "confirmed_level3_counts.csv"))
  stopifnot(nrow(confirmed) == 58L, sum(confirmed$confirmed_cells) == 256134L,
    setequal(confirmed$level3, global_level3_labels()), !any(c("Erythroblast1", "Erythroblast2") %in% global_level3_labels()))
  equals(global_erythroid_labels(c("Erythroblast1", "Erythroblast2", "Erythrocyte", "Microglia1")), c(rep("Erythrocyte",3), "Microglia1"))
})

test("NVU source two-step calculation includes all tissues in its denominator", {
  main <- data.frame(Day = rep(c("D0", "D1"), c(10, 20)), level2 = c(rep("Astrocyte",2),rep("T Cell",8),rep("Pericyte",5),rep("T Cell",15)), Group = "Brain")
  main$Group[c(3,4,16,17)] <- "Skull"
  local <- data.frame(Day = c("D0", "D0", "D1", "D1", "D1"), level3 = c("DAA", "homeostatic_Astro", "SMC", "SMC", "VLMCs"))
  r <- calculate_nvu_all_cell_proportions(main, local)
  equals(r$totals$n_all, c(10L,20L)); equals(r$totals$three_type_prop, c(.2,.25))
  sums <- stats::aggregate(r$composition$plot_height, list(Day=r$composition$Day), sum)
  equals(sums$x[match(c("D0","D1"),sums$Day)],c(.2,.25))
  equals(r$composition$plot_height[r$composition$Day=="D1" & r$composition$level3=="SMC"],(2/3)*.25)
  equals(r$composition$plot_height[r$composition$Day=="D0" & r$composition$level3=="SMC"],0)
  errors(calculate_nvu_all_cell_proportions(main,local[local$Day=="D0",]),"lacks local NVU")
})

test("author gene-subtype choices are used directly and excluded types are not renamed", {
  choices <- read_csv(project_path("config", "manual_gene_subtypes.csv"))
  stopifnot(nrow(choices) == 175L, length(unique(paste(choices$key,choices$panel_id))) == 9L,
    all(choices$selection == "author_manual"))
  mac <- choices[choices$key=="Macrophage_Brain" & choices$gene=="Gpnmb",]
  equals(sort(unique(mac$level3)),c("Macrophage2","Macrophage5"))
  summary <- data.frame(gene="Gpnmb",level3=c("Macrophage1","Macrophage2","Macrophage5"),Day="D1",mean=c(99,2,3))
  filtered <- manual_gene_subset(summary,mac)
  equals(sort(filtered$mean),c(2,3))
  # A higher expression in another subtype must not override the author's choice.
  stopifnot(!"Macrophage1" %in% filtered$level3)
})

test("directional CellChat counts sum to strict category total including zero days", {
  tab <- data.frame(Day=c("D0","D0","D1"),Source_Class=c("Myeloid","Lymphoid","Myeloid"),Target_Class=c("Lymphoid","Myeloid","Lymphoid"),Strength=c(.2,.3,.4))
  r <- communication_direction_metrics(tab,day_order,c("Myeloid","Lymphoid"))
  equals(r$Count[r$Direction=="Total"],c(2L,1L,0L,0L))
  equals(r$Weight[r$Direction=="Total"],c(.5,.4,0,0))
  equals(r$Count[r$Direction=="Myeloid->Lymphoid"],c(1L,1L,0L,0L))
  r <- communication_direction_metrics(tab[FALSE,],day_order,c("Myeloid","Lymphoid"))
  stopifnot(nrow(r)==12L,all(r$Count==0L),all(r$Weight==0))
})

test("sender-receiver source frequency counts entries across days and keeps roles separate", {
  tab <- data.frame(Source=c("A","A","B"),Target=c("B","B","A"),Day=c("D0","D1","D1"))
  r <- communication_role_frequencies(tab)
  equals(r$frequency[r$Cell=="A" & r$Role=="Sender"],2L)
  equals(r$frequency[r$Cell=="A" & r$Role=="Receiver"],1L)
  stopifnot(sum(r$frequency)==2*nrow(tab),nrow(communication_role_frequencies(tab[FALSE,]))==0L)
  by_day <- communication_role_frequencies(tab, by_day = TRUE)
  stopifnot(nrow(by_day)==6L,all(by_day$frequency==1L),sum(by_day$frequency)==2*nrow(tab),"Day" %in% names(communication_role_frequencies(tab[FALSE,],by_day=TRUE)))
})

test("GO filtering is reproducible without a historical selected CSV", {
  tab <- data.frame(group=c("A","A","A","B"),Description=c("first","second","reject","other"),p.adjust=c(.01,.02,.05,.03),Count=c(2,3,100,5),ONTOLOGY="BP")
  selected <- select_go_display_terms(tab,top_n=1L)
  equals(selected$Description,c("first","other"))
  equals(select_go_display_terms(tab,selected_descriptions="second")$Description,"second")
  stopifnot(nrow(select_go_display_terms(tab[FALSE,]))==0L)
})

test("all module-pair gene exports preserve manually declared trend correspondence", {
  b <- data.frame(gene=c("g1","g2","g3"),cluster=c("C1","C1","C2"))
  s <- data.frame(gene=c("g1","g2","g4"),cluster=c("C1","C2","C2"))
  p <- data.frame(brain="C1",skull="C1")
  r <- all_module_pair_gene_overlaps(b,s,p)
  equals(sort(r$gene),c("g1","g2")); stopifnot(r$is_corresponding[r$gene=="g1"],!r$is_corresponding[r$gene=="g2"])
})

test("final_main with global erythrocytes imports local trajectory states without changing global counts", {
  old_root<-TBI_ROOT;old_cfg<-cfg;old_need<-need
  TBI_ROOT<-file.path(output,"global_local_ery_fixture");need<-function(...) invisible(NULL)
  for(folder in c("config","inputs/objects","results/objects")) ensure_dir(project_path(folder))
  tryCatch({
    labels<-unlist(subtype_orders,use.names=FALSE)
    parents<-unlist(lapply(names(subtype_orders),function(n) rep(lineage_parent[[n]][1],length(subtype_orders[[n]]))),use.names=FALSE)
    fixture<-object(labels);fixture$level2<-parents;fixture$Group<-ifelse(parents=="Erythrocyte","Skull","Brain")
    fixture$Day<-"D0";fixture$Lib<-paste0("TBI_",fixture$Group)
    local<-subset(fixture,cells=colnames(fixture)[fixture$level2=="Erythrocyte"])
    local$subtype<-local$level3;saveRDS(local,project_path("inputs/objects/Erythroid.rds"))
    fixture$level3<-global_erythroid_labels(fixture$level3);saveRDS(fixture,input_path(cfg$input_main))
    write_csv(data.frame(lineage="Erythroid",path="inputs/objects/Erythroid.rds",label_column="subtype",annotation_state="final"),input_path(cfg$subtype_manifest))
    write_csv(data.frame(cell_id=character(),reason=character()),project_path("config/manual_excluded_cells.csv"))
    cfg$input_mode<-"final_main";cfg$main_source<-"configured";run_annotation()
    final<-readRDS(main_final_path());ery<-readRDS(project_path("results/objects/Erythroid_final.rds"))
    equals(ncol(final),ncol(fixture));stopifnot(all(final$level3[final$level2=="Erythrocyte"]=="Erythrocyte"),setequal(ery$level3,subtype_orders$Erythroid))
  },finally={TBI_ROOT<-old_root;cfg<-old_cfg;need<-old_need})
})

test("erythroid getter validates global metadata while retaining local subtype labels", {
  old_root<-TBI_ROOT;old_main<-get_main;old_need<-need
  TBI_ROOT<-file.path(output,"ery_getter_fixture");ensure_dir(project_path("results/objects"));need<-function(...) invisible(NULL)
  tryCatch({
    local<-object(subtype_orders$Erythroid);local$Lib<-"TBI_S";local$Day<-"D0";local$Group<-"Skull";local$level2<-"Erythrocyte"
    fixture<-local;fixture$level3<-"Erythrocyte";saveRDS(fixture,main_final_path())
    local@misc$tbi_main_md5<-hash_input(main_final_path());saveRDS(local,project_path("results/objects/Erythroid_final.rds"));get_main<-function() fixture
    replace_assay<-function(expr){if(missing(expr))return(expr);if(identical(expr,quote(Seurat::DefaultAssay)))return(as.name("fixture_assay"));if(is.call(expr))for(i in seq_along(expr))expr[[i]]<-replace_assay(expr[[i]]);expr}
    `fixture_assay<-`<-function(object,value)object
    getter<-get_lineage;body(getter)<-replace_assay(body(getter));environment(getter)<-environment()
    equals(getter("Erythroid")$level3,subtype_orders$Erythroid)
    fixture$level3[1]<-"Microglia1";errors(getter("Erythroid"),"level3 metadata disagree")
  },finally={TBI_ROOT<-old_root;get_main<-old_main;need<-old_need})
})

test("CellChat manual core filters preserve either-cell and bidirectional cross-set rules", {
  tab<-data.frame(Source=c("A","B","A","X"),Target=c("X","A","B","Y"))
  choice<-data.frame(rule="either",group="A",level3="A")
  equals(nrow(filter_author_cellchat_core(tab,choice)),3L)
  choice<-data.frame(rule="cross",group=c("A","B"),level3=c("A","B"))
  selected<-filter_author_cellchat_core(tab,choice)
  equals(paste(selected$Source,selected$Target),c("B A","A B"))
  config<-read_csv(project_path("config","paper_cellchat_core_cells.csv"))
  stopifnot(nrow(config)==51L,all(config$level3 %in% global_level3_labels()))
})

tab <- do.call(rbind, results); write_csv(tab, file.path(output_base, "regression_results.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_base, "test_sessionInfo.txt"))
cat(sum(tab$pass), "/", nrow(tab), "tests passed\n")
if (!all(tab$pass)) quit(status = 1L)
