# 2026_TBI_SkullBrain_Atlas
Analysis code for the single-cell multi-omics atlas of the skull–brain immune compartments in traumatic brain injury

# TBI Brain-Skull Single-Cell Analysis Scripts

This repository contains the organized R analysis code for the TBI brain–skull study. Scripts are numbered 0–37 in manuscript order. Each numbered `.R` file has a matching R-kernel `.ipynb` with the same code; run either format, not both. The implementation is in `functions/`; the short numbered files are analysis entry points.

## Contents

- `0–37`: 38 analysis entry points, with matching notebooks.
- `functions/`: 10 shared R modules; all are required.
- `config/`: 16 configuration and author-selection files; all are required. These CSVs are analysis inputs, not disposable results.
- `inputs/README.md`: required input objects, metadata and matrix layout.
- `qa/regression_tests.R`: 58 core regression checks, including small real mgcv fits. It does not rerun the study data.
- `script_index_and_dependencies.csv`: figure/table coverage and execution dependencies.
- `Figure1_S1_output_dimensions.csv`, `panel_output_mapping.csv`: output mapping and Figure1/S1 export dimensions.
- `SOFTWARE_VERSIONS.md`: versions actually supplied from the analysis server; this is not a complete environment lock.

## Environment and dependencies

Use an R environment compatible with the analysis server (R 4.2.2). Confirmed server package versions are DoubletFinder 2.0.3, org.Mm.eg.db 3.16.0 and pheatmap 1.0.12. The code uses Seurat APIs including `slot` access and was not validated against arbitrary future package versions. Full historical versions of all packages are not available; no `renv.lock` has been fabricated.

Required packages checked by entry 0: Seurat, DoubletFinder, harmony, clusterProfiler, org.Mm.eg.db, mgcv, CellChat, pheatmap, dplyr, tidyr, ggplot2, patchwork, scales, openxlsx, ComplexUpset, slingshot, monocle (Monocle2), SingleCellExperiment, S4Vectors, Biobase, VGAM, Matrix and future. `scattermore` is optional for some rasterized displays. Packages must be installed separately; these scripts do not install or upgrade them.

Use UTF-8. The `.ipynb` files require an R kernel. When running notebooks outside the project working directory, set `TBI_PROJECT_ROOT` to the repository's absolute path before running a cell. For RStudio, source a numbered `.R` file; `Rscript` entry files determine their project root from their own location. In an existing R session, remove or update any old `TBI_ROOT` variable before switching repositories.

## Run from existing annotated objects (default)

1. Supply the main Seurat object, seven annotated lineage objects and optional reference inputs as described in `inputs/README.md`.
2. Edit `config/settings.R` and `config/subtype_inputs.csv`. Defaults use project-relative paths; absolute paths may be configured locally. Keep `input_mode="subtype_objects"`, `main_source="configured"`, `subtype_mode="annotated"` for the default route.
3. Run entry 0 for dependency and input preflight, then entry 3 to synchronize final annotation, then entry 4 to build the full CellChat objects. Entry 0 reports missing packages; inspect its CSV before continuing.
4. Run the desired figure scripts according to `script_index_and_dependencies.csv`. Entry 27 must precede 28; entry 29 fits each time-point trajectory separately. Entry 36 creates global CellChat summaries/Table S2. Run 37 after the complete figure/table workflow; it expects all stages 3–36 and Tables S1–S11.
5. Outputs are generated under `results/`, including runtime session information, input fingerprints and actual counts. Configuration or upstream input changes can invalidate cached results; rerun the affected stages.

`input_mode="final_main"` additionally requires the local three-state Skull erythroid labels (either in `erythroid_subtype` or the configured local erythroid object). Global 58-type annotation alone is insufficient for the erythroid trajectory.

## Run from matrices (optional)

The distributed `config/sample_manifest.csv` is an incomplete template: its Day values must be filled from actual sample metadata. Do not infer sham/naive time assignments from filenames. Entry 0 can export `results/audit/sample_manifest_from_actual_object.csv` from an existing main object for review.

After completing the manifest and placing matrices, run 1 (QC/DoubletFinder), then 2 (global integration). Entry 2 saves an unannotated object and deliberately stops for marker-based review unless `accept_recomputed_cluster_map=TRUE`. Review the cluster map before approving it. Use `main_source="generated"` to consume the resulting level2 object in entry 3. Rebuilding lineage clusters is an optional review route via `subtype_mode="recluster_review"`; inspect marker results and maps rather than assuming old cluster numbers transfer.

No DNBC4tools or SoupX implementation is included. Public-data download/preprocessing is outside this bundle; integration of the supplied public-data objects/matrices is included. GSE192616 is the configured Skull reference object used in entry 6, not a claim that every public study input comes from that accession.

## Analysis conventions

- Seeds are 1234. Harmony uses `Lib` as the batch variable, defaults for global integration, and `max.iter.harmony=40` for lineage subclustering. Lineage clustering resolutions remain lineage-specific.
- TBI-prefixed libraries each pool four mice. Independent libraries/pools, not individual cells or pooled mice, are the units for proportion t tests.
- Proportion tests use two-sided unpaired Student tests (`var.equal=TRUE`). Standard Brain subtype plots compare against D0 and show raw p; BH-adjusted values are additionally exported. S1H intentionally uses adjacent-time comparisons.
- CellChat uses all cells within eligible subtypes (at least 50 cells), no downsampling, `nboot=100`, and an explicit strict `p<0.05` mask for exported interaction tables and derived panels. Diagnostic exports compare the strict mask with package aggregate counts; exact-threshold interactions can differ.
- tDEG candidates are the union of significant DEGs from all available pairwise time-point comparisons within each lineage/compartment. Thresholds are adjusted p<0.05 and absolute log2FC>0.5, with FindMarkers Wilcoxon, min.pct=0.1 and logfc.threshold=0.25. Candidates enter GAM fitting without an HVG intersection or additional GAM-p filter. Only successful fits with fitted range strictly greater than Q25 are retained.
- tDEG curves are gene-wise min–max scaled. Internal 95% bounds are scaled and clipped to [0,1], then averaged within modules. Their averaged width supplies an effective display SE for approximate 50%/75% display bands; these are not formal module-mean confidence intervals. The 95% band is not plotted.
- Current GO code uses BP enrichment, BH correction and `p.adjust` for plotted significance. Historical manuscript figures may have used raw p; labels and figures must be aligned with the actual rerun. No historical selected-GO CSV is needed.
- Manually selected representative genes/subtypes, CellChat display subsets and 24 trend-module correspondences are preserved in configuration. Global annotation has 58 types; Erythroblast1/2/Erythrocyte are local erythroid trajectory states. Confirmed cell counts are comparison targets, not values forced into new results.
- Entry 6 uses Spearman correlation, the union of top 200 markers per cluster by avg_log2FC, then intersects expression matrices. Supply the original marker CSV for historical comparability; otherwise markers are recomputed. The heatmap is column-scaled while raw coefficients are also exported.

## Verification and data availability

On 2026-11-09, 50 R files parsed, all 38 R/notebook pairs matched, and all 58 core checks passed. Real study matrices/RDS were unavailable locally, so this is code-level verification, not proof of full manuscript reproduction. Real model, figure, count and package-compatibility checks must be run in the analysis environment; entry 37 supplies final checks. The supplied bundle contains no study matrices, Seurat RDS or final CellChat objects.

Before publication, provide the actual data-access information and complete runtime package versions. These details have not been invented. Historical exploratory notebooks, old audit records and machine-specific path crosswalks are retained locally and excluded from this publication folder.

## English release and spatial-filtering provenance

This is the English-language function-based edition. All 38 numbered R entries and their 38 notebooks use English filenames; relative input/output paths are English. Source-notebook metadata in configuration uses English aliases for unbundled historical files, not additional runtime inputs. The bilingual replacement map is kept outside this publication folder for local use. Greek characters in biological labels, such as gamma-delta T cells, are retained intentionally.

Historical UMAP-coordinate filtering is not automatically reconstructed by these function-based entries. Entry 3 can inherit cell-ID filtering from supplied lineage objects and synchronize those cells with the main object. However, `config/manual_excluded_cells.csv` is currently empty, and the optional reclustering route does not reproduce historical coordinate rectangles. Without the actual RDS inputs, prior spatial filtering cannot be verified. Other cell types not represented by the seven configured lineage objects also require checking in the main input or an explicit cell-ID exclusion list.

Do not apply old coordinate thresholds to a newly computed UMAP. Use the original embedding and verified retained/excluded IDs. Additional B-cell level4 and CD8 subtyping branches must not automatically be treated as global level3 exclusions. This translation does not enable new filtering or change analysis parameters.

The author selected the MIT License for the GitHub repository. Retain its existing LICENSE file and verify the copyright holder/year. This folder does not fabricate a replacement copyright statement.

