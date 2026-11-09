# Required local inputs

Input data are not included. Place files at the paths below or edit `config/settings.R` and `config/subtype_inputs.csv`. Paths are relative to the repository unless absolute. Keep gene symbols and globally unique cell IDs consistent across objects.

## Main object

Default: `inputs/objects/main_integrated.rds` (alternatively `integrated_level2.rds`). A Seurat object with an RNA assay and normalized `data`; counts are needed for analyses that use raw counts. Metadata row names must match cell IDs. Required fields: `Lib` (or `sample_batch_library`), `Day` (D0/D1/D7/M6), `Group` (Brain/Skull), and valid `level2` labels. `Source`, library pooling and default control labels are normalized by the study-specific code. TBI-prefixed libraries are In-house and pool four mice; other libraries are Public. Each Lib must map to exactly one Day and Group. M6 is Public Brain only.

## Annotated lineage objects

Seven paths and label columns are specified in `config/subtype_inputs.csv`: Macrophage, Neutrophil, Microglia, Bcell, Tcell, Erythroid and NVU. Cell IDs must match the main object. Objects must retain the RNA assay and required embeddings; the local erythroid object must include its UMAP. Erythroid is Skull-only with Erythroblast1/Erythroblast2/Erythrocyte states; NVU is Brain-only. Do not replace local erythroid labels with global Erythrocyte labels.

`annotation_state=final` expects final subtype labels; Macrophage defaults to `auto`, which distinguishes supported legacy/final schemes. Check `subtype_cluster_maps.csv` and explicit exclusions when using raw cluster labels. The scripts reject ambiguous or unapproved mappings instead of silently inventing labels.

## Reference validation (entry 6)

Default reference: `inputs/reference/GSE192616_annotated.rds`, annotated metadata column `levelx`. It must contain the 16 matched reference labels expected by `run_reference()` and the embedding used for display. Optional historical marker file: `inputs/reference/skull_markers.csv`, columns `gene`, `cluster`, `avg_log2FC`. Without that CSV, markers are recomputed and the exact historical gene set is not guaranteed.

## Matrices (optional entries 1–2)

Complete all 31 rows of `config/sample_manifest.csv` using actual metadata, including currently blank Day values. For each library, place exactly one matrix.mtx(.gz), barcodes.tsv(.gz), and features.tsv(.gz) or genes.tsv(.gz) under the configured matrix directory; GEO filename prefixes are supported. `gene_column` selects the feature column. Existing annotated-object analysis does not require these matrices or the completed raw-QC manifest.

Real input files and analysis outputs are ignored by `.gitignore`; share data through the study's actual data repository and document its accession separately.
