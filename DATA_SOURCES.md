# Input data: what is distributed, what is not, and which scripts use it

All analyses start from the DIA-NN outputs of the 60 dried blood spot specimens (jPOST JPST005044 / ProteomeXchange
PXD086615). Specimen identifiers: `WW_` healthy adults, `WS_` heterozygous carriers, `SS_` LoF/LoF, `SM_` LoF/Mis,
`DD_` disease controls.

## Distributed in `data/`

| File | Content | Used by |
|---|---|---|
| `raw_data.xlsx` (sheet `data`) | DIA-NN protein-group quantities after cross-run normalization, non-imputed; 7,474 protein groups × 60 specimens, with accession, gene symbol, peptide and unique-peptide counts; zero = not quantified. Identical to Supplementary Table S8. | every proteomic analysis (01–05, 07) |
| `sample_info.csv` | `SampleID`, `group` (WW / WS / SS / SM / DD), `batch` (Run1 / Run2 / Run3) | every proteomic analysis |
| `DBS-Tohoku.pg_matrix.tsv` | DIA-NN MaxLFQ protein-group matrix (alternative normalization) | `04_SLC25A13_measurement/slc25a13_normalization_sensitivity.R` |

## Not distributed (supply through `CITRIN_EXTERNAL_DIR`)

| File | Content | Why not distributed | Used by |
|---|---|---|---|
| `family_map.csv` | `SampleID`, `family_id` (parent–child and sibling relationships) | family structure of a rare-disease cohort; available from the corresponding author on request | family-grouped cross-validation and family-cluster bootstraps (see README) |
| `DBS-Tohoku.pr_matrix.tsv` | DIA-NN precursor matrix of the 60 specimens | part of the jPOST deposit | precursor-level detection (Fig. 2B), peptide uniqueness and detection floor (Tables S25, S26), precursor evidence of panel proteins (Table S28) |
| `raw/250404_GC_Human_UP000005640_E20644.fasta` | UniProt human reference proteome used for the DIA-NN library | public (UniProt UP000005640, downloaded 2025-04-04); part of the jPOST deposit | `peptide_uniqueness_check.R`, `variant_peptide_check.R`, `detection_floor_check.R` |
| `DBS_nACTB_analysis_data_20250808.csv` | ACTB-normalized matrix exported in August 2025; used only for a cross-check | superseded by `raw_data.xlsx` | `03_sensitivity/two_stage/slc25a13_detection_60.R` (optional cross-check) |
| raw files / mzML of five specimens (WW_6, WS_10, SS_13, SM_5, DD_17) and the DIA-NN report (`.parquet`) | fragment-ion chromatograms | jPOST deposit (raw files); convert with ThermoRawFileParser 2.0.0 | `slc25a13_fragment_xic.py`, `slc25a13_xic_figure_final.py` (Fig. 2C, Fig. S4); `results/run_mapping.csv` maps specimens to raw files |
| `counts_featureCounts_20261103.txt`, `info_20261103.xlsx`, `panel_symbols_from_xlsx.csv` | PBMC RNA-seq gene-level counts (featureCounts, GENCODE v44), sample sheet (genotype, age, sex, platform, batch) and the 792-gene capture panel | consent restrictions on per-specimen sequencing data; processed data on request under the conditions stated in the manuscript | `06_RNAseq/Analysis_resubmission_RNAseq_20261103.R` |

## Regenerable intermediate outputs that are not distributed

Per-fold caches and per-run folders of the nested ensemble (`outer_cache/`), of the single-model multi-seed runs, of the
protein-exclusion refits and of the 100 label permutations are not distributed (several GB; hours to days of compute).
Their aggregated outputs are in the `results/` folder of each analysis, and the launchers regenerate them into
`results_rerun/` or into the folder given by the corresponding environment variable.

## Checksums

MD5 sums of all distributed files are in `MANIFEST.md5`. MD5 sums of the non-distributed inputs used for the reported
results:

| File | MD5 |
|---|---|
| `family_map.csv` | 3bdf0c31760ebdbf0e763078588146a5 |
| `DBS-Tohoku.pr_matrix.tsv` | 0ec56d68b5c819d2fecc44d232a1b9ec |
| `raw/250404_GC_Human_UP000005640_E20644.fasta` | 6e661b5fbd55f905efbdf5e20ed2b746 |
| `DBS_nACTB_analysis_data_20250808.csv` | a602c6e173a095a33dc33ae304e6668d |
| `counts_featureCounts_20261103.txt` | 5171d5dca5eb9d0131d7e157de278e49 |
| `info_20261103.xlsx` | 94d9a31ec4e70fd19d4e83ea8807e00e |
| `panel_symbols_from_xlsx.csv` | 04305289cb2d31d9100f87b67dddf60f |
