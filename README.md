# Dengue Caribbean Genomics

> **Tracking dengue virus evolution across the Caribbean — with a focus on Jamaica — using publicly available whole-genome sequences.**

---

## Why This Project Exists

Dengue fever is endemic across the Caribbean. Jamaica has experienced recurrent outbreaks driven by all four dengue serotypes (DENV-1–4), yet the country's genomic surveillance capacity remains limited. Without knowing *which lineages* are circulating, *where they came from*, and *how they relate* to regional and global strains, public-health responses are reactive rather than anticipatory.

This project closes that gap by:

1. Downloading all available Caribbean dengue whole-genome sequences from NCBI GenBank.
2. Aligning them to reference genomes for each serotype.
3. Identifying serotype composition and the dominant lineages for Jamaica and its neighbours.
4. Reconstructing a maximum-likelihood phylogeny to trace introduction routes and outbreak clusters.

All scripts are fully reproducible, and raw data are retrieved programmatically rather than stored in the repository.

---

## What We Found (In Progress)

| Serotype | Caribbean sequences | Jamaican sequences | Dominant clade |
| -------- | ------------------- | ------------------ | -------------- |
| DENV-1   |                     |                    | Genotype V     |
| DENV-2   |                     |                    | Asian/American |
| DENV-3   |                     |                    | Genotype III   |
| DENV-4   |                     |                    | Genotype II    |

Key findings:
- 

---

## What This Means for Jamaica



---

## Repository Structure

```
dengue-caribbean-genomics/
├── README.md                         ← project narrative (this file)
├── data/
│   ├── raw/                          ← retrieval logs and programmatically generated FASTA outputs
│   └── metadata/
│       └── dengue_global_genomes.csv ← cleaned metadata for all sequences used
├── scripts/
│   ├── 01_download.py                ← fetches sequences from NCBI Entrez
│   ├── 02_align.sh                   ← reference-based alignment with minimap2/MAFFT
│   ├── 03_serotype_analysis.R        ← serotype counts, lineage assignment, figures
│   └── 04_phylogenetics.R            ← IQ-TREE ML phylogeny, annotated tree plots
├── results/
│   └── figures/                      ← publication-ready PNG/SVG outputs
└── report/
    └── dengue_caribbean_report.Rmd   ← full reproducible report (knit to HTML/PDF)
```

---

## How to Reproduce

```bash
# 1. Run the full pipeline from the repository root
bash scripts/run_pipeline.sh

# Or run stages individually:
python scripts/01_download.py
bash scripts/02_align.sh
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" scripts/03_serotype_analysis.R
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" scripts/04_phylogenetics.R
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" -e "rmarkdown::render('report/dengue_caribbean_report.Rmd')"
```

Dependencies are listed in `report/dengue_caribbean_report.Rmd` session-info block and can be installed via `renv::restore()` if an `renv.lock` is present.

---

## Data Sources

- **NCBI GenBank** — complete and near-complete dengue genome sequences, accession metadata in `data/raw/retrieval_metadata.txt`
- **ViPR / BV-BRC** — supplemental Caribbean sequences not yet in GenBank
- **PAHO surveillance reports** — epidemiological context for outbreak years

---

## Citation / Acknowledgements

If you use this analysis, please cite the NCBI accessions listed in `data/metadata/dengue_global_genomes.csv` and acknowledge the original submitting laboratories. This work is part of ongoing Caribbean infectious-disease genomics research.

---

*Last updated: May 2026*
