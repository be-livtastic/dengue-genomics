## =============================================================================
## 03_serotype_analysis.R
## Serotype composition and lineage assignment for Caribbean dengue sequences.
##
## Inputs:
##   data/metadata/dengue_global_genomes.csv  – cleaned sequence metadata
##   data/aligned/DENV*_aligned.fasta         – per-serotype alignments
##
## Outputs:
##   results/figures/serotype_composition.png
##   results/figures/year_serotype_heatmap.png
##   results/figures/jamaica_timeline.png
##   results/tables/lineage_summary.csv
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

theme_set(theme_bw(base_size = 12))

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load metadata
# ---------------------------------------------------------------------------
meta <- read_csv("data/metadata/dengue_global_genomes.csv",
  show_col_types = FALSE
)

# Expected columns: accession, serotype, country, year, lineage, length_bp
message("Sequences loaded: ", nrow(meta))
message("Countries represented: ", n_distinct(meta$country))

if (!"accession" %in% names(meta) && "assembly_accession" %in% names(meta)) {
  meta <- meta |> mutate(accession = assembly_accession)
}

if (!"year" %in% names(meta) && "collection_year" %in% names(meta)) {
  meta <- meta |> mutate(year = suppressWarnings(as.integer(collection_year)))
}

if (!"lineage" %in% names(meta)) {
  meta <- meta |>
    mutate(
      lineage = case_when(
        serotype == "DENV-1" ~ "Genotype V",
        serotype == "DENV-2" ~ "Asian/American",
        serotype == "DENV-3" ~ "Genotype III",
        serotype == "DENV-4" ~ "Genotype II",
        TRUE ~ NA_character_
      )
    )
}

cluster_memberships_path <- "results/tables/phylo_cluster_memberships.csv"
if (file.exists(cluster_memberships_path)) {
  cluster_memberships <- read_csv(cluster_memberships_path, show_col_types = FALSE) |>
    select(accession, cluster_id, cluster_type)

  meta <- meta |>
    left_join(cluster_memberships, by = "accession") |>
    mutate(
      lineage_detail = if_else(
        !is.na(cluster_id),
        paste(lineage, cluster_id, sep = " | "),
        lineage
      ),
      lineage_source = if_else(
        !is.na(cluster_id),
        paste0("phylogeny:", cluster_type),
        "serotype_fallback"
      )
    )
} else {
  meta <- meta |>
    mutate(
      lineage_detail = lineage,
      lineage_source = "serotype_fallback"
    )
}

# ---------------------------------------------------------------------------
# 2. Caribbean subset flag
# ---------------------------------------------------------------------------
caribbean_countries <- c(
  "Jamaica", "Dominican Republic", "Puerto Rico", "Trinidad and Tobago",
  "Barbados", "Haiti", "Cuba", "Martinique", "Guadeloupe", "Bahamas",
  "Belize", "Honduras", "Cayman Islands", "Antigua and Barbuda", "Grenada"
)

meta <- meta |>
  mutate(
    is_caribbean = country %in% caribbean_countries,
    is_jamaica = country == "Jamaica"
  )

carib <- meta |> filter(is_caribbean)
message("Caribbean sequences: ", nrow(carib))
message("Jamaican sequences: ", sum(carib$is_jamaica))

# ---------------------------------------------------------------------------
# 3. Serotype composition — Caribbean bar chart
# ---------------------------------------------------------------------------
p_serotype <- carib |>
  count(serotype, country) |>
  mutate(country = fct_reorder(country, n, sum)) |>
  ggplot(aes(x = n, y = country, fill = serotype)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "DENV-1" = "#E63946", "DENV-2" = "#457B9D",
      "DENV-3" = "#2A9D8F", "DENV-4" = "#F4A261"
    )
  ) +
  labs(
    title    = "Dengue serotype composition — Caribbean sequences",
    subtitle = "All complete genomes available in NCBI GenBank",
    x        = "Number of sequences",
    y        = NULL,
    fill     = "Serotype"
  )

ggsave("results/figures/serotype_composition.png",
  p_serotype,
  width = 9, height = 6, dpi = 300
)

# ---------------------------------------------------------------------------
# 4. Year × serotype heatmap for Jamaica
# ---------------------------------------------------------------------------
jam_heat <- carib |>
  filter(is_jamaica) |>
  count(year, serotype) |>
  complete(year = full_seq(year, 1), serotype, fill = list(n = 0))

p_heat <- ggplot(jam_heat, aes(x = year, y = serotype, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient(low = "white", high = "#E63946", name = "Sequences") +
  scale_x_continuous(breaks = seq(2000, 2025, 5)) +
  labs(
    title    = "Jamaican dengue sequences by year and serotype",
    subtitle = "Zero cells indicate no sequences deposited (not necessarily absence)",
    x        = "Collection year",
    y        = "Serotype"
  )

ggsave("results/figures/year_serotype_heatmap.png",
  p_heat,
  width = 10, height = 4, dpi = 300
)

# ---------------------------------------------------------------------------
# 5. Jamaica temporal trend — stacked area
# ---------------------------------------------------------------------------
p_timeline <- jam_heat |>
  filter(n > 0) |>
  ggplot(aes(x = year, y = n, fill = serotype)) +
  geom_area(alpha = 0.85) +
  scale_fill_manual(
    values = c(
      "DENV-1" = "#E63946", "DENV-2" = "#457B9D",
      "DENV-3" = "#2A9D8F", "DENV-4" = "#F4A261"
    )
  ) +
  scale_x_continuous(breaks = seq(2000, 2025, 2)) +
  labs(
    title = "Jamaica: dengue genomic sequences over time",
    x     = "Year",
    y     = "Number of sequences",
    fill  = "Serotype"
  )

ggsave("results/figures/jamaica_timeline.png",
  p_timeline,
  width = 10, height = 5, dpi = 300
)

# ---------------------------------------------------------------------------
# 6. Lineage summary table
# ---------------------------------------------------------------------------
lineage_summary <- carib |>
  count(serotype, lineage_detail, sort = TRUE) |>
  group_by(serotype) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  ungroup()

write_csv(lineage_summary, "results/tables/lineage_summary.csv")
message("Lineage summary written to results/tables/lineage_summary.csv")
message("Lineage detail sourced from phylogenetic clusters when available.")

# ---------------------------------------------------------------------------
# 7. Print session info for reproducibility
# ---------------------------------------------------------------------------
message("\n--- Session info ---")
sessionInfo()
