## =============================================================================
## 04_phylogenetics.R
## Maximum-likelihood phylogenetic analysis of Caribbean dengue sequences.
##
## Inputs:
##   data/aligned/DENV*_aligned.fasta         – per-serotype alignments
##   data/metadata/dengue_global_genomes.csv  – sequence metadata for tip labels
##
## Requires IQ-TREE2 installed and available in PATH.
##
## Outputs:
##   results/trees/DENV*_ml.treefile
##   results/figures/DENV*_phylo.png
##   results/figures/tanglegram_DENV12.png
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggtree)
  library(treeio)
  library(patchwork)
})

theme_set(theme_tree2())

dir.create("results/trees",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

ALIGN_DIR <- "data/aligned"
TREE_DIR  <- "results/trees"
SEROTYPES <- 1:4

# ---------------------------------------------------------------------------
# Helper: run IQ-TREE2 on a FASTA alignment
# ---------------------------------------------------------------------------
run_iqtree <- function(aln_file, out_prefix, model = "GTR+G", threads = 4) {
  cmd <- sprintf(
    "iqtree2 -s %s -m %s -B 1000 --prefix %s -T %d --quiet",
    aln_file, model, out_prefix, threads
  )
  message("[IQ-TREE2] Running: ", cmd)
  ret <- system(cmd)
  if (ret != 0) stop("IQ-TREE2 failed for ", aln_file)
  paste0(out_prefix, ".treefile")
}

# ---------------------------------------------------------------------------
# 1. Load metadata for tip annotation
# ---------------------------------------------------------------------------
meta <- read_csv("data/metadata/dengue_global_genomes.csv",
                 show_col_types = FALSE) |>
  mutate(
    tip_label = accession,
    is_jamaica = country == "Jamaica"
  )

caribbean_countries <- c(
  "Jamaica", "Dominican Republic", "Puerto Rico", "Trinidad and Tobago",
  "Barbados", "Haiti", "Cuba", "Martinique", "Guadeloupe", "Bahamas"
)

# ---------------------------------------------------------------------------
# 2. Build ML tree per serotype
# ---------------------------------------------------------------------------
tree_files <- list()

for (ST in SEROTYPES) {
  aln_file   <- file.path(ALIGN_DIR, sprintf("DENV%d_aligned.fasta", ST))
  out_prefix <- file.path(TREE_DIR,  sprintf("DENV%d_ml", ST))

  if (!file.exists(aln_file)) {
    warning("Alignment not found, skipping DENV-", ST, ": ", aln_file)
    next
  }

  treefile <- run_iqtree(aln_file, out_prefix)
  tree_files[[ST]] <- treefile
  message("DENV-", ST, " tree saved: ", treefile)
}

# ---------------------------------------------------------------------------
# 3. Plot and annotate each tree
# ---------------------------------------------------------------------------
plot_dengue_tree <- function(treefile, metadata, serotype, caribbean) {
  tree <- read.tree(treefile)

  # Merge tip metadata
  tip_meta <- tibble(label = tree$tip.label) |>
    left_join(metadata, by = c("label" = "accession")) |>
    mutate(
      region = case_when(
        country %in% caribbean ~ "Caribbean",
        is.na(country)         ~ "Unknown",
        TRUE                   ~ "Other"
      ),
      highlight = if_else(country == "Jamaica", "Jamaica", region)
    )

  p <- ggtree(tree, layout = "rectangular", size = 0.25) %<+% tip_meta +
    geom_tippoint(aes(colour = highlight), size = 1.2, alpha = 0.85) +
    scale_colour_manual(
      values = c("Jamaica"   = "#E63946",
                 "Caribbean" = "#457B9D",
                 "Other"     = "#AAAAAA",
                 "Unknown"   = "#CCCCCC"),
      name = "Origin"
    ) +
    geom_treescale(x = 0, y = 1, fontsize = 3) +
    labs(title = sprintf("DENV-%d maximum-likelihood phylogeny", serotype),
         subtitle = "Tips coloured by geographic origin") +
    theme(legend.position = "right",
          plot.title      = element_text(face = "bold"))

  p
}

for (ST in SEROTYPES) {
  if (is.null(tree_files[[ST]])) next

  p <- plot_dengue_tree(
    treefile  = tree_files[[ST]],
    metadata  = meta,
    serotype  = ST,
    caribbean = caribbean_countries
  )

  outfile <- sprintf("results/figures/DENV%d_phylo.png", ST)
  ggsave(outfile, p, width = 12, height = 16, dpi = 300)
  message("Saved: ", outfile)
}

# ---------------------------------------------------------------------------
# 4. Session info
# ---------------------------------------------------------------------------
message("\n--- Session info ---")
sessionInfo()
