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
##   results/tables/phylo_cluster_summary.csv
##   results/tables/phylo_cluster_memberships.csv
##   results/tables/phylo_cluster_overview.csv
##   results/figures/phylo_cluster_overview.png
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
})

theme_set(theme_bw(base_size = 12))

dir.create("results/trees", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

ALIGN_DIR <- "data/aligned"
TREE_DIR <- "results/trees"
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

normalize_metadata <- function(meta) {
  if (!"accession" %in% names(meta) && "assembly_accession" %in% names(meta)) {
    meta <- meta |>
      mutate(accession = assembly_accession)
  }

  if (!"year" %in% names(meta) && "collection_year" %in% names(meta)) {
    meta <- meta |>
      mutate(year = suppressWarnings(as.integer(collection_year)))
  }

  meta
}

get_descendant_tips <- function(tree, node) {
  child_nodes <- tree$edge[tree$edge[, 1] == node, 2]
  tips <- character(0)

  for (child in child_nodes) {
    if (child <= Ntip(tree)) {
      tips <- c(tips, tree$tip.label[child])
    } else {
      tips <- c(tips, get_descendant_tips(tree, child))
    }
  }

  unique(tips)
}

summarise_phylo_clusters <- function(tree, metadata, serotype, caribbean_countries) {
  tip_meta <- tibble(label = tree$tip.label) |>
    left_join(metadata, by = c("label" = "accession")) |>
    mutate(
      country = coalesce(country, "Unknown"),
      year = suppressWarnings(as.integer(year))
    )

  jamaica_labels <- tip_meta |>
    filter(country == "Jamaica") |>
    pull(label)

  empty_result <- list(
    cluster_summary = tibble(),
    cluster_memberships = tibble(),
    tip_meta = tip_meta
  )

  if (length(jamaica_labels) == 0) {
    return(empty_result)
  }

  internal_nodes <- seq.int(Ntip(tree) + 1L, Ntip(tree) + tree$Nnode)
  node_descendants <- setNames(vector("list", length(internal_nodes)), as.character(internal_nodes))
  node_has_jamaica <- setNames(logical(length(internal_nodes)), as.character(internal_nodes))

  for (node in internal_nodes) {
    descendant_tips <- get_descendant_tips(tree, node)
    node_descendants[[as.character(node)]] <- descendant_tips
    node_has_jamaica[[as.character(node)]] <- any(descendant_tips %in% jamaica_labels)
  }

  minimal_nodes <- integer(0)
  for (node in internal_nodes) {
    child_nodes <- tree$edge[tree$edge[, 1] == node, 2]
    child_nodes <- child_nodes[child_nodes > Ntip(tree)]
    child_has_jamaica <- length(child_nodes) > 0 && any(vapply(as.character(child_nodes), function(child_node) isTRUE(node_has_jamaica[[child_node]]), logical(1)))

    if (isTRUE(node_has_jamaica[[as.character(node)]]) && !child_has_jamaica) {
      minimal_nodes <- c(minimal_nodes, node)
    }
  }

  if (length(minimal_nodes) == 0) {
    return(empty_result)
  }

  summary_rows <- vector("list", length(minimal_nodes))
  member_rows <- vector("list", length(minimal_nodes))

  for (i in seq_along(minimal_nodes)) {
    node <- minimal_nodes[[i]]
    cluster_id <- sprintf("DENV%d_C%02d", serotype, i)
    descendant_tips <- node_descendants[[as.character(node)]]

    tip_info <- tip_meta |>
      filter(label %in% descendant_tips) |>
      mutate(cluster_id = cluster_id)

    jamaica_info <- tip_info |>
      filter(country == "Jamaica")

    other_countries <- sort(unique(setdiff(na.omit(tip_info$country), "Jamaica")))
    cluster_type <- case_when(
      nrow(jamaica_info) >= 2 && length(other_countries) == 0 ~ "jamaica_outbreak_cluster",
      nrow(jamaica_info) >= 1 && length(other_countries) > 0 ~ "introduction_cluster",
      nrow(jamaica_info) == 1 && length(other_countries) == 0 ~ "singleton_jamaica_tip",
      TRUE ~ "mixed_cluster"
    )

    node_support <- NA_real_
    node_index <- node - Ntip(tree)
    if (!is.null(tree$node.label) && node_index >= 1 && node_index <= length(tree$node.label)) {
      node_support <- suppressWarnings(as.numeric(tree$node.label[node_index]))
    }

    summary_rows[[i]] <- tibble(
      serotype = paste0("DENV-", serotype),
      cluster_id = cluster_id,
      node = node,
      node_support = node_support,
      n_tips = nrow(tip_info),
      n_jamaica = nrow(jamaica_info),
      n_caribbean = sum(tip_info$country %in% caribbean_countries, na.rm = TRUE),
      cluster_type = cluster_type,
      countries = paste(sort(unique(na.omit(tip_info$country))), collapse = "; "),
      introduction_countries = paste(other_countries, collapse = "; "),
      jamaica_accessions = paste(sort(unique(jamaica_info$label)), collapse = "; "),
      jamaica_year_min = if (nrow(jamaica_info) > 0) suppressWarnings(min(jamaica_info$year, na.rm = TRUE)) else NA_real_,
      jamaica_year_max = if (nrow(jamaica_info) > 0) suppressWarnings(max(jamaica_info$year, na.rm = TRUE)) else NA_real_
    )

    member_rows[[i]] <- tip_info |>
      transmute(
        serotype = paste0("DENV-", serotype),
        cluster_id = cluster_id,
        cluster_type = cluster_type,
        accession = label,
        country = country,
        year = year,
        is_jamaica = country == "Jamaica"
      )
  }

  list(
    cluster_summary = bind_rows(summary_rows),
    cluster_memberships = bind_rows(member_rows),
    tip_meta = tip_meta
  )
}

plot_dengue_tree <- function(treefile, metadata, serotype, caribbean_countries) {
  tree <- read.tree(treefile)
  tip_meta <- tibble(label = tree$tip.label) |>
    left_join(metadata, by = c("label" = "accession")) |>
    mutate(
      country = coalesce(country, "Unknown"),
      region = case_when(
        country == "Jamaica" ~ "Jamaica",
        country %in% caribbean_countries ~ "Caribbean",
        is.na(country) ~ "Unknown",
        TRUE ~ "Other"
      )
    )

  tip_colours <- c(
    Jamaica = "#E63946",
    Caribbean = "#457B9D",
    Other = "#AAAAAA",
    Unknown = "#CCCCCC"
  )[tip_meta$region]

  outfile <- sprintf("results/figures/DENV%d_phylo.png", serotype)
  png(outfile, width = 1800, height = 1600, res = 180)
  par(mar = c(2, 1, 4, 10), xpd = NA)
  plot(
    tree,
    show.tip.label = FALSE,
    cex = 0.55,
    no.margin = TRUE,
    main = sprintf("DENV-%d maximum-likelihood phylogeny", serotype),
    sub = "Tips coloured by geographic origin"
  )
  tiplabels(pch = 19, col = tip_colours, cex = 0.7)
  legend(
    "topright",
    inset = c(-0.22, 0),
    legend = c("Jamaica", "Caribbean", "Other", "Unknown"),
    pch = 19,
    col = c("#E63946", "#457B9D", "#AAAAAA", "#CCCCCC"),
    bty = "n",
    xpd = NA,
    title = "Origin"
  )
  axisPhylo()
  dev.off()
}

# ---------------------------------------------------------------------------
# 1. Load metadata for tip annotation
# ---------------------------------------------------------------------------
meta <- read_csv("data/metadata/dengue_global_genomes.csv", show_col_types = FALSE) |>
  normalize_metadata()

if (!"accession" %in% names(meta)) {
  stop("Metadata must contain accession values before phylogenetic analysis")
}

caribbean_countries <- c(
  "Jamaica", "Dominican Republic", "Puerto Rico", "Trinidad and Tobago",
  "Barbados", "Haiti", "Cuba", "Martinique", "Guadeloupe", "Bahamas"
)

# ---------------------------------------------------------------------------
# 2. Build ML tree per serotype
# ---------------------------------------------------------------------------
tree_files <- list()
cluster_summaries <- list()
cluster_memberships <- list()

for (ST in SEROTYPES) {
  aln_file <- file.path(ALIGN_DIR, sprintf("DENV%d_aligned.fasta", ST))
  out_prefix <- file.path(TREE_DIR, sprintf("DENV%d_ml", ST))

  if (!file.exists(aln_file)) {
    warning("Alignment not found, skipping DENV-", ST, ": ", aln_file)
    next
  }

  treefile <- run_iqtree(aln_file, out_prefix)
  tree_files[[as.character(ST)]] <- treefile
  message("DENV-", ST, " tree saved: ", treefile)

  tree <- read.tree(treefile)
  cluster_info <- summarise_phylo_clusters(
    tree = tree,
    metadata = meta,
    serotype = ST,
    caribbean_countries = caribbean_countries
  )

  if (nrow(cluster_info$cluster_summary) > 0) {
    cluster_summaries[[length(cluster_summaries) + 1]] <- cluster_info$cluster_summary
    cluster_memberships[[length(cluster_memberships) + 1]] <- cluster_info$cluster_memberships
  }

  plot_dengue_tree(
    treefile = treefile,
    metadata = meta,
    serotype = ST,
    caribbean_countries = caribbean_countries
  )
}

cluster_summary <- bind_rows(cluster_summaries)
cluster_membership_table <- bind_rows(cluster_memberships)

if (nrow(cluster_summary) > 0) {
  write_csv(cluster_summary, "results/tables/phylo_cluster_summary.csv")
  write_csv(cluster_membership_table, "results/tables/phylo_cluster_memberships.csv")

  cluster_overview <- cluster_summary |>
    count(serotype, cluster_type, name = "n_clusters") |>
    arrange(serotype, desc(n_clusters))

  write_csv(cluster_overview, "results/tables/phylo_cluster_overview.csv")

  cluster_plot <- cluster_overview |>
    ggplot(aes(x = serotype, y = n_clusters, fill = cluster_type)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Jamaica introduction and outbreak cluster overview",
      x = "Serotype",
      y = "Cluster count",
      fill = "Cluster type"
    )

  ggsave(
    "results/figures/phylo_cluster_overview.png",
    cluster_plot,
    width = 10,
    height = 5,
    dpi = 300
  )

  message("Cluster summary written to results/tables/phylo_cluster_summary.csv")
}

message("\n--- Cluster overview ---")
if (nrow(cluster_summary) > 0) {
  print(
    cluster_summary |>
      group_by(serotype) |>
      summarise(
        n_clusters = n(),
        n_outbreak_clusters = sum(cluster_type == "jamaica_outbreak_cluster"),
        n_introduction_clusters = sum(cluster_type == "introduction_cluster"),
        n_singleton_jamaica = sum(cluster_type == "singleton_jamaica_tip"),
        .groups = "drop"
      )
  )
} else {
  message("No Jamaica-containing clusters were detected in the available trees.")
}

# ---------------------------------------------------------------------------
# 4. Session info
# ---------------------------------------------------------------------------
message("\n--- Session info ---")
sessionInfo()
