#!/usr/bin/env bash
# =============================================================================
# 02_align.sh
# Align Caribbean dengue sequences to serotype-specific reference genomes.
#
# Strategy:
#   1. Split the multi-FASTA by serotype (based on header keywords / BLAST)
#   2. Align each serotype set to its reference using minimap2 (fast) or
#      MAFFT --auto (for smaller curated sets)
#   3. Trim alignment ends and export trimmed FASTA for downstream analyses
#
# Requirements:
#   - minimap2   (https://github.com/lh3/minimap2)
#   - samtools
#   - MAFFT      (https://mafft.cbrc.jp/alignment/software/)
#   - seqkit     (optional, for splitting by header)
#
# Reference accessions (NCBI):
#   DENV-1: NC_001477   DENV-2: NC_001474
#   DENV-3: NC_001475   DENV-4: NC_002640
# =============================================================================

set -euo pipefail

RAW_DIR="data/raw"
ALIGN_DIR="data/aligned"
REF_DIR="data/references"

mkdir -p "$ALIGN_DIR" "$REF_DIR"

SEROTYPES=(1 2 3 4)
declare -A REFS=(
  [1]="NC_001477"
  [2]="NC_001474"
  [3]="NC_001475"
  [4]="NC_002640"
)

# ---------------------------------------------------------------------------
# Step 1: Download reference genomes if not present
# ---------------------------------------------------------------------------
for ST in "${SEROTYPES[@]}"; do
  ACC="${REFS[$ST]}"
  REF_FILE="$REF_DIR/DENV${ST}_${ACC}.fasta"
  if [[ ! -f "$REF_FILE" ]]; then
    echo "[$(date +%T)] Fetching reference DENV-$ST ($ACC)..."
    efetch -db nucleotide -id "$ACC" -format fasta > "$REF_FILE"
  else
    echo "[$(date +%T)] Reference DENV-$ST already present."
  fi
done

# ---------------------------------------------------------------------------
# Step 2: Classify sequences by serotype using minimap2 vs all four refs
# ---------------------------------------------------------------------------
INPUT_FASTA="$RAW_DIR/dengue_caribbean_genomes.fasta"

echo "[$(date +%T)] Classifying sequences by serotype..."
for ST in "${SEROTYPES[@]}"; do
  ACC="${REFS[$ST]}"
  REF_FILE="$REF_DIR/DENV${ST}_${ACC}.fasta"

  # Map all sequences; keep only those that align (primary, non-supplementary)
  minimap2 -a --secondary=no "$REF_FILE" "$INPUT_FASTA" 2>/dev/null \
    | samtools view -F 4 -F 2048 \
    | awk '{print $1}' \
    | sort -u \
    > "$ALIGN_DIR/DENV${ST}_ids.txt"

  COUNT=$(wc -l < "$ALIGN_DIR/DENV${ST}_ids.txt")
  echo "  DENV-$ST: $COUNT sequences mapped"
done

# ---------------------------------------------------------------------------
# Step 3: Extract per-serotype FASTAs and align with MAFFT
# ---------------------------------------------------------------------------
for ST in "${SEROTYPES[@]}"; do
  ACC="${REFS[$ST]}"
  REF_FILE="$REF_DIR/DENV${ST}_${ACC}.fasta"
  IDS_FILE="$ALIGN_DIR/DENV${ST}_ids.txt"
  SUBSET_FASTA="$ALIGN_DIR/DENV${ST}_subset.fasta"
  ALIGNED_FASTA="$ALIGN_DIR/DENV${ST}_aligned.fasta"

  echo "[$(date +%T)] Aligning DENV-$ST sequences..."

  # Extract sequences matching the ID list
  seqkit grep -f "$IDS_FILE" "$INPUT_FASTA" > "$SUBSET_FASTA"

  # Prepend reference so it anchors the alignment
  cat "$REF_FILE" "$SUBSET_FASTA" > "${SUBSET_FASTA%.fasta}_with_ref.fasta"

  # MAFFT alignment (--auto selects algorithm based on dataset size)
  mafft --auto --thread -1 \
    "${SUBSET_FASTA%.fasta}_with_ref.fasta" \
    > "$ALIGNED_FASTA"

  echo "[$(date +%T)]  → $ALIGNED_FASTA"
done

echo "[$(date +%T)] Alignment complete. Files in $ALIGN_DIR/"
