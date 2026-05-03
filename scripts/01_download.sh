#!/usr/bin/env bash
# =============================================================================
# 01_download.sh
# Download Caribbean dengue whole-genome sequences from NCBI GenBank using
# Entrez Direct (EDirect).
#
# Requirements:
#   - NCBI EDirect installed (https://www.ncbi.nlm.nih.gov/books/NBK179288/)
#   - Optional: set NCBI_API_KEY in your environment to raise rate limits
#
# Output:
#   data/raw/dengue_caribbean_genomes.fasta   — raw multi-FASTA
#   data/raw/retrieval_metadata.txt           — accession + submission info
# =============================================================================

set -euo pipefail

OUTDIR="data/raw"
mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Search query: Caribbean countries, complete genomes, all four serotypes
# ---------------------------------------------------------------------------
QUERY='dengue virus[organism] AND ("complete genome"[title] OR "complete sequence"[title])
       AND (Jamaica[All Fields] OR "Dominican Republic"[All Fields]
            OR "Puerto Rico"[All Fields] OR Trinidad[All Fields]
            OR Barbados[All Fields] OR Haiti[All Fields]
            OR Cuba[All Fields] OR Martinique[All Fields]
            OR Guadeloupe[All Fields] OR Bahamas[All Fields])'

echo "[$(date +%T)] Searching NCBI Nucleotide..."
esearch -db nucleotide -query "$QUERY" \
  | efetch -format acc \
  > "$OUTDIR/accessions.txt"

ACC_COUNT=$(wc -l < "$OUTDIR/accessions.txt")
echo "[$(date +%T)] Found $ACC_COUNT accessions."

# ---------------------------------------------------------------------------
# Download FASTA sequences (batched to respect NCBI rate limits)
# ---------------------------------------------------------------------------
echo "[$(date +%T)] Downloading FASTA sequences..."
esearch -db nucleotide -query "$QUERY" \
  | efetch -format fasta \
  > "$OUTDIR/dengue_caribbean_genomes.fasta"

# ---------------------------------------------------------------------------
# Download GenBank records for metadata extraction
# ---------------------------------------------------------------------------
echo "[$(date +%T)] Downloading GenBank records for metadata..."
esearch -db nucleotide -query "$QUERY" \
  | efetch -format gb \
  > "$OUTDIR/dengue_caribbean_genomes.gb"

# ---------------------------------------------------------------------------
# Write retrieval metadata (timestamp, query, counts)
# ---------------------------------------------------------------------------
cat > "$OUTDIR/retrieval_metadata.txt" <<EOF
Retrieval date : $(date -u +"%Y-%m-%dT%H:%M:%SZ")
Database       : NCBI Nucleotide
Query          : $QUERY
Accessions     : $ACC_COUNT
FASTA output   : $OUTDIR/dengue_caribbean_genomes.fasta
GenBank output : $OUTDIR/dengue_caribbean_genomes.gb
EDirect version: $(edirect --version 2>/dev/null || echo "unknown")
EOF

echo "[$(date +%T)] Done. Metadata written to $OUTDIR/retrieval_metadata.txt"
