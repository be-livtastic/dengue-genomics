#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

resolve_python() {
  if command -v python >/dev/null 2>&1; then
    PYTHON_CMD=(python)
  elif command -v python3 >/dev/null 2>&1; then
    PYTHON_CMD=(python3)
  elif command -v py >/dev/null 2>&1; then
    PYTHON_CMD=(py -3)
  else
    echo "Missing Python launcher: python, python3, or py" >&2
    exit 1
  fi
}

resolve_rscript() {
  if [[ -n "${R_SCRIPT_BIN:-}" ]]; then
    if [[ ! -x "$R_SCRIPT_BIN" ]]; then
      echo "R_SCRIPT_BIN is set but not executable: $R_SCRIPT_BIN" >&2
      exit 1
    fi
    R_CMD=("$R_SCRIPT_BIN")
    return
  fi

  if command -v Rscript >/dev/null 2>&1; then
    R_CMD=(Rscript)
    return
  fi

  if [[ -x "C:/Program Files/R/R-4.5.2/bin/Rscript.exe" ]]; then
    R_CMD=("C:/Program Files/R/R-4.5.2/bin/Rscript.exe")
    return
  fi

  if [[ -x "/usr/bin/Rscript" ]]; then
    R_CMD=(/usr/bin/Rscript)
    return
  fi

  echo "Missing Rscript executable. Set R_SCRIPT_BIN to a valid Rscript path." >&2
  exit 1
}

resolve_python
resolve_rscript

for cmd in efetch minimap2 samtools seqkit mafft iqtree2; do
  require_command "$cmd"
done

mkdir -p data/raw data/metadata data/aligned data/references results/figures results/tables results/trees

echo "[1/5] Downloading sequences"
"${PYTHON_CMD[@]}" scripts/01_download.py

echo "[2/5] Aligning sequences"
bash scripts/02_align.sh

echo "[3/5] Running serotype and lineage analysis"
"${R_CMD[@]}" scripts/03_serotype_analysis.R

echo "[4/5] Running phylogenetic analysis"
"${R_CMD[@]}" scripts/04_phylogenetics.R

echo "[5/5] Rendering report"
"${R_CMD[@]}" -e "rmarkdown::render('report/dengue_caribbean_report.Rmd')"

echo "Pipeline complete."
