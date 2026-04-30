#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
#  install_truth_set.sh
#
#  Downloads the HG002 GIAB truth set into the locations expected by the
#  pipeline config (config/config.yaml):
#
#    SNV  → input_data/truth/HG002/SNV/   (NISTv4.2.1, GRCh38)
#    SV   → input_data/truth/HG002/SV/    (CMRG v1.00, GRCh38)
#
#  Usage:
#    bash install_truth_set.sh           # download both SNV and SV
#    bash install_truth_set.sh --snv     # SNV only
#    bash install_truth_set.sh --sv      # SV only
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Source URLs ───────────────────────────────────────────────────────────────
SNV_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
SV_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant"

# ── Target directories (must match config/config.yaml) ───────────────────────
SNV_DIR="input_data/truth/HG002/SNV"
SV_DIR="input_data/truth/HG002/SV"

# ── Parse flags ───────────────────────────────────────────────────────────────
DO_SNV=true
DO_SV=true
if [[ $# -gt 0 ]]; then
    case "$1" in
        --snv) DO_SV=false ;;
        --sv)  DO_SNV=false ;;
        *) echo "Usage: $0 [--snv|--sv]"; exit 1 ;;
    esac
fi

# ── Helper: download a single file, skip if already present ──────────────────
download() {
    local url="$1"
    local dest="$2"
    if [[ -f "$dest" ]]; then
        echo "  [skip] $(basename "$dest") already exists"
        return
    fi
    echo "  [get]  $(basename "$dest")"
    wget --quiet --show-progress --continue -O "$dest" "$url"
}

# ── SNV truth set ─────────────────────────────────────────────────────────────
if $DO_SNV; then
    echo ""
    echo "=== SNV truth set (NISTv4.2.1 GRCh38) ==="
    mkdir -p "$SNV_DIR"
    download "${SNV_BASE}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
             "${SNV_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    download "${SNV_BASE}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
             "${SNV_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
    download "${SNV_BASE}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
             "${SNV_DIR}/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
    echo "  SNV files ready in ${SNV_DIR}/"
fi

# ── SV truth set ──────────────────────────────────────────────────────────────
if $DO_SV; then
    echo ""
    echo "=== SV truth set (CMRG v1.00 GRCh38) ==="
    mkdir -p "$SV_DIR"
    download "${SV_BASE}/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz" \
             "${SV_DIR}/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz"
    download "${SV_BASE}/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi" \
             "${SV_DIR}/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi"
    download "${SV_BASE}/HG002_GRCh38_CMRG_SV_v1.00.bed" \
             "${SV_DIR}/HG002_GRCh38_CMRG_SV_v1.00.bed"
    echo "  SV files ready in ${SV_DIR}/"
fi

echo ""
echo "Done. Run the pipeline with:  snakemake --cores 4"
