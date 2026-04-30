#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
#  install_truth_set.sh
#
#  Downloads the HG002 GIAB truth sets into the locations expected by the
#  pipeline config (config/config.yaml):
#
#    SNV       → input_data/truth/HG002/SNV/        (NISTv4.2.1, GRCh38)
#    SV CMRG   → input_data/truth/HG002/SV/         (CMRG v1.00, GRCh38)
#                273 medically relevant genes, ~250 SVs
#    SV Tier1  → input_data/truth/HG002/SV_Tier1/   (NIST SV v0.6, GRCh38)
#                Genome-wide, ~10k SVs — better coverage for short-read callers
#
#  Usage:
#    bash install_truth_set.sh                    # SNV + SV CMRG (default)
#    bash install_truth_set.sh --snv              # SNV only
#    bash install_truth_set.sh --sv               # SV CMRG only
#    bash install_truth_set.sh --sv-tier1         # SV Tier 1 only
#    bash install_truth_set.sh --all              # SNV + SV CMRG + SV Tier 1
#    bash install_truth_set.sh --snv --sv-tier1   # combine flags freely
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Source URLs ───────────────────────────────────────────────────────────────
SNV_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
SV_CMRG_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant"
SV_TIER1_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716"

# ── Target directories (must match config/config.yaml) ───────────────────────
SNV_DIR="input_data/truth/HG002/SNV"
SV_DIR="input_data/truth/HG002/SV"
SV_TIER1_DIR="input_data/truth/HG002/SV_Tier1"

# ── Parse flags ───────────────────────────────────────────────────────────────
DO_SNV=false
DO_SV=false
DO_SV_TIER1=false

if [[ $# -eq 0 ]]; then
    DO_SNV=true
    DO_SV=true
else
    for flag in "$@"; do
        case "$flag" in
            --snv)      DO_SNV=true ;;
            --sv)       DO_SV=true ;;
            --sv-tier1) DO_SV_TIER1=true ;;
            --all)      DO_SNV=true; DO_SV=true; DO_SV_TIER1=true ;;
            *) echo "Usage: $0 [--snv] [--sv] [--sv-tier1] [--all]"; exit 1 ;;
        esac
    done
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

# ── SV truth set: CMRG (medically relevant genes, ~250 SVs) ──────────────────
if $DO_SV; then
    echo ""
    echo "=== SV truth set: CMRG v1.00 GRCh38 (273 medically relevant genes) ==="
    mkdir -p "$SV_DIR"
    download "${SV_CMRG_BASE}/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz" \
             "${SV_DIR}/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz"
    download "${SV_CMRG_BASE}/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi" \
             "${SV_DIR}/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi"
    download "${SV_CMRG_BASE}/HG002_GRCh38_CMRG_SV_v1.00.bed" \
             "${SV_DIR}/HG002_GRCh38_CMRG_SV_v1.00.bed"
    echo "  SV CMRG files ready in ${SV_DIR}/"
fi

# ── SV truth set: T2T Q100 Draft Benchmark (genome-wide, 2024-07-18) ──────────
# Single VCF covers both small and structural variants; separate BEDs per type.
# For SV benchmarking use the _stvar BED; for SNV use _smvar BED.
if $DO_SV_TIER1; then
    echo ""
    echo "=== SV truth set: GIAB T2T Q100 Draft Benchmark GRCh38 (2024-07-18) ==="
    mkdir -p "$SV_TIER1_DIR"
    download "${SV_TIER1_BASE}/GRCh38_HG2-T2TQ100-V1.1.vcf.gz" \
             "${SV_TIER1_DIR}/GRCh38_HG2-T2TQ100-V1.1.vcf.gz"
    download "${SV_TIER1_BASE}/GRCh38_HG2-T2TQ100-V1.1.vcf.gz.tbi" \
             "${SV_TIER1_DIR}/GRCh38_HG2-T2TQ100-V1.1.vcf.gz.tbi"
    download "${SV_TIER1_BASE}/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed" \
             "${SV_TIER1_DIR}/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"
    echo "  SV Tier 1 files ready in ${SV_TIER1_DIR}/"
    echo ""
    echo "  To use this benchmark, update config/config.yaml truth.sv block:"
    echo "    vcf: ${SV_TIER1_DIR}/GRCh38_HG2-T2TQ100-V1.1.vcf.gz"
    echo "    bed: ${SV_TIER1_DIR}/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"
fi

echo ""
echo "Done. Run the pipeline with:  snakemake --cores 4"
