#!/bin/bash

#SBATCH --partition=shortterm
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH --mem=32GB
#SBATCH --job-name=vcf-bench
#SBATCH --output=logs/slurm_%j_%u_%N.out
#SBATCH --error=logs/slurm_%j_%u_%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alihassan1697@gmail.com

# ============================================================
# Variant Calling Benchmark — SLURM submission script
# ============================================================
#
# Usage:
#   sbatch submit_job.sh [configfile] [OPTIONS]
#
# Options:
#   --pipeline NAME   Key name for CLI-supplied files — used as output folder
#                     name and in plot titles (default: cli_input)
#   --snv PATH        SNV VCF (.vcf.gz) to benchmark
#   --sv  PATH        SV  VCF (.vcf.gz) to benchmark
#
# Examples:
#   # Use paths defined in config/config.yaml:
#   sbatch submit_job.sh
#   sbatch submit_job.sh config/config_test.yaml
#
#   # Benchmark a specific SNV+SV VCF pair without editing the config:
#   sbatch submit_job.sh config/config.yaml \
#       --pipeline varient_piper_v2 \
#       --snv      /data/calls/HG002.snv.vcf.gz \
#       --sv       /data/calls/HG002.sv.vcf.gz
#
#   # SNV only (omit --sv → no SV benchmarking):
#   sbatch submit_job.sh config/config.yaml \
#       --pipeline varient_piper_v2 \
#       --snv      /data/calls/HG002.snv.vcf.gz
#
# Override SLURM resources at submission time:
#   sbatch --cpus-per-task=16 --mem=64GB --time=1-00:00:00 submit_job.sh config/config.yaml
# ============================================================

set -euo pipefail

# ── Argument parsing ──────────────────────────────────────────────────────────

CONFIGFILE="config/config.yaml"
PIPELINE_NAME=""
SNV_PATH=""
SV_PATH=""

# First positional argument (if not a flag) is the config file
if [[ $# -ge 1 && "${1:0:1}" != "-" ]]; then
    CONFIGFILE="$1"
    shift
fi

# Remaining named flags
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pipeline|-p) PIPELINE_NAME="$2"; shift 2 ;;
        --snv)         SNV_PATH="$2";      shift 2 ;;
        --sv)          SV_PATH="$2";       shift 2 ;;
        *) echo "ERROR: unknown option '$1'" >&2; exit 1 ;;
    esac
done

# ── Environment setup ─────────────────────────────────────────────────────────

module load singularity/v4.1.3

MINIFORGE_PATH="/work/hassan/hassan/miniforge"
source "${MINIFORGE_PATH}/etc/profile.d/conda.sh"
conda activate snakemake

mkdir -p logs

# Tee all output into a timestamped log file alongside the SLURM .out/.err
TIMESTAMP=$(date +%Y%m%dT%H%M%S)
RUN_LOG="logs/run_${TIMESTAMP}_${SLURM_JOB_ID:-interactive}.log"
exec > >(tee -a "$RUN_LOG") 2>&1

# ── Header ────────────────────────────────────────────────────────────────────

echo "=============================================="
echo "  Variant Calling Benchmark (omics cluster)"
echo "=============================================="
echo "Job ID:            ${SLURM_JOB_ID:-interactive}"
echo "Node:              $(hostname)"
echo "Working directory: $(pwd)"
echo "Snakemake version: $(snakemake --version)"
echo "Date:              $(date)"
echo "CPUs:              ${SLURM_CPUS_PER_TASK:-8}"
echo "Memory:            ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "Config:            $CONFIGFILE"
[[ -n "$SNV_PATH" ]] && echo "SNV input:         $SNV_PATH"
[[ -n "$SV_PATH"  ]] && echo "SV  input:         $SV_PATH"
echo "Run log:           $RUN_LOG"
echo "=============================================="

# ── Build CLI override config (if --snv / --sv were given) ───────────────────
#
# When files are supplied on the command line a minimal YAML is written to a
# temp file and passed as a second --configfile.  Snakemake merges configfiles
# left-to-right, so values here override whatever is in config/config.yaml.
#
# The override sets ONLY the pipeline(s) derived from the CLI args; all other
# top-level keys (truth, truvari, reference, …) are inherited from the base
# config.

EXTRA_CONFIGFILE=""

if [[ -n "$SNV_PATH" || -n "$SV_PATH" ]]; then

    # Validate that provided paths exist and are indexed
    for PATH_VAR in "$SNV_PATH" "$SV_PATH"; do
        [[ -z "$PATH_VAR" ]] && continue
        if [[ ! -f "$PATH_VAR" ]]; then
            echo "ERROR: file not found: $PATH_VAR" >&2
            exit 1
        fi
        if [[ ! -f "${PATH_VAR}.tbi" ]]; then
            echo "ERROR: index not found: ${PATH_VAR}.tbi  (run tabix -p vcf first)" >&2
            exit 1
        fi
    done

    [[ -z "$PIPELINE_NAME" ]] && PIPELINE_NAME="cli_input"

    EXTRA_CONFIGFILE=$(mktemp /tmp/bench_override_XXXXXX.yaml)
    trap "rm -f $EXTRA_CONFIGFILE" EXIT

    {
        echo "# Auto-generated override — do not edit"
        echo "pipelines:"
        echo "  ${PIPELINE_NAME}:"
        [[ -n "$SNV_PATH" ]] && echo "    snv: \"${SNV_PATH}\""
        [[ -n "$SV_PATH"  ]] && echo "    sv:  \"${SV_PATH}\""
        echo ""

        # Only register the pipeline for the variant types that were provided
        if [[ -n "$SNV_PATH" ]]; then
            echo "snv_pipelines:"
            echo "  - ${PIPELINE_NAME}"
        fi
        if [[ -n "$SV_PATH" ]]; then
            echo "sv_pipelines:"
            echo "  - ${PIPELINE_NAME}"
        fi
    } > "$EXTRA_CONFIGFILE"

    echo "--- CLI override config ---"
    cat "$EXTRA_CONFIGFILE"
    echo "---------------------------"
fi

# ── Assemble Snakemake configfile arguments ───────────────────────────────────

CONFIGFILE_ARGS=(--configfile "$CONFIGFILE")
[[ -n "$EXTRA_CONFIGFILE" ]] && CONFIGFILE_ARGS+=(--configfile "$EXTRA_CONFIGFILE")

# ── Unlock in case a previous job was killed or timed out ────────────────────

mkdir -p .snakemake/locks
snakemake "${CONFIGFILE_ARGS[@]}" --unlock 2>/dev/null || true

# ── Run ───────────────────────────────────────────────────────────────────────

snakemake \
    "${CONFIGFILE_ARGS[@]}" \
    --use-conda \
    --conda-prefix /work/hassan/hassan/snakemake-conda \
    --cores "${SLURM_CPUS_PER_TASK:-8}" \
    --rerun-incomplete \
    --nolock

EXIT_CODE=$?

# ── Summary ───────────────────────────────────────────────────────────────────

echo "=============================================="
if [[ "$EXIT_CODE" -eq 0 ]]; then
    echo "Pipeline finished successfully"
else
    echo "Pipeline FAILED with exit code $EXIT_CODE"
    echo ""
    echo "--- Last 50 lines of Snakemake internal log ---"
    SMLOG=$(ls -t .snakemake/log/*.snakemake.log 2>/dev/null | head -1)
    if [[ -n "$SMLOG" ]]; then
        echo "Snakemake log: $SMLOG"
        tail -50 "$SMLOG"
    fi
fi
echo "=============================================="
echo "Full run log: $RUN_LOG"
echo "=============================================="

exit "$EXIT_CODE"
