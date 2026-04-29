"""
calculate_metrics.py
─────────────────────
Count variants from bcftools isec output files and compute:
  - sensitivity  = TP / (TP + FN)   [how much of truth does the query find?]
  - precision    = TP / (TP + FP)   [how much of the query is in truth?]
  - F1           = harmonic mean of sensitivity and precision

bcftools isec file mapping
  0000.vcf  only in truth  → FN (missed by query)
  0001.vcf  only in query  → FP (not in truth)
  0002.vcf  shared (truth side) → TP
"""

import csv
import json
import os
import subprocess
import sys


BCFTOOLS    = snakemake.params.bcftools
TRUTH_LABEL = snakemake.params.truth_label
QUERY_LABEL = snakemake.params.query_label


def count_vcf_variants(vcf_path: str) -> int:
    """Return number of non-header variant lines in a VCF (handles empty files)."""
    result = subprocess.run(
        [BCFTOOLS, "view", "-H", vcf_path],
        capture_output=True,
        text=True,
        check=True,
    )
    lines = [l for l in result.stdout.splitlines() if l.strip()]
    return len(lines)


# ── Snakemake interface ───────────────────────────────────────────────────────

log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path), exist_ok=True)

with open(log_path, "w") as log:
    try:
        fn = count_vcf_variants(snakemake.input.only_a)    # only in truth
        fp = count_vcf_variants(snakemake.input.only_b)    # only in query
        tp = count_vcf_variants(snakemake.input.shared_a)  # shared

        total_truth = tp + fn   # total in truth pipeline
        total_query = tp + fp   # total in query pipeline

        sensitivity = tp / total_truth if total_truth > 0 else 0.0
        precision   = tp / total_query if total_query > 0 else 0.0
        f1 = (
            2 * precision * sensitivity / (precision + sensitivity)
            if (precision + sensitivity) > 0
            else 0.0
        )

        metrics = {
            "pipeline": str(snakemake.wildcards.pipeline),
            "variant_type": "SNV",
            "truth_pipeline": TRUTH_LABEL,
            "query_pipeline": QUERY_LABEL,
            "counts": {
                f"only_{TRUTH_LABEL}": fn,
                f"only_{QUERY_LABEL}": fp,
                "shared":              tp,
                f"total_{TRUTH_LABEL}": total_truth,
                f"total_{QUERY_LABEL}": total_query,
            },
            "metrics": {
                "sensitivity": round(sensitivity, 6),
                "precision":   round(precision,   6),
                "f1":          round(f1,           6),
            },
        }

        log.write(f"Pipeline:         {snakemake.wildcards.pipeline}\n")
        log.write(f"Truth:            {TRUTH_LABEL}\n")
        log.write(f"Query:            {QUERY_LABEL}\n")
        log.write(f"FN (only truth):  {fn:,}\n")
        log.write(f"FP (only query):  {fp:,}\n")
        log.write(f"TP (shared):      {tp:,}\n")
        log.write(f"Total RD:         {total_truth:,}\n")
        log.write(f"Total VP:         {total_query:,}\n")
        log.write(f"Sensitivity:      {sensitivity:.4%}\n")
        log.write(f"Precision:        {precision:.4%}\n")
        log.write(f"F1 score:         {f1:.6f}\n")

        os.makedirs(os.path.dirname(snakemake.output.metrics), exist_ok=True)
        with open(snakemake.output.metrics, "w") as fh:
            json.dump(metrics, fh, indent=2)

        # Write flat CSV for downstream use / cohort aggregation
        with open(snakemake.output.csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "sample", "variant_type", "truth_pipeline", "query_pipeline",
                "fn", "fp", "tp", "total_truth", "total_query",
                "sensitivity", "precision", "f1",
            ])
            writer.writerow([
                metrics["pipeline"], metrics["variant_type"],
                TRUTH_LABEL, QUERY_LABEL,
                fn, fp, tp, total_truth, total_query,
                round(sensitivity, 6), round(precision, 6), round(f1, 6),
            ])

        print(json.dumps(metrics, indent=2))

    except Exception as exc:
        log.write(f"ERROR: {exc}\n")
        raise
