"""
plot_chromosome.py
───────────────────
Two-panel per-chromosome figure:
  Panel A – Stacked bar: private-RD | shared | private-VP
  Panel B – Per-chromosome sensitivity (varient_piper recall of rare_disease)
"""

import csv
import os
import subprocess
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ── Helpers ───────────────────────────────────────────────────────────────────

BCFTOOLS    = snakemake.params.bcftools
TRUTH_LABEL = snakemake.params.truth_label
QUERY_LABEL = snakemake.params.query_label


def get_chrom_counts(vcf_path: str) -> dict:
    """Return {chrom: count} for all variant records in a VCF."""
    result = subprocess.run(
        [BCFTOOLS, "query", "-f", "%CHROM\n", vcf_path],
        capture_output=True, text=True, check=True,
    )
    counts: dict = defaultdict(int)
    for line in result.stdout.splitlines():
        chrom = line.strip()
        if chrom:
            counts[chrom] += 1
    return dict(counts)


def chrom_sort_key(chrom: str):
    """Sort chromosomes in standard human-genome order."""
    c = chrom.lstrip("cChHrR").lstrip("hr")   # strip "chr" prefix variants
    # Try numeric
    try:
        return (0, int(c), "")
    except ValueError:
        pass
    order = {"X": (1, 0), "Y": (1, 1), "M": (1, 2), "MT": (1, 2)}
    return order.get(c.upper(), (2, 0, c))


def bar_color(sensitivity: float) -> str:
    if sensitivity >= 0.99:
        return "#27AE60"
    elif sensitivity >= 0.95:
        return "#F39C12"
    return "#E74C3C"


# ── Snakemake interface ───────────────────────────────────────────────────────

log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path), exist_ok=True)

with open(log_path, "w") as log:
    try:
        sample = str(snakemake.wildcards.sample)

        fn_counts = get_chrom_counts(snakemake.input.only_a)    # only in truth
        fp_counts = get_chrom_counts(snakemake.input.only_b)    # only in query
        tp_counts = get_chrom_counts(snakemake.input.shared_a)  # shared

        all_chroms = (
            set(fn_counts) | set(fp_counts) | set(tp_counts)
        )
        chroms = sorted(all_chroms, key=chrom_sort_key)
        n = len(chroms)

        fn_vals  = np.array([fn_counts.get(c, 0) for c in chroms])
        tp_vals  = np.array([tp_counts.get(c, 0) for c in chroms])
        fp_vals  = np.array([fp_counts.get(c, 0) for c in chroms])

        total_truth = tp_vals + fn_vals
        with np.errstate(invalid="ignore", divide="ignore"):
            sensitivities = np.where(total_truth > 0, tp_vals / total_truth, np.nan)

        # ── Figure ────────────────────────────────────────────────────────────

        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(max(12, n * 0.55), 11), sharex=True
        )
        fig.suptitle(
            f"Per-Chromosome SNV Concordance  ·  {TRUTH_LABEL}  vs  {QUERY_LABEL}\nSample: {sample}",
            fontsize=13, fontweight="bold",
        )

        x = np.arange(n)
        width = 0.7

        # Panel A – stacked bar
        ax1.bar(x, fn_vals,  width, label=f"Only {TRUTH_LABEL}", color="#E74C3C", alpha=0.85)
        ax1.bar(x, tp_vals,  width, bottom=fn_vals,
                label="Shared",              color="#2ECC71", alpha=0.85)
        ax1.bar(x, fp_vals,  width, bottom=fn_vals + tp_vals,
                label=f"Only {QUERY_LABEL}", color="#3498DB", alpha=0.85)

        ax1.set_ylabel("Variant count", fontsize=11)
        ax1.set_title("Variant counts per chromosome (stacked)", fontsize=11)
        ax1.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda v, _: f"{int(v):,}")
        )
        ax1.legend(loc="upper right", fontsize=9)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)

        # Panel B – per-chromosome sensitivity
        bar_colors = [
            bar_color(s) if not np.isnan(s) else "#CCCCCC"
            for s in sensitivities
        ]
        ax2.bar(x, np.nan_to_num(sensitivities), width,
                color=bar_colors, edgecolor="black", linewidth=0.3)

        ax2.axhline(0.99, color="red",    linestyle="--", linewidth=1.2,
                    alpha=0.7, label="99% threshold")
        ax2.axhline(0.95, color="orange", linestyle="--", linewidth=1.2,
                    alpha=0.7, label="95% threshold")
        ax2.set_ylim(0, 1.08)
        ax2.set_ylabel("Sensitivity (recall)", fontsize=11)
        ax2.set_title(
            f"Per-chromosome sensitivity  ({QUERY_LABEL} recall of {TRUTH_LABEL})",
            fontsize=11,
        )
        ax2.set_xticks(x)
        ax2.set_xticklabels(chroms, rotation=45, ha="right", fontsize=8)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)

        # Colour legend
        legend_patches = [
            mpatches.Patch(color="#27AE60", label="≥ 99%"),
            mpatches.Patch(color="#F39C12", label="95 – 99%"),
            mpatches.Patch(color="#E74C3C", label="< 95%"),
            plt.Line2D([0], [0], color="red",    linestyle="--", label="99% line"),
            plt.Line2D([0], [0], color="orange",  linestyle="--", label="95% line"),
        ]
        ax2.legend(handles=legend_patches, fontsize=8, loc="lower right")

        # ── Annotate chromosomes with low sensitivity ─────────────────────────
        for xi, (chrom, s) in enumerate(zip(chroms, sensitivities)):
            if not np.isnan(s) and s < 0.95:
                ax2.annotate(
                    f"{s:.1%}",
                    xy=(xi, s + 0.01),
                    ha="center", va="bottom", fontsize=7, color="red",
                    fontweight="bold",
                )

        # ── Save ──────────────────────────────────────────────────────────────

        plt.tight_layout()
        os.makedirs(os.path.dirname(snakemake.output.plot), exist_ok=True)
        plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
        plt.close()

        # CSV alongside the plot
        with open(snakemake.output.csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "chrom", "only_truth", "shared", "only_query",
                "total_truth", "total_query", "sensitivity",
            ])
            for i, chrom in enumerate(chroms):
                sens = float(sensitivities[i]) if not np.isnan(sensitivities[i]) else ""
                writer.writerow([
                    chrom,
                    int(fn_vals[i]), int(tp_vals[i]), int(fp_vals[i]),
                    int(tp_vals[i] + fn_vals[i]), int(tp_vals[i] + fp_vals[i]),
                    sens,
                ])

        log.write(f"Chromosomes plotted: {chroms}\n")
        log.write(f"Plot saved: {snakemake.output.plot}\n")
        print(f"Chromosome plot saved to {snakemake.output.plot}")

    except Exception as exc:
        log.write(f"ERROR: {exc}\n")
        raise
