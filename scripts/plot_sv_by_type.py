"""
plot_sv_by_type.py
───────────────────
Per-SVTYPE concordance breakdown from Truvari TP/FP/FN VCFs.

Extracts SVTYPE from INFO field, normalises aliases (DUP:TANDEM → DUP, etc.),
then generates a two-panel figure:
  Panel A – grouped bar: TP / FP / FN counts per SVTYPE
  Panel B – per-SVTYPE sensitivity and precision bars
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

BCFTOOLS    = snakemake.params.bcftools
TRUTH_LABEL = snakemake.params.truth_label
QUERY_LABEL = snakemake.params.query_label

# Canonical SVTYPE order and normalisation map
SVTYPE_ORDER = ["DEL", "INS", "DUP", "INV", "BND", "OTHER"]
ALIASES = {
    "DUP:TANDEM": "DUP",
    "DUP:INV":    "DUP",
    "CNV":        "DUP",
    "TRA":        "BND",
    "CTX":        "BND",
}


def normalise(raw: str) -> str:
    """Map raw SVTYPE string to a canonical label."""
    t = ALIASES.get(raw, raw)
    return t if t in SVTYPE_ORDER[:-1] else "OTHER"


def infer_from_alt(alt: str, ref: str) -> str:
    """Infer SV type when INFO/SVTYPE is absent (sequence-resolved alleles)."""
    # Symbolic allele: <DEL>, <INS:ME>, etc.
    if alt.startswith("<") and alt.endswith(">"):
        return normalise(alt[1:-1].split(":")[0])
    # Sequence allele: classify by length difference
    diff = len(alt) - len(ref)
    if diff <= -50:
        return "DEL"
    if diff >= 50:
        return "INS"
    return "OTHER"


def get_svtype_counts(vcf_path: str) -> dict:
    """Return {normalised_svtype: count} for a VCF.

    Tries INFO/SVTYPE first.  Falls back to ALT-based inference when the
    field is absent (e.g. sequence-resolved benchmark VCFs like CMRG).
    """
    counts: dict = defaultdict(int)

    # ── Primary: INFO/SVTYPE ──────────────────────────────────────────────────
    r = subprocess.run(
        [BCFTOOLS, "query", "-f", "%INFO/SVTYPE\n", vcf_path],
        capture_output=True, text=True,
    )
    if r.returncode == 0:
        for line in r.stdout.splitlines():
            raw = line.strip()
            if raw and raw != ".":
                counts[normalise(raw)] += 1

    # ── Fallback: infer from ALT + REF alleles ────────────────────────────────
    if not counts:
        r2 = subprocess.run(
            [BCFTOOLS, "query", "-f", "%ALT\t%REF\n", vcf_path],
            capture_output=True, text=True,
        )
        if r2.returncode == 0:
            for line in r2.stdout.splitlines():
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                counts[infer_from_alt(parts[0], parts[1])] += 1

    return dict(counts)


log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path), exist_ok=True)

with open(log_path, "w") as log:
    try:
        tp_counts = get_svtype_counts(snakemake.input.tp_base)
        fp_counts = get_svtype_counts(snakemake.input.fp)
        fn_counts = get_svtype_counts(snakemake.input.fn)
        pipeline  = str(snakemake.wildcards.pipeline)

        # Only show types that appear in at least one file
        present = sorted(
            {t for t in SVTYPE_ORDER if tp_counts.get(t, 0) + fp_counts.get(t, 0) + fn_counts.get(t, 0) > 0},
            key=SVTYPE_ORDER.index,
        )

        tp_vals  = np.array([tp_counts.get(t, 0) for t in present])
        fp_vals  = np.array([fp_counts.get(t, 0) for t in present])
        fn_vals  = np.array([fn_counts.get(t, 0) for t in present])

        total_base = tp_vals + fn_vals
        total_comp = tp_vals + fp_vals

        with np.errstate(invalid="ignore", divide="ignore"):
            sensitivity = np.where(total_base > 0, tp_vals / total_base, np.nan)
            precision   = np.where(total_comp > 0, tp_vals / total_comp, np.nan)

        log.write(f"SV types present: {present}\n")
        for i, t in enumerate(present):
            log.write(
                f"  {t}: TP={tp_vals[i]}  FP={fp_vals[i]}  FN={fn_vals[i]}"
                f"  Sens={sensitivity[i]:.3f}  Prec={precision[i]:.3f}\n"
            )

        # ── Figure ────────────────────────────────────────────────────────────

        n = len(present)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(10, n * 2), 11))
        fig.suptitle(
            f"Per-SVTYPE Concordance  ·  Pipeline: {pipeline}\n"
            f"{TRUTH_LABEL} (truth) vs {QUERY_LABEL} (query)",
            fontsize=13, fontweight="bold",
        )

        x     = np.arange(n)
        w     = 0.26

        # Panel A – grouped bar (TP / FP / FN)
        ax1.bar(x - w, tp_vals, w, label="TP (shared)",                    color="#2ECC71", edgecolor="black", linewidth=0.4)
        ax1.bar(x,     fp_vals, w, label=f"FP (only {QUERY_LABEL})",       color="#3498DB", edgecolor="black", linewidth=0.4)
        ax1.bar(x + w, fn_vals, w, label=f"FN (only {TRUTH_LABEL})",       color="#E74C3C", edgecolor="black", linewidth=0.4)

        # Annotate bars
        for xi, (tp, fp, fn) in enumerate(zip(tp_vals, fp_vals, fn_vals)):
            for offset, val, color in [(-w, tp, "#1E8449"), (0, fp, "#1A5276"), (w, fn, "#922B21")]:
                if val > 0:
                    ax1.text(xi + offset, val + max(tp_vals.max(), 1) * 0.01,
                             f"{val:,}", ha="center", va="bottom", fontsize=8, color=color, fontweight="bold")

        ax1.set_xticks(x)
        ax1.set_xticklabels(present, fontsize=11)
        ax1.set_ylabel("Number of SVs", fontsize=11)
        ax1.set_title("SV counts per type", fontsize=11)
        ax1.legend(fontsize=9, loc="upper right")
        ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)

        # Panel B – sensitivity + precision per type
        def metric_color(v):
            if np.isnan(v):    return "#CCCCCC"
            if v >= 0.99:      return "#27AE60"
            if v >= 0.95:      return "#F39C12"
            return "#E74C3C"

        sens_colors = [metric_color(v) for v in sensitivity]
        prec_colors = [metric_color(v) for v in precision]

        ax2.bar(x - w / 2, np.nan_to_num(sensitivity), w,
                color=sens_colors, edgecolor="black", linewidth=0.4, label="Sensitivity")
        ax2.bar(x + w / 2, np.nan_to_num(precision),   w,
                color=prec_colors, edgecolor="black", linewidth=0.4,
                label="Precision", alpha=0.7, hatch="//")

        # Annotate values
        for xi, (s, p) in enumerate(zip(sensitivity, precision)):
            if not np.isnan(s):
                ax2.text(xi - w / 2, s + 0.01, f"{s:.0%}",
                         ha="center", va="bottom", fontsize=8, fontweight="bold")
            if not np.isnan(p):
                ax2.text(xi + w / 2, p + 0.01, f"{p:.0%}",
                         ha="center", va="bottom", fontsize=8, fontweight="bold")

        ax2.axhline(0.99, color="red",    linestyle="--", linewidth=1.2, alpha=0.6)
        ax2.axhline(0.95, color="orange", linestyle="--", linewidth=1.2, alpha=0.6)
        ax2.set_ylim(0, 1.15)
        ax2.set_xticks(x)
        ax2.set_xticklabels(present, fontsize=11)
        ax2.set_ylabel("Score", fontsize=11)
        ax2.set_title("Sensitivity & Precision per SVTYPE", fontsize=11)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)

        legend_patches = [
            mpatches.Patch(facecolor="#27AE60",  label="≥ 99%"),
            mpatches.Patch(facecolor="#F39C12",  label="95 – 99%"),
            mpatches.Patch(facecolor="#E74C3C",  label="< 95%"),
            plt.Line2D([0], [0], color="red",    linestyle="--", label="99% line"),
            plt.Line2D([0], [0], color="orange", linestyle="--", label="95% line"),
            mpatches.Patch(facecolor="grey",     label="Sensitivity (solid)"),
            mpatches.Patch(facecolor="grey",     hatch="//", label="Precision (hatched)"),
        ]
        ax2.legend(handles=legend_patches, fontsize=8, loc="lower right", ncol=2)

        plt.tight_layout()
        os.makedirs(os.path.dirname(snakemake.output.plot), exist_ok=True)
        plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
        plt.close()

        # CSV alongside the plot
        with open(snakemake.output.csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "svtype", "tp", "fp", "fn",
                "total_truth", "total_query", "sensitivity", "precision",
            ])
            for i, t in enumerate(present):
                sens = float(sensitivity[i]) if not np.isnan(sensitivity[i]) else ""
                prec = float(precision[i])   if not np.isnan(precision[i])   else ""
                writer.writerow([
                    t,
                    int(tp_vals[i]), int(fp_vals[i]), int(fn_vals[i]),
                    int(tp_vals[i] + fn_vals[i]), int(tp_vals[i] + fp_vals[i]),
                    sens, prec,
                ])

        log.write(f"Plot saved: {snakemake.output.plot}\n")

    except Exception as exc:
        log.write(f"ERROR: {exc}\n")
        raise
