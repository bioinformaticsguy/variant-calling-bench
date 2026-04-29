"""
plot_concordance.py
────────────────────
Two-panel summary figure:
  Panel A – Variant count breakdown  (only-RD | shared | only-VP)
  Panel B – Concordance metrics bar  (sensitivity, precision, F1)
"""

import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ── Load metrics ─────────────────────────────────────────────────────────────

log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path), exist_ok=True)

with open(log_path, "w") as log:
    try:
        with open(snakemake.input.metrics) as fh:
            data = json.load(fh)

        counts  = data["counts"]
        metrics = data["metrics"]
        sample  = data["sample"]
        truth   = data["truth_pipeline"]
        query   = data["query_pipeline"]

        # ── Figure layout ─────────────────────────────────────────────────────

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(
            f"SNV Concordance  ·  {truth}  vs  {query}\nSample: {sample}",
            fontsize=13,
            fontweight="bold",
            y=1.02,
        )

        # ── Panel A: Variant count breakdown ──────────────────────────────────

        categories = [f"Only\n{truth}", "Shared", f"Only\n{query}"]
        values     = [
            counts[f"only_{truth}"],
            counts["shared"],
            counts[f"only_{query}"],
        ]
        colors = ["#E74C3C", "#2ECC71", "#3498DB"]

        bars = ax1.bar(
            categories, values, color=colors,
            edgecolor="black", linewidth=0.5, width=0.5,
        )

        max_val = max(values) if values else 1
        for bar, val in zip(bars, values):
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max_val * 0.012,
                f"{val:,}",
                ha="center", va="bottom", fontsize=10, fontweight="bold",
            )

        ax1.set_ylabel("Number of Variants", fontsize=11)
        ax1.set_title("Variant Count Breakdown", fontsize=11, pad=10)
        ax1.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{int(x):,}")
        )
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)

        ax1.text(
            0.98, 0.97,
            (
                f"{truth} total:  {counts[f'total_{truth}']:,}\n"
                f"{query} total: {counts[f'total_{query}']:,}"
            ),
            transform=ax1.transAxes, ha="right", va="top",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#FDFDE4", alpha=0.85),
        )

        # ── Panel B: Metrics ──────────────────────────────────────────────────

        metric_labels  = ["Sensitivity\n(Recall)", "Precision", "F1 Score"]
        metric_values  = [
            metrics["sensitivity"],
            metrics["precision"],
            metrics["f1"],
        ]

        def _bar_color(val: float) -> str:
            if val >= 0.99:
                return "#27AE60"
            elif val >= 0.95:
                return "#F39C12"
            return "#E74C3C"

        bar_colors = [_bar_color(v) for v in metric_values]

        bars2 = ax2.bar(
            metric_labels, metric_values,
            color=bar_colors, edgecolor="black", linewidth=0.5, width=0.5,
        )

        for bar, val in zip(bars2, metric_values):
            ax2.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.006,
                f"{val:.4f}\n({val:.2%})",
                ha="center", va="bottom", fontsize=10, fontweight="bold",
            )

        ax2.set_ylim(0, 1.15)
        ax2.set_ylabel("Score", fontsize=11)
        ax2.set_title("Concordance Metrics", fontsize=11, pad=10)
        ax2.axhline(y=0.99, color="red",    linestyle="--", alpha=0.5, linewidth=1.2,
                    label="99% threshold")
        ax2.axhline(y=0.95, color="orange", linestyle="--", alpha=0.5, linewidth=1.2,
                    label="95% threshold")
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

        # ── Save ──────────────────────────────────────────────────────────────

        plt.tight_layout()
        os.makedirs(os.path.dirname(snakemake.output.plot), exist_ok=True)
        plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
        plt.close()

        # CSV alongside the plot
        with open(snakemake.output.csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "sample", "truth_pipeline", "query_pipeline",
                "fn", "fp", "tp", "total_truth", "total_query",
                "sensitivity", "precision", "f1",
            ])
            writer.writerow([
                sample, truth, query,
                counts[f"only_{truth}"], counts[f"only_{query}"], counts["shared"],
                counts[f"total_{truth}"], counts[f"total_{query}"],
                metrics["sensitivity"], metrics["precision"], metrics["f1"],
            ])

        log.write(f"Plot saved: {snakemake.output.plot}\n")
        print(f"Summary plot saved to {snakemake.output.plot}")

    except Exception as exc:
        log.write(f"ERROR: {exc}\n")
        raise
