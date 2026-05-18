"""
plot_sv_summary.py
───────────────────
Overall SV concordance from truvari summary.json.

Two-panel figure (mirrors the SNV summary plot):
  Panel A – Variant count breakdown: FN | TP | FP
  Panel B – Concordance metrics: Sensitivity / Precision / F1
"""

import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path), exist_ok=True)

with open(log_path, "w") as log:
    try:
        with open(snakemake.input.summary) as fh:
            tv = json.load(fh)

        # Truvari key names differ slightly across versions — handle both
        tp  = tv.get("TP-base", tv.get("TP", 0))
        fp  = tv.get("FP", 0)
        fn  = tv.get("FN", 0)
        rec = tv.get("recall",    tv.get("Recall",    0))
        pre = tv.get("precision", tv.get("Precision", 0))
        f1  = tv.get("f1",        tv.get("F1",        0))
        base_total  = tv.get("base cnt", tp + fn)
        comp_total  = tv.get("comp cnt", tp + fp)
        pipeline = getattr(
            snakemake.wildcards,
            "pipeline",
            f"ggtyped_cert_{getattr(snakemake.wildcards, 'threshold', 'NA')}",
        )
        pipeline = str(pipeline)
        truth_label = snakemake.params.truth_label
        query_label = snakemake.params.query_label

        log.write(f"TP={tp}  FP={fp}  FN={fn}\n")
        log.write(f"Recall={rec:.4f}  Precision={pre:.4f}  F1={f1:.4f}\n")

        # ── Figure ────────────────────────────────────────────────────────────

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(
            f"SV Concordance  ·  {truth_label}  vs  {query_label}\n"
            f"Pipeline: {pipeline}  |  Truvari: refdist={snakemake.config['truvari']['refdist']}bp  "
            f"pctsize={snakemake.config['truvari']['pctsize']}  sizemin={snakemake.config['truvari']['sizemin']}bp",
            fontsize=12, fontweight="bold", y=1.03,
        )

        # Panel A – counts
        categories = [f"Only\n{truth_label}\n(FN)", "Shared\n(TP)", f"Only\n{query_label}\n(FP)"]
        values     = [fn, tp, fp]
        colors     = ["#E74C3C", "#2ECC71", "#3498DB"]

        bars = ax1.bar(categories, values, color=colors,
                       edgecolor="black", linewidth=0.5, width=0.5)
        max_val = max(values) if values else 1
        for bar, val in zip(bars, values):
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max_val * 0.012,
                f"{val:,}", ha="center", va="bottom", fontsize=10, fontweight="bold",
            )

        ax1.set_ylabel("Number of SVs", fontsize=11)
        ax1.set_title("SV Count Breakdown", fontsize=11, pad=10)
        ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.text(
            0.98, 0.97,
            f"{truth_label} total:  {base_total:,}\n{query_label} total: {comp_total:,}",
            transform=ax1.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#FDFDE4", alpha=0.85),
        )

        # Panel B – metrics
        def bar_color(v):
            return "#27AE60" if v >= 0.99 else ("#F39C12" if v >= 0.95 else "#E74C3C")

        metric_labels  = ["Sensitivity\n(Recall)", "Precision", "F1 Score"]
        metric_values  = [rec, pre, f1]

        bars2 = ax2.bar(metric_labels, metric_values,
                        color=[bar_color(v) for v in metric_values],
                        edgecolor="black", linewidth=0.5, width=0.5)
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
        ax2.axhline(0.99, color="red",    linestyle="--", alpha=0.5, linewidth=1.2)
        ax2.axhline(0.95, color="orange", linestyle="--", alpha=0.5, linewidth=1.2)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)

        legend_patches = [
            mpatches.Patch(color="#27AE60", label="≥ 99%"),
            mpatches.Patch(color="#F39C12", label="95 – 99%"),
            mpatches.Patch(color="#E74C3C", label="< 95%"),
            plt.Line2D([0], [0], color="red",    linestyle="--", label="99% line"),
            plt.Line2D([0], [0], color="orange",  linestyle="--", label="95% line"),
        ]
        ax2.legend(handles=legend_patches, fontsize=8, loc="lower right")

        plt.tight_layout()
        os.makedirs(os.path.dirname(snakemake.output.plot), exist_ok=True)
        plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
        plt.close()

        # CSV alongside the plot
        with open(snakemake.output.csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "pipeline", "truth_pipeline", "query_pipeline",
                "fn", "tp", "fp", "base_total", "comp_total",
                "sensitivity", "precision", "f1",
            ])
            writer.writerow([
                pipeline, truth_label, query_label,
                fn, tp, fp, base_total, comp_total,
                round(rec, 6), round(pre, 6), round(f1, 6),
            ])

        log.write(f"Plot saved: {snakemake.output.plot}\n")

    except Exception as exc:
        log.write(f"ERROR: {exc}\n")
        raise
