"""
Combined SV truth-set comparison plot.

Reads Truvari summary.json files for all configured SV truth sets and query
callsets, writes one CSV, and creates a single comparison figure.
"""

import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def metric(summary, *names, default=0):
    for name in names:
        if name in summary:
            return summary[name]
    return default


def short_label(label):
    return (
        label.replace("HG002 ", "")
        .replace("SV v1.00", "")
        .replace("SV Tier1", "Tier1")
        .strip()
    )


log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path), exist_ok=True)

with open(log_path, "w") as log:
    try:
        rows = []
        for idx, summary_path in enumerate(snakemake.input.summaries):
            with open(summary_path) as fh:
                tv = json.load(fh)

            tp = int(metric(tv, "TP-base", "TP", default=0))
            fp = int(metric(tv, "FP", default=0))
            fn = int(metric(tv, "FN", default=0))
            sensitivity = float(metric(tv, "recall", "Recall", default=0))
            precision = float(metric(tv, "precision", "Precision", default=0))
            f1 = float(metric(tv, "f1", "F1", default=0))
            truth_total = int(metric(tv, "base cnt", default=tp + fn))
            query_total = int(metric(tv, "comp cnt", default=tp + fp))

            row = {
                "truthset": snakemake.params.truthsets[idx],
                "truth_label": snakemake.params.truth_labels[idx],
                "query_label": snakemake.params.query_labels[idx],
                "query_kind": snakemake.params.query_kinds[idx],
                "threshold": snakemake.params.thresholds[idx],
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "truth_total": truth_total,
                "query_total": query_total,
                "sensitivity": sensitivity,
                "precision": precision,
                "f1": f1,
            }
            rows.append(row)

        rows.sort(key=lambda r: (r["truthset"], r["query_kind"], r["query_label"]))
        os.makedirs(os.path.dirname(snakemake.output.csv), exist_ok=True)

        with open(snakemake.output.csv, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "truthset",
                    "truth_label",
                    "query_label",
                    "query_kind",
                    "threshold",
                    "tp",
                    "fp",
                    "fn",
                    "truth_total",
                    "query_total",
                    "sensitivity",
                    "precision",
                    "f1",
                ],
            )
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {
                        **row,
                        "sensitivity": f"{row['sensitivity']:.6f}",
                        "precision": f"{row['precision']:.6f}",
                        "f1": f"{row['f1']:.6f}",
                    }
                )

        labels = [
            f"{short_label(row['truth_label'])} | {row['query_label']}"
            for row in rows
        ]
        y = np.arange(len(rows))
        height = max(7, 0.45 * len(rows) + 2)
        fig, axes = plt.subplots(1, 4, figsize=(18, height), sharey=True)
        fig.suptitle("SV benchmark comparison across HG002 truth sets", fontsize=14, fontweight="bold")

        metric_specs = [
            ("Sensitivity", "sensitivity", "#2E86AB"),
            ("Precision", "precision", "#4CAF50"),
            ("F1", "f1", "#7E57C2"),
        ]
        for ax, (title, key, color) in zip(axes[:3], metric_specs):
            values = [row[key] for row in rows]
            ax.barh(y, values, color=color, alpha=0.85)
            ax.set_title(title)
            ax.set_xlim(0, 1.05)
            ax.grid(axis="x", alpha=0.25)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            for yi, value in zip(y, values):
                ax.text(min(value + 0.015, 1.01), yi, f"{value:.3f}", va="center", fontsize=8)

        tp = np.array([row["tp"] for row in rows])
        fp = np.array([row["fp"] for row in rows])
        fn = np.array([row["fn"] for row in rows])
        axes[3].barh(y, fn, color="#E74C3C", alpha=0.85, label="FN")
        axes[3].barh(y, tp, left=fn, color="#2ECC71", alpha=0.85, label="TP")
        axes[3].barh(y, fp, left=fn + tp, color="#3498DB", alpha=0.85, label="FP")
        axes[3].set_title("Counts")
        axes[3].grid(axis="x", alpha=0.25)
        axes[3].legend(loc="lower right", fontsize=8)
        axes[3].spines["top"].set_visible(False)
        axes[3].spines["right"].set_visible(False)

        axes[0].set_yticks(y)
        axes[0].set_yticklabels(labels, fontsize=9)
        axes[0].invert_yaxis()

        for ax in axes:
            ax.tick_params(axis="y", length=0)

        plt.tight_layout()
        plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
        plt.close()

        for row in rows:
            log.write(
                f"{row['truthset']} {row['query_label']} TP={row['tp']} FP={row['fp']} FN={row['fn']} "
                f"sensitivity={row['sensitivity']:.4f} precision={row['precision']:.4f} f1={row['f1']:.4f}\n"
            )

    except Exception as exc:
        log.write(f"ERROR: {exc}\n")
        raise
