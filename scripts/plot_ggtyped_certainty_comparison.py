"""
Compare GGtyped certainty thresholds from Truvari summary.json files.

Specificity needs true negatives, which Truvari summary output does not define
for variant records. The CSV therefore keeps a specificity column as NA and
uses precision as the FP-aware metric in the comparison plot.
"""

import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def metric(summary, *names, default=0):
    for name in names:
        if name in summary:
            return summary[name]
    return default


log_path = str(snakemake.log[0])
os.makedirs(os.path.dirname(log_path), exist_ok=True)

with open(log_path, "w") as log:
    try:
        rows = []
        for threshold, summary_path in zip(snakemake.params.thresholds, snakemake.input.summaries):
            with open(summary_path) as fh:
                tv = json.load(fh)

            tp = metric(tv, "TP-base", "TP", default=0)
            fp = metric(tv, "FP", default=0)
            fn = metric(tv, "FN", default=0)
            sensitivity = metric(tv, "recall", "Recall", default=0)
            precision = metric(tv, "precision", "Precision", default=0)
            f1 = metric(tv, "f1", "F1", default=0)
            base_total = metric(tv, "base cnt", default=tp + fn)
            comp_total = metric(tv, "comp cnt", default=tp + fp)

            rows.append(
                {
                    "threshold": float(threshold),
                    "shared_svs": int(tp),
                    "false_positives": int(fp),
                    "false_negatives": int(fn),
                    "truth_total": int(base_total),
                    "query_total": int(comp_total),
                    "sensitivity": float(sensitivity),
                    "precision": float(precision),
                    "specificity": "NA",
                    "f1": float(f1),
                }
            )

        rows.sort(key=lambda row: row["threshold"])
        os.makedirs(os.path.dirname(snakemake.output.csv), exist_ok=True)

        with open(snakemake.output.csv, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "threshold",
                    "shared_svs",
                    "false_positives",
                    "false_negatives",
                    "truth_total",
                    "query_total",
                    "sensitivity",
                    "precision",
                    "specificity",
                    "f1",
                ],
            )
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {
                        **row,
                        "threshold": f"{row['threshold']:g}",
                        "sensitivity": f"{row['sensitivity']:.6f}",
                        "precision": f"{row['precision']:.6f}",
                        "f1": f"{row['f1']:.6f}",
                    }
                )

        thresholds = [row["threshold"] for row in rows]
        sensitivity = [row["sensitivity"] for row in rows]
        precision = [row["precision"] for row in rows]
        f1 = [row["f1"] for row in rows]
        shared = [row["shared_svs"] for row in rows]
        fp = [row["false_positives"] for row in rows]
        fn = [row["false_negatives"] for row in rows]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9), sharex=True)
        fig.suptitle(
            f"{snakemake.params.query_label} SV certainty threshold comparison\n"
            f"Truth: {snakemake.params.truth_label}",
            fontsize=13,
            fontweight="bold",
        )

        ax1.plot(thresholds, sensitivity, marker="o", linewidth=2, label="Sensitivity")
        ax1.plot(thresholds, precision, marker="s", linewidth=2, label="Precision")
        ax1.plot(thresholds, f1, marker="^", linewidth=2, label="F1")
        ax1.set_ylabel("Score")
        ax1.set_ylim(0, 1.05)
        ax1.grid(axis="y", alpha=0.25)
        ax1.legend(loc="best")
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)

        ax2.plot(thresholds, shared, marker="o", linewidth=2, label="Shared SVs")
        ax2.plot(thresholds, fp, marker="s", linewidth=2, label="False positives")
        ax2.plot(thresholds, fn, marker="^", linewidth=2, label="False negatives")
        ax2.set_xlabel("GGtyped certainty threshold")
        ax2.set_ylabel("SV count")
        ax2.grid(axis="y", alpha=0.25)
        ax2.legend(loc="best")
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)

        ax2.set_xticks(thresholds)
        ax2.set_xticklabels([f"{t:g}" for t in thresholds])

        plt.tight_layout()
        plt.savefig(snakemake.output.plot, dpi=150, bbox_inches="tight")
        plt.close()

        for row in rows:
            log.write(
                f"threshold={row['threshold']:g} TP={row['shared_svs']} FP={row['false_positives']} "
                f"FN={row['false_negatives']} sensitivity={row['sensitivity']:.4f} "
                f"precision={row['precision']:.4f} f1={row['f1']:.4f}\n"
            )

    except Exception as exc:
        log.write(f"ERROR: {exc}\n")
        raise
