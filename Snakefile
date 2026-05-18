# ─────────────────────────────────────────────────────────────────────────────
#  Variant Calling Benchmark Pipeline
#
#  Truth set : HG002 GIAB benchmark (configured in config/config.yaml)
#  Query     : one or more pipelines listed under `pipelines` in the config
#
#  Usage:
#    snakemake --cores 4
#    snakemake --configfile config/config_test.yaml --cores 4   # test run
# ─────────────────────────────────────────────────────────────────────────────

configfile: "config/config.yaml"

wildcard_constraints:
    pipeline="[^/]+",
    callset="[^/]+",
    truthset="[^/]+",
    threshold="[0-9.]+",

include: "rules/common.smk"
include: "rules/snv_bench.smk"
include: "rules/sv_bench.smk"
include: "rules/plots.smk"


rule all:
    input:
        # ── SNV outputs (one set per SNV pipeline) ────────────────────────────
        expand(f"{RESULTS}/{{pipeline}}/metrics/snv_metrics.json",             pipeline=SNV_PIPELINES),
        expand(f"{RESULTS}/{{pipeline}}/metrics/snv_metrics.csv",              pipeline=SNV_PIPELINES),
        expand(f"{RESULTS}/{{pipeline}}/plots/snv_concordance_summary.png",    pipeline=SNV_PIPELINES),
        expand(f"{RESULTS}/{{pipeline}}/plots/snv_concordance_summary.csv",    pipeline=SNV_PIPELINES),
        expand(f"{RESULTS}/{{pipeline}}/plots/snv_chromosome_concordance.png", pipeline=SNV_PIPELINES),
        expand(f"{RESULTS}/{{pipeline}}/plots/snv_chromosome_concordance.csv", pipeline=SNV_PIPELINES),
        # ── SV outputs (one set per SV pipeline) ─────────────────────────────
        expand(f"{RESULTS}/{{pipeline}}/sv/truvari/summary.json",              pipeline=SV_PIPELINES),
        expand(f"{RESULTS}/{{pipeline}}/plots/sv_concordance_summary.png",     pipeline=SV_PIPELINES),
        expand(f"{RESULTS}/{{pipeline}}/plots/sv_concordance_summary.csv",     pipeline=SV_PIPELINES),
        expand(f"{RESULTS}/{{pipeline}}/plots/sv_concordance_by_type.png",     pipeline=SV_PIPELINES),
        expand(f"{RESULTS}/{{pipeline}}/plots/sv_concordance_by_type.csv",     pipeline=SV_PIPELINES),
        # ── GGtyped certainty-threshold SV comparison ───────────────────────
        expand(
            f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/summary.json",
            threshold=GGTYPED_THRESHOLDS if GGTYPED_ENABLED else [],
        ),
        expand(
            f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/plots/sv_concordance_summary.csv",
            threshold=GGTYPED_THRESHOLDS if GGTYPED_ENABLED else [],
        ),
        expand(
            f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/plots/sv_concordance_summary.png",
            threshold=GGTYPED_THRESHOLDS if GGTYPED_ENABLED else [],
        ),
        f"{RESULTS}/ggtyped/certainty_thresholds/combined_sv_metrics.csv" if GGTYPED_ENABLED else [],
        f"{RESULTS}/ggtyped/certainty_thresholds/combined_sv_metrics.png" if GGTYPED_ENABLED else [],
        # ── Multi-truth SV comparisons: CMRG vs SV_Tier1, etc. ──────────────
        expand(
            f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/summary.json",
            truthset=SV_TRUTHSETS,
            pipeline=SV_PIPELINES,
        ),
        expand(
            f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/plots/sv_concordance_summary.png",
            truthset=SV_TRUTHSETS,
            pipeline=SV_PIPELINES,
        ),
        expand(
            f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/plots/sv_concordance_summary.csv",
            truthset=SV_TRUTHSETS,
            pipeline=SV_PIPELINES,
        ),
        expand(
            f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/plots/sv_concordance_by_type.png",
            truthset=SV_TRUTHSETS,
            pipeline=SV_PIPELINES,
        ),
        expand(
            f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/plots/sv_concordance_by_type.csv",
            truthset=SV_TRUTHSETS,
            pipeline=SV_PIPELINES,
        ),
        expand(
            f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/summary.json",
            truthset=SV_TRUTHSETS if GGTYPED_ENABLED else [],
            threshold=GGTYPED_THRESHOLDS if GGTYPED_ENABLED else [],
        ),
        expand(
            f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/plots/sv_concordance_summary.png",
            truthset=SV_TRUTHSETS if GGTYPED_ENABLED else [],
            threshold=GGTYPED_THRESHOLDS if GGTYPED_ENABLED else [],
        ),
        expand(
            f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/plots/sv_concordance_summary.csv",
            truthset=SV_TRUTHSETS if GGTYPED_ENABLED else [],
            threshold=GGTYPED_THRESHOLDS if GGTYPED_ENABLED else [],
        ),
        f"{RESULTS}/sv_truths/combined_sv_truthset_metrics.csv" if SV_TRUTHSETS else [],
        f"{RESULTS}/sv_truths/combined_sv_truthset_metrics.png" if SV_TRUTHSETS else [],
