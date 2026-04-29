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
