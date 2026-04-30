# ─────────────────────────────────────────────────────────────────────────────
#  plots.smk — all visualisation rules
#
#  SNV plots:
#    plot_summary    – stacked bar (counts) + metrics bar chart
#    plot_chromosome – per-chromosome stacked bar + sensitivity bar
#
#  SV plots:
#    plot_sv_summary  – overall SV sensitivity / precision / F1
#    plot_sv_by_type  – per-SVTYPE (DEL/INS/DUP/INV/BND) concordance
# ─────────────────────────────────────────────────────────────────────────────


rule plot_summary:
    """Two-panel SNV summary: variant count breakdown + concordance metrics."""
    input:
        metrics=f"{RESULTS}/{{pipeline}}/metrics/snv_metrics.json",
    output:
        plot=f"{RESULTS}/{{pipeline}}/plots/snv_concordance_summary.png",
        csv =f"{RESULTS}/{{pipeline}}/plots/snv_concordance_summary.csv",
    conda:
        f"{ENVS}/plotting.yaml"
    log:
        f"{RESULTS}/{{pipeline}}/logs/plot_summary.log",
    script:
        f"{SCRIPTS}/plot_concordance.py"


rule plot_chromosome:
    """Per-chromosome stacked bar + sensitivity bar chart."""
    input:
        only_a  =f"{RESULTS}/{{pipeline}}/isec/snv/0000.vcf",
        only_b  =f"{RESULTS}/{{pipeline}}/isec/snv/0001.vcf",
        shared_a=f"{RESULTS}/{{pipeline}}/isec/snv/0002.vcf",
    output:
        plot=f"{RESULTS}/{{pipeline}}/plots/snv_chromosome_concordance.png",
        csv =f"{RESULTS}/{{pipeline}}/plots/snv_chromosome_concordance.csv",
    params:
        bcftools   =config["bcftools_bin"],
        truth_label=config["truth"]["name"],
        query_label=get_pipeline_label,
    conda:
        f"{ENVS}/plotting.yaml"
    log:
        f"{RESULTS}/{{pipeline}}/logs/plot_chromosome.log",
    script:
        f"{SCRIPTS}/plot_chromosome.py"


rule plot_sv_summary:
    """Overall SV concordance metrics from truvari summary.json."""
    input:
        summary=f"{RESULTS}/{{pipeline}}/sv/truvari/summary.json",
    output:
        plot=f"{RESULTS}/{{pipeline}}/plots/sv_concordance_summary.png",
        csv =f"{RESULTS}/{{pipeline}}/plots/sv_concordance_summary.csv",
    params:
        truth_label=config["truth"]["name"],
        query_label=get_pipeline_label,
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/{{pipeline}}/logs/plot_sv_summary.log",
    script:
        f"{SCRIPTS}/plot_sv_summary.py"


rule plot_sv_by_type:
    """Per-SVTYPE (DEL/INS/DUP/INV/BND) concordance breakdown."""
    input:
        tp_base=f"{RESULTS}/{{pipeline}}/sv/truvari/tp-base.vcf.gz",
        fp     =f"{RESULTS}/{{pipeline}}/sv/truvari/fp.vcf.gz",
        fn     =f"{RESULTS}/{{pipeline}}/sv/truvari/fn.vcf.gz",
    output:
        plot=f"{RESULTS}/{{pipeline}}/plots/sv_concordance_by_type.png",
        csv =f"{RESULTS}/{{pipeline}}/plots/sv_concordance_by_type.csv",
    params:
        bcftools   =config["bcftools_bin"],
        truth_label=config["truth"]["name"],
        query_label=get_pipeline_label,
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/{{pipeline}}/logs/plot_sv_by_type.log",
    script:
        f"{SCRIPTS}/plot_sv_by_type.py"
