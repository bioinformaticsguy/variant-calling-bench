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


rule plot_ggtyped_sv_summary:
    """Overall SV concordance for one GGtyped certainty threshold."""
    input:
        summary=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/summary.json",
    output:
        plot=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/plots/sv_concordance_summary.png",
        csv=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/plots/sv_concordance_summary.csv",
    params:
        truth_label=config["truth"]["name"],
        query_label=lambda wc: f"{config['ggtyped'].get('name', 'GGtyped')} cert>={wc.threshold}",
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/logs/plot_sv_summary.log",
    script:
        f"{SCRIPTS}/plot_sv_summary.py"


rule plot_ggtyped_certainty_comparison:
    """Combined CSV and plot comparing all GGtyped certainty thresholds."""
    input:
        summaries=expand(
            f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/summary.json",
            threshold=GGTYPED_THRESHOLDS,
        ),
    output:
        csv=f"{RESULTS}/ggtyped/certainty_thresholds/combined_sv_metrics.csv",
        plot=f"{RESULTS}/ggtyped/certainty_thresholds/combined_sv_metrics.png",
    params:
        thresholds=GGTYPED_THRESHOLDS,
        truth_label=config["truth"]["name"],
        query_label=config.get("ggtyped", {}).get("name", "GGtyped"),
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/ggtyped/certainty_thresholds/logs/plot_certainty_comparison.log",
    script:
        f"{SCRIPTS}/plot_ggtyped_certainty_comparison.py"


rule plot_sv_truthset_summary:
    """Overall SV concordance for one pipeline against one SV truth set."""
    input:
        summary=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/summary.json",
    output:
        plot=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/plots/sv_concordance_summary.png",
        csv=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/plots/sv_concordance_summary.csv",
    params:
        truth_label=get_sv_truth_label,
        query_label=get_pipeline_label,
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/logs/plot_sv_summary.log",
    script:
        f"{SCRIPTS}/plot_sv_summary.py"


rule plot_sv_truthset_by_type:
    """Per-SVTYPE concordance for one pipeline against one SV truth set."""
    input:
        tp_base=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/tp-base.vcf.gz",
        fp=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/fp.vcf.gz",
        fn=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/fn.vcf.gz",
    output:
        plot=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/plots/sv_concordance_by_type.png",
        csv=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/plots/sv_concordance_by_type.csv",
    params:
        bcftools=config["bcftools_bin"],
        truth_label=get_sv_truth_label,
        query_label=get_pipeline_label,
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/logs/plot_sv_by_type.log",
    script:
        f"{SCRIPTS}/plot_sv_by_type.py"


rule plot_ggtyped_sv_truthset_summary:
    """Overall SV concordance for one GGtyped threshold against one SV truth set."""
    input:
        summary=f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/summary.json",
    output:
        plot=f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/plots/sv_concordance_summary.png",
        csv=f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/plots/sv_concordance_summary.csv",
    params:
        truth_label=get_sv_truth_label,
        query_label=lambda wc: f"{config['ggtyped'].get('name', 'GGtyped')} cert>={wc.threshold}",
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/logs/plot_sv_summary.log",
    script:
        f"{SCRIPTS}/plot_sv_summary.py"


rule plot_combined_sv_truthset_comparison:
    """Combined CSV and plot comparing all SV truth sets and query callsets."""
    input:
        summaries=[
            (
                f"{RESULTS}/sv_truths/{item['truthset']}/{item['query']}/sv/truvari/summary.json"
                if item["kind"] == "pipeline"
                else f"{RESULTS}/sv_truths/{item['truthset']}/ggtyped/certainty_thresholds/{item['threshold']}/sv/truvari/summary.json"
            )
            for item in ALL_SV_TRUTH_COMPARISONS
        ],
    output:
        csv=f"{RESULTS}/sv_truths/combined_sv_truthset_metrics.csv",
        plot=f"{RESULTS}/sv_truths/combined_sv_truthset_metrics.png",
    params:
        truthsets=[item["truthset"] for item in ALL_SV_TRUTH_COMPARISONS],
        truth_labels=[
            SV_TRUTH_CONFIGS[item["truthset"]].get("name", item["truthset"])
            for item in ALL_SV_TRUTH_COMPARISONS
        ],
        query_labels=[
            (
                item["query"]
                if item["kind"] == "pipeline"
                else f"{config.get('ggtyped', {}).get('name', 'GGtyped')} cert>={item['threshold']}"
            )
            for item in ALL_SV_TRUTH_COMPARISONS
        ],
        query_kinds=[item["kind"] for item in ALL_SV_TRUTH_COMPARISONS],
        thresholds=[item["threshold"] for item in ALL_SV_TRUTH_COMPARISONS],
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_truths/logs/plot_combined_sv_truthset_comparison.log",
    script:
        f"{SCRIPTS}/plot_combined_sv_truthset_comparison.py"


rule plot_sv_benchmark_summary:
    """Overall SV concordance for one generic truth/query benchmark."""
    input:
        summary=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/summary.json",
    output:
        plot=f"{RESULTS}/sv_benchmarks/{{benchmark}}/plots/sv_concordance_summary.png",
        csv=f"{RESULTS}/sv_benchmarks/{{benchmark}}/plots/sv_concordance_summary.csv",
    params:
        truth_label=get_sv_benchmark_truth_label,
        query_label=get_sv_benchmark_query_label,
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_benchmarks/{{benchmark}}/logs/plot_sv_summary.log",
    script:
        f"{SCRIPTS}/plot_sv_summary.py"


rule plot_sv_benchmark_by_type:
    """Per-SVTYPE concordance for one generic truth/query benchmark."""
    input:
        tp_base=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/tp-base.vcf.gz",
        fp=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/fp.vcf.gz",
        fn=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/fn.vcf.gz",
    output:
        plot=f"{RESULTS}/sv_benchmarks/{{benchmark}}/plots/sv_concordance_by_type.png",
        csv=f"{RESULTS}/sv_benchmarks/{{benchmark}}/plots/sv_concordance_by_type.csv",
    params:
        bcftools=config["bcftools_bin"],
        truth_label=get_sv_benchmark_truth_label,
        query_label=get_sv_benchmark_query_label,
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_benchmarks/{{benchmark}}/logs/plot_sv_by_type.log",
    script:
        f"{SCRIPTS}/plot_sv_by_type.py"


rule plot_combined_sv_benchmark_comparison:
    """Combined CSV and plot for all configured generic SV benchmarks."""
    input:
        summaries=expand(
            f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/summary.json",
            benchmark=SV_BENCHMARKS,
        ),
    output:
        csv=f"{RESULTS}/sv_benchmarks/combined_sv_benchmark_metrics.csv",
        plot=f"{RESULTS}/sv_benchmarks/combined_sv_benchmark_metrics.png",
    params:
        benchmarks=SV_BENCHMARKS,
        benchmark_labels=[
            SV_BENCHMARK_CONFIGS[name].get("name", name)
            for name in SV_BENCHMARKS
        ],
        truth_labels=[
            SV_BENCHMARK_CONFIGS[name]["truth"].get("name", "truth")
            for name in SV_BENCHMARKS
        ],
        query_labels=[
            SV_BENCHMARK_CONFIGS[name]["query"].get("name", "query")
            for name in SV_BENCHMARKS
        ],
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_benchmarks/logs/plot_combined_sv_benchmark_comparison.log",
    script:
        f"{SCRIPTS}/plot_combined_sv_benchmark_comparison.py"
