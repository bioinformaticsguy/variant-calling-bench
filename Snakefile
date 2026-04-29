# ─────────────────────────────────────────────────────────────────────────────
#  Variant Concordance Pipeline  (SNV + SV)
#
#  SNV steps:
#    1. filter_snvs   – keep SNVs only (optionally PASS-only)
#    2. normalize     – split multiallelics (-m-any); left-align if ref given
#    3. bcftools_isec – split into: only-A, only-B, shared-A, shared-B
#    4. calc_metrics  – sensitivity / precision / F1 → JSON
#    5. plot_summary  – stacked bar + metrics bar chart
#    6. plot_chrom    – per-chromosome stacked bar + sensitivity bar
#
#  SV steps:
#    7. filter_svs    – optionally restrict to PASS SVs
#    8. truvari_bench – fuzzy breakpoint matching via Truvari
#    9. plot_sv_summary  – overall SV sensitivity / precision / F1
#   10. plot_sv_by_type – per-SVTYPE (DEL/INS/DUP/INV/BND) concordance
# ─────────────────────────────────────────────────────────────────────────────

configfile: "config.yaml"

SNV_SAMPLES = config["snv_samples"]       # all samples with SNV comparison
SV_SAMPLES  = config["sv_samples"]        # subset that also have SV data
RESULTS     = config.get("results_dir", "results")
PIPELINES   = ["rare_disease", "varient_piper"]
# Plots now live inside the results tree: results/{sample}/plots/


def folder_label(path):
    """Return the top-level folder name from input_data/<folder>/...
    e.g. 'input_data/benchmarks/HG002/...' → 'benchmarks'
         'input_data/rare_disease/DNA_01/...' → 'rare_disease'
    """
    return path.split("/")[1]

def truth_label(wc):
    return folder_label(config["input_vcfs"]["rare_disease"][wc.sample])

def query_label(wc):
    return folder_label(config["input_vcfs"]["varient_piper"][wc.sample])

def sv_truth_label(wc):
    return folder_label(config["input_svs"]["rare_disease"][wc.sample])

def sv_query_label(wc):
    return folder_label(config["input_svs"]["varient_piper"][wc.sample])


# ── Target ────────────────────────────────────────────────────────────────────

rule all:
    input:
        # SNV outputs (all SNV samples)
        expand(f"{RESULTS}/{{sample}}/plots/snv_concordance_summary.png",   sample=SNV_SAMPLES),
        expand(f"{RESULTS}/{{sample}}/plots/snv_concordance_summary.csv",   sample=SNV_SAMPLES),
        expand(f"{RESULTS}/{{sample}}/plots/snv_chromosome_concordance.png", sample=SNV_SAMPLES),
        expand(f"{RESULTS}/{{sample}}/plots/snv_chromosome_concordance.csv", sample=SNV_SAMPLES),
        expand(f"{RESULTS}/{{sample}}/metrics/snv_metrics.json",            sample=SNV_SAMPLES),
        expand(f"{RESULTS}/{{sample}}/metrics/snv_metrics.csv",             sample=SNV_SAMPLES),
        # SV outputs (only samples with SV VCFs)
        expand(f"{RESULTS}/{{sample}}/plots/sv_concordance_summary.png",    sample=SV_SAMPLES),
        expand(f"{RESULTS}/{{sample}}/plots/sv_concordance_summary.csv",    sample=SV_SAMPLES),
        expand(f"{RESULTS}/{{sample}}/plots/sv_concordance_by_type.png",    sample=SV_SAMPLES),
        expand(f"{RESULTS}/{{sample}}/plots/sv_concordance_by_type.csv",    sample=SV_SAMPLES),
        expand(f"{RESULTS}/{{sample}}/sv/truvari/summary.json",             sample=SV_SAMPLES),


# ── 1. Filter to SNVs ─────────────────────────────────────────────────────────

rule filter_snvs:
    """Keep SNVs only; optionally restrict to PASS and/or high-confidence BED regions."""
    input:
        vcf=lambda wc: config["input_vcfs"][wc.pipeline][wc.sample],
    output:
        vcf=f"{RESULTS}/{{sample}}/filtered/{{pipeline}}_snv.vcf.gz",
        tbi=f"{RESULTS}/{{sample}}/filtered/{{pipeline}}_snv.vcf.gz.tbi",
    params:
        pass_flag=lambda wc: "-f PASS,." if config.get("filter_pass", False) else "",
        bed_flag =lambda wc: (
            f"-R {config['high_confidence_bed'][wc.sample]}"
            if config.get("high_confidence_bed", {}).get(wc.sample)
            else ""
        ),
    conda:
        "/root/miniforge3/envs/dsd_snv_pipeline"
    log:
        f"logs/{{sample}}/filter_snvs_{{pipeline}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools view -v snps {params.pass_flag} {params.bed_flag} {input.vcf} 2>{log} \
            | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


# ── 2. Normalize ──────────────────────────────────────────────────────────────

rule normalize:
    """Split multiallelics; left-align if reference is provided."""
    input:
        vcf=f"{RESULTS}/{{sample}}/filtered/{{pipeline}}_snv.vcf.gz",
        tbi=f"{RESULTS}/{{sample}}/filtered/{{pipeline}}_snv.vcf.gz.tbi",
    output:
        vcf=f"{RESULTS}/{{sample}}/normalized/{{pipeline}}_snv_norm.vcf.gz",
        tbi=f"{RESULTS}/{{sample}}/normalized/{{pipeline}}_snv_norm.vcf.gz.tbi",
    params:
        ref_flag=lambda wc: (
            f"-f {config['reference']}" if config.get("reference") else ""
        ),
    conda:
        "/root/miniforge3/envs/dsd_snv_pipeline"
    log:
        f"logs/{{sample}}/normalize_{{pipeline}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools norm -m-any {params.ref_flag} {input.vcf} 2>{log} \
            | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


# ── 3. bcftools isec ──────────────────────────────────────────────────────────
#
#  Output files (bcftools isec convention):
#    0000.vcf  – private to rare_disease   (FN: missed by varient_piper)
#    0001.vcf  – private to varient_piper  (FP: not in rare_disease)
#    0002.vcf  – shared, as in rare_disease  (TP: truth side)
#    0003.vcf  – shared, as in varient_piper (TP: query side)
#    README.txt – mapping of files to callsets

rule bcftools_isec:
    """Find private and shared variants between the two normalised callsets."""
    input:
        rare_disease=f"{RESULTS}/{{sample}}/normalized/rare_disease_snv_norm.vcf.gz",
        varient_piper=f"{RESULTS}/{{sample}}/normalized/varient_piper_snv_norm.vcf.gz",
    output:
        only_a   =f"{RESULTS}/{{sample}}/isec/snv/0000.vcf",
        only_b   =f"{RESULTS}/{{sample}}/isec/snv/0001.vcf",
        shared_a =f"{RESULTS}/{{sample}}/isec/snv/0002.vcf",
        shared_b =f"{RESULTS}/{{sample}}/isec/snv/0003.vcf",
        readme   =f"{RESULTS}/{{sample}}/isec/snv/README.txt",
    params:
        outdir=f"{RESULTS}/{{sample}}/isec/snv",
    conda:
        "/root/miniforge3/envs/dsd_snv_pipeline"
    log:
        f"logs/{{sample}}/bcftools_isec.log",
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools isec -p {params.outdir} \
            {input.rare_disease} {input.varient_piper} 2>{log}
        """


# ── 4. Calculate metrics ──────────────────────────────────────────────────────

rule calculate_metrics:
    """Compute sensitivity, precision, F1 from isec output."""
    input:
        only_a  =f"{RESULTS}/{{sample}}/isec/snv/0000.vcf",
        only_b  =f"{RESULTS}/{{sample}}/isec/snv/0001.vcf",
        shared_a=f"{RESULTS}/{{sample}}/isec/snv/0002.vcf",
    output:
        metrics=f"{RESULTS}/{{sample}}/metrics/snv_metrics.json",
        csv    =f"{RESULTS}/{{sample}}/metrics/snv_metrics.csv",
    params:
        bcftools   =config["bcftools_bin"],
        truth_label=truth_label,
        query_label=query_label,
    conda:
        "envs/plotting.yaml"
    log:
        f"logs/{{sample}}/calculate_metrics.log",
    script:
        "scripts/calculate_metrics.py"


# ── 5. Summary concordance plot ───────────────────────────────────────────────

rule plot_summary:
    """Stacked bar (counts) + concordance metrics bar chart."""
    input:
        metrics=f"{RESULTS}/{{sample}}/metrics/snv_metrics.json",
    output:
        plot=f"{RESULTS}/{{sample}}/plots/snv_concordance_summary.png",
        csv =f"{RESULTS}/{{sample}}/plots/snv_concordance_summary.csv",
    conda:
        "envs/plotting.yaml"
    log:
        f"logs/{{sample}}/plot_summary.log",
    script:
        "scripts/plot_concordance.py"


# ── 6. Per-chromosome concordance plot ───────────────────────────────────────

rule plot_chromosome:
    """Per-chromosome stacked bar + sensitivity bar chart."""
    input:
        only_a  =f"{RESULTS}/{{sample}}/isec/snv/0000.vcf",
        only_b  =f"{RESULTS}/{{sample}}/isec/snv/0001.vcf",
        shared_a=f"{RESULTS}/{{sample}}/isec/snv/0002.vcf",
    output:
        plot=f"{RESULTS}/{{sample}}/plots/snv_chromosome_concordance.png",
        csv =f"{RESULTS}/{{sample}}/plots/snv_chromosome_concordance.csv",
    params:
        bcftools   =config["bcftools_bin"],
        truth_label=truth_label,
        query_label=query_label,
    conda:
        "envs/plotting.yaml"
    log:
        f"logs/{{sample}}/plot_chromosome.log",
    script:
        "scripts/plot_chromosome.py"


# ═════════════════════════════════════════════════════════════════════════════
#  SV RULES
# ═════════════════════════════════════════════════════════════════════════════

# ── 7. Filter SVs ─────────────────────────────────────────────────────────────

rule filter_svs:
    """Optionally restrict SVs to PASS variants before Truvari."""
    input:
        vcf=lambda wc: config["input_svs"][wc.pipeline][wc.sample],
    output:
        vcf=f"{RESULTS}/{{sample}}/sv/filtered/{{pipeline}}_sv.vcf.gz",
        tbi=f"{RESULTS}/{{sample}}/sv/filtered/{{pipeline}}_sv.vcf.gz.tbi",
    params:
        pass_flag=lambda wc: "-f PASS,." if config.get("filter_pass", False) else "",
    conda:
        "/root/miniforge3/envs/dsd_snv_pipeline"
    log:
        f"logs/{{sample}}/filter_svs_{{pipeline}}.log",
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools view {params.pass_flag} {input.vcf} 2>{log} \
            | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """




# ── 8. Truvari bench ──────────────────────────────────────────────────────────
#
#  Output files (truvari convention):
#    tp-base.vcf.gz  – matched SVs as seen in base (rare_disease)   → TP
#    tp-comp.vcf.gz  – matched SVs as seen in comp (varient_piper)  → TP
#    fp.vcf.gz       – comp-only SVs                                 → FP
#    fn.vcf.gz       – base-only SVs                                 → FN
#    summary.json    – overall precision / recall / F1

rule truvari_bench:
    """Fuzzy-match SVs between pipelines using Truvari (allows breakpoint slop)."""
    input:
        base=f"{RESULTS}/{{sample}}/sv/filtered/rare_disease_sv.vcf.gz",
        comp=f"{RESULTS}/{{sample}}/sv/filtered/varient_piper_sv.vcf.gz",
    output:
        summary =f"{RESULTS}/{{sample}}/sv/truvari/summary.json",
        tp_base =f"{RESULTS}/{{sample}}/sv/truvari/tp-base.vcf.gz",
        tp_comp =f"{RESULTS}/{{sample}}/sv/truvari/tp-comp.vcf.gz",
        fp      =f"{RESULTS}/{{sample}}/sv/truvari/fp.vcf.gz",
        fn      =f"{RESULTS}/{{sample}}/sv/truvari/fn.vcf.gz",
    params:
        outdir   =f"{RESULTS}/{{sample}}/sv/truvari",
        refdist  =config["truvari"]["refdist"],
        pctsize  =config["truvari"]["pctsize"],
        pctseq   =config["truvari"]["pctseq"],
        sizemin  =config["truvari"]["sizemin"],
        ref_flag =lambda wc: f"-f {config['reference']}" if config.get("reference") else "",
        bed_flag =lambda wc: (
            f"--includebed {config['sv_high_confidence_bed'][wc.sample]}"
            if config.get("sv_high_confidence_bed", {}).get(wc.sample)
            else ""
        ),
    conda:
        "envs/truvari.yaml"
    log:
        f"logs/{{sample}}/truvari_bench.log",
    shell:
        """
        mkdir -p $(dirname {log})
        rm -rf {params.outdir}
        truvari bench \
            -b {input.base} \
            -c {input.comp} \
            {params.ref_flag} \
            {params.bed_flag} \
            -o {params.outdir} \
            --refdist {params.refdist} \
            --pctsize {params.pctsize} \
            --pctseq  {params.pctseq} \
            --sizemin {params.sizemin} 2>{log}
        """


# ── 9. SV summary plot ────────────────────────────────────────────────────────

rule plot_sv_summary:
    """Overall SV concordance metrics from truvari summary.json."""
    input:
        summary=f"{RESULTS}/{{sample}}/sv/truvari/summary.json",
    output:
        plot=f"{RESULTS}/{{sample}}/plots/sv_concordance_summary.png",
        csv =f"{RESULTS}/{{sample}}/plots/sv_concordance_summary.csv",
    params:
        truth_label=sv_truth_label,
        query_label=sv_query_label,
    conda:
        "envs/truvari.yaml"
    log:
        f"logs/{{sample}}/plot_sv_summary.log",
    script:
        "scripts/plot_sv_summary.py"


# ── 10. SV per-type plot ──────────────────────────────────────────────────────

rule plot_sv_by_type:
    """Per-SVTYPE (DEL/INS/DUP/INV/BND) concordance breakdown."""
    input:
        tp_base=f"{RESULTS}/{{sample}}/sv/truvari/tp-base.vcf.gz",
        fp     =f"{RESULTS}/{{sample}}/sv/truvari/fp.vcf.gz",
        fn     =f"{RESULTS}/{{sample}}/sv/truvari/fn.vcf.gz",
    output:
        plot=f"{RESULTS}/{{sample}}/plots/sv_concordance_by_type.png",
        csv =f"{RESULTS}/{{sample}}/plots/sv_concordance_by_type.csv",
    params:
        bcftools=config["bcftools_bin"],
    conda:
        "envs/truvari.yaml"
    log:
        f"logs/{{sample}}/plot_sv_by_type.log",
    script:
        "scripts/plot_sv_by_type.py"
