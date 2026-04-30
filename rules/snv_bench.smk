# ─────────────────────────────────────────────────────────────────────────────
#  snv_bench.smk — SNV benchmarking rules
#
#  Steps:
#    1. filter_snvs    – keep SNVs only; PASS filter on query (not truth)
#    2. normalize_snvs – split multiallelics; left-align if reference given
#    3. bcftools_isec  – split into FN / FP / TP; BED restriction applied here
#    4. calculate_metrics – sensitivity / precision / F1 → JSON + CSV
#
#  Intermediate outputs (filtered/, normalized/) are marked temp() and deleted
#  automatically after use — only isec, metrics and plots are kept.
# ─────────────────────────────────────────────────────────────────────────────


rule filter_snvs:
    """Keep SNVs only. Apply PASS filter to query; truth is never PASS-filtered."""
    input:
        vcf=get_snv_vcf,
    output:
        vcf=temp(f"{RESULTS}/{{callset}}/filtered/snv.vcf.gz"),
        tbi=temp(f"{RESULTS}/{{callset}}/filtered/snv.vcf.gz.tbi"),
    params:
        pass_flag=lambda wc: (
            "-f PASS,."
            if wc.callset != "truth" and config.get("filter_pass", False)
            else ""
        ),
    conda:
        f"{ENVS}/bcftools.yaml"
    log:
        f"{RESULTS}/{{callset}}/logs/filter_snvs.log",
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools view -v snps {params.pass_flag} {input.vcf} 2>{log} \
            | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule normalize_snvs:
    """Split multiallelics (-m-any); left-align against reference if provided."""
    input:
        vcf=f"{RESULTS}/{{callset}}/filtered/snv.vcf.gz",
        tbi=f"{RESULTS}/{{callset}}/filtered/snv.vcf.gz.tbi",
    output:
        vcf=temp(f"{RESULTS}/{{callset}}/normalized/snv_norm.vcf.gz"),
        tbi=temp(f"{RESULTS}/{{callset}}/normalized/snv_norm.vcf.gz.tbi"),
    params:
        ref_flag=lambda wc: (
            f"-f {config['reference']}" if config.get("reference") else ""
        ),
    conda:
        f"{ENVS}/bcftools.yaml"
    log:
        f"{RESULTS}/{{callset}}/logs/normalize_snvs.log",
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools norm -m-any {params.ref_flag} {input.vcf} 2>{log} \
            | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule bcftools_isec:
    """Intersect truth and query; restrict to HG002 high-confidence BED if set.

    bcftools isec output convention:
      0000.vcf – private to truth   (FN: missed by query)
      0001.vcf – private to query   (FP: not in truth)
      0002.vcf – shared, truth side (TP)
      0003.vcf – shared, query side (TP)
    """
    input:
        truth=f"{RESULTS}/truth/normalized/snv_norm.vcf.gz",
        query=f"{RESULTS}/{{pipeline}}/normalized/snv_norm.vcf.gz",
    output:
        only_truth  =f"{RESULTS}/{{pipeline}}/isec/snv/0000.vcf",
        only_query  =f"{RESULTS}/{{pipeline}}/isec/snv/0001.vcf",
        shared_truth=f"{RESULTS}/{{pipeline}}/isec/snv/0002.vcf",
        shared_query=f"{RESULTS}/{{pipeline}}/isec/snv/0003.vcf",
        readme      =f"{RESULTS}/{{pipeline}}/isec/snv/README.txt",
    params:
        outdir  =f"{RESULTS}/{{pipeline}}/isec/snv",
        bed_flag=lambda wc: (
            f"-R {config['truth']['snv']['bed']}"
            if config["truth"]["snv"].get("bed")
            else ""
        ),
    conda:
        f"{ENVS}/bcftools.yaml"
    log:
        f"{RESULTS}/{{pipeline}}/logs/bcftools_isec.log",
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools isec -p {params.outdir} {params.bed_flag} \
            {input.truth} {input.query} 2>{log}
        """


rule calculate_metrics:
    """Compute sensitivity, precision, F1 from isec output."""
    input:
        only_a  =f"{RESULTS}/{{pipeline}}/isec/snv/0000.vcf",
        only_b  =f"{RESULTS}/{{pipeline}}/isec/snv/0001.vcf",
        shared_a=f"{RESULTS}/{{pipeline}}/isec/snv/0002.vcf",
    output:
        metrics=f"{RESULTS}/{{pipeline}}/metrics/snv_metrics.json",
        csv    =f"{RESULTS}/{{pipeline}}/metrics/snv_metrics.csv",
    params:
        bcftools   =config["bcftools_bin"],
        truth_label=config["truth"]["name"],
        query_label=get_pipeline_label,
    conda:
        f"{ENVS}/plotting.yaml"
    log:
        f"{RESULTS}/{{pipeline}}/logs/calculate_metrics.log",
    script:
        "scripts/calculate_metrics.py"
