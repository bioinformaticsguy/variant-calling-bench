# ─────────────────────────────────────────────────────────────────────────────
#  sv_bench.smk — SV benchmarking rules (Truvari)
#
#  Steps:
#    1. filter_svs    – optionally restrict query to PASS SVs
#    2. truvari_bench – fuzzy breakpoint matching; BED restriction from truth
#
#  Wildcard note:
#    {callset}  = "truth" or a pipeline name (filter)
#    {pipeline} = a pipeline name only        (truvari)
# ─────────────────────────────────────────────────────────────────────────────


rule filter_svs:
    """Optionally restrict to PASS SVs. Truth is never PASS-filtered."""
    input:
        vcf=get_sv_vcf,
    output:
        vcf=f"{RESULTS}/{{callset}}/sv/filtered/sv.vcf.gz",
        tbi=f"{RESULTS}/{{callset}}/sv/filtered/sv.vcf.gz.tbi",
    params:
        pass_flag=lambda wc: (
            "-f PASS,."
            if wc.callset != "truth" and config.get("filter_pass", False)
            else ""
        ),
    conda:
        f"{ENVS}/bcftools.yaml"
    log:
        f"logs/{{callset}}/filter_svs.log",
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools view {params.pass_flag} {input.vcf} 2>{log} \
            | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule truvari_bench:
    """Fuzzy-match SVs between truth and query using Truvari.

    Truvari output convention:
      tp-base.vcf.gz – matched SVs as seen in base (truth)   → TP
      tp-comp.vcf.gz – matched SVs as seen in comp (query)   → TP
      fp.vcf.gz      – comp-only SVs                          → FP
      fn.vcf.gz      – base-only SVs                          → FN
      summary.json   – overall precision / recall / F1
    """
    input:
        base=f"{RESULTS}/truth/sv/filtered/sv.vcf.gz",
        comp=f"{RESULTS}/{{pipeline}}/sv/filtered/sv.vcf.gz",
    output:
        summary=f"{RESULTS}/{{pipeline}}/sv/truvari/summary.json",
        tp_base=f"{RESULTS}/{{pipeline}}/sv/truvari/tp-base.vcf.gz",
        tp_comp=f"{RESULTS}/{{pipeline}}/sv/truvari/tp-comp.vcf.gz",
        fp     =f"{RESULTS}/{{pipeline}}/sv/truvari/fp.vcf.gz",
        fn     =f"{RESULTS}/{{pipeline}}/sv/truvari/fn.vcf.gz",
    params:
        outdir  =f"{RESULTS}/{{pipeline}}/sv/truvari",
        refdist =config["truvari"]["refdist"],
        pctsize =config["truvari"]["pctsize"],
        pctseq  =config["truvari"]["pctseq"],
        sizemin =config["truvari"]["sizemin"],
        ref_flag=lambda wc: (
            f"-f {config['reference']}" if config.get("reference") else ""
        ),
        bed_flag=lambda wc: (
            f"--includebed {config['truth']['sv']['bed']}"
            if config["truth"]["sv"].get("bed")
            else ""
        ),
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"logs/{{pipeline}}/truvari_bench.log",
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
