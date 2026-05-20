# ─────────────────────────────────────────────────────────────────────────────
#  sv_bench.smk — SV benchmarking rules (Truvari)
#
#  Steps:
#    1. filter_svs    – optionally restrict query to PASS SVs
#    2. truvari_bench – fuzzy breakpoint matching; BED restriction from truth
#
#  Intermediate outputs (sv/filtered/) are marked temp() and deleted after use.
# ─────────────────────────────────────────────────────────────────────────────


rule filter_svs:
    """Optionally restrict to PASS SVs. Truth is never PASS-filtered."""
    input:
        vcf=get_sv_vcf,
    output:
        vcf=temp(f"{RESULTS}/{{callset}}/sv/filtered/sv.vcf.gz"),
        tbi=temp(f"{RESULTS}/{{callset}}/sv/filtered/sv.vcf.gz.tbi"),
    params:
        pass_flag=lambda wc: (
            "-f PASS,."
            if wc.callset != "truth" and config.get("filter_pass", False)
            else ""
        ),
    conda:
        f"{ENVS}/bcftools.yaml"
    log:
        f"{RESULTS}/{{callset}}/logs/filter_svs.log",
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
        base    =f"{RESULTS}/truth/sv/filtered/sv.vcf.gz",
        base_tbi=f"{RESULTS}/truth/sv/filtered/sv.vcf.gz.tbi",
        comp    =f"{RESULTS}/{{pipeline}}/sv/filtered/sv.vcf.gz",
        comp_tbi=f"{RESULTS}/{{pipeline}}/sv/filtered/sv.vcf.gz.tbi",
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
        f"{RESULTS}/{{pipeline}}/logs/truvari_bench.log",
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


rule filter_ggtyped_certainty:
    """Filter GGtyped annotated SVs by genotype certainty threshold."""
    input:
        vcf=lambda wc: config["ggtyped"]["vcf"],
    output:
        vcf=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/filtered/sv.vcf.gz",
        tbi=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/filtered/sv.vcf.gz.tbi",
    params:
        threshold=lambda wc: wc.threshold,
        keep_genotypes=lambda wc: ",".join(config["ggtyped"].get("keep_genotypes", ["REF/VAR", "VAR/VAR"])),
    conda:
        f"{ENVS}/bcftools.yaml"
    log:
        f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/logs/filter_ggtyped_certainty.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.vcf})
        python {SCRIPTS}/filter_ggtyped_vcf.py \
            --input {input.vcf} \
            --threshold {params.threshold} \
            --keep-genotypes {params.keep_genotypes} 2>{log} \
            | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule truvari_bench_ggtyped_certainty:
    """Benchmark one GGtyped certainty threshold against the SV truth set."""
    input:
        base=f"{RESULTS}/truth/sv/filtered/sv.vcf.gz",
        base_tbi=f"{RESULTS}/truth/sv/filtered/sv.vcf.gz.tbi",
        comp=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/filtered/sv.vcf.gz",
        comp_tbi=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/filtered/sv.vcf.gz.tbi",
    output:
        summary=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/summary.json",
        tp_base=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/tp-base.vcf.gz",
        tp_comp=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/tp-comp.vcf.gz",
        fp=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/fp.vcf.gz",
        fn=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/fn.vcf.gz",
    params:
        outdir=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari",
        refdist=config["truvari"]["refdist"],
        pctsize=config["truvari"]["pctsize"],
        pctseq=config["truvari"]["pctseq"],
        sizemin=config["truvari"]["sizemin"],
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
        f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/logs/truvari_bench.log",
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


rule filter_sv_truthset:
    """Prepare one configured SV truth set for multi-truth benchmarking."""
    input:
        vcf=get_sv_truth_vcf,
    output:
        vcf=temp(f"{RESULTS}/sv_truths/{{truthset}}/truth/sv/filtered/sv.vcf.gz"),
        tbi=temp(f"{RESULTS}/sv_truths/{{truthset}}/truth/sv/filtered/sv.vcf.gz.tbi"),
    conda:
        f"{ENVS}/bcftools.yaml"
    log:
        f"{RESULTS}/sv_truths/{{truthset}}/truth/logs/filter_svs.log",
    shell:
        """
        mkdir -p $(dirname {log})
        bcftools view {input.vcf} 2>{log} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule truvari_bench_sv_truthset:
    """Benchmark a standard SV query pipeline against one configured truth set."""
    input:
        base=f"{RESULTS}/sv_truths/{{truthset}}/truth/sv/filtered/sv.vcf.gz",
        base_tbi=f"{RESULTS}/sv_truths/{{truthset}}/truth/sv/filtered/sv.vcf.gz.tbi",
        comp=f"{RESULTS}/{{pipeline}}/sv/filtered/sv.vcf.gz",
        comp_tbi=f"{RESULTS}/{{pipeline}}/sv/filtered/sv.vcf.gz.tbi",
    output:
        summary=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/summary.json",
        tp_base=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/tp-base.vcf.gz",
        tp_comp=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/tp-comp.vcf.gz",
        fp=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/fp.vcf.gz",
        fn=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari/fn.vcf.gz",
    params:
        outdir=f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/sv/truvari",
        refdist=config["truvari"]["refdist"],
        pctsize=config["truvari"]["pctsize"],
        pctseq=config["truvari"]["pctseq"],
        sizemin=config["truvari"]["sizemin"],
        ref_flag=lambda wc: (
            f"-f {config['reference']}" if config.get("reference") else ""
        ),
        bed_flag=lambda wc: (
            f"--includebed {get_sv_truth_bed(wc)}" if get_sv_truth_bed(wc) else ""
        ),
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_truths/{{truthset}}/{{pipeline}}/logs/truvari_bench.log",
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


rule truvari_bench_ggtyped_sv_truthset:
    """Benchmark one GGtyped threshold against one configured SV truth set."""
    input:
        base=f"{RESULTS}/sv_truths/{{truthset}}/truth/sv/filtered/sv.vcf.gz",
        base_tbi=f"{RESULTS}/sv_truths/{{truthset}}/truth/sv/filtered/sv.vcf.gz.tbi",
        comp=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/filtered/sv.vcf.gz",
        comp_tbi=f"{RESULTS}/ggtyped/certainty_thresholds/{{threshold}}/sv/filtered/sv.vcf.gz.tbi",
    output:
        summary=f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/summary.json",
        tp_base=f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/tp-base.vcf.gz",
        tp_comp=f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/tp-comp.vcf.gz",
        fp=f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/fp.vcf.gz",
        fn=f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari/fn.vcf.gz",
    params:
        outdir=f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/sv/truvari",
        refdist=config["truvari"]["refdist"],
        pctsize=config["truvari"]["pctsize"],
        pctseq=config["truvari"]["pctseq"],
        sizemin=config["truvari"]["sizemin"],
        ref_flag=lambda wc: (
            f"-f {config['reference']}" if config.get("reference") else ""
        ),
        bed_flag=lambda wc: (
            f"--includebed {get_sv_truth_bed(wc)}" if get_sv_truth_bed(wc) else ""
        ),
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_truths/{{truthset}}/ggtyped/certainty_thresholds/{{threshold}}/logs/truvari_bench.log",
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


rule filter_sv_benchmark_truth:
    """Filter the truth VCF for one generic SV benchmark."""
    input:
        vcf=get_sv_benchmark_truth_vcf,
    output:
        vcf=temp(f"{RESULTS}/sv_benchmarks/{{benchmark}}/truth/filtered/sv.vcf.gz"),
        tbi=temp(f"{RESULTS}/sv_benchmarks/{{benchmark}}/truth/filtered/sv.vcf.gz.tbi"),
    params:
        sample=lambda wc: get_sv_benchmark_sample(wc, "truth"),
        sample_flag=lambda wc: (
            f"--sample {get_sv_benchmark_sample(wc, 'truth')}"
            if get_sv_benchmark_sample(wc, "truth")
            else ""
        ),
        pass_flag=lambda wc: (
            "--pass-only"
            if get_sv_benchmark_filter_value(
                wc,
                "truth_pass_only",
                get_sv_benchmark_filter_value(wc, "pass_only", False),
            )
            else ""
        ),
        min_size=lambda wc: get_sv_benchmark_filter_value(wc, "min_size", config["truvari"]["sizemin"]),
        svtypes=lambda wc: ",".join(get_sv_benchmark_filter_value(wc, "svtypes", [])),
        svtypes_flag=lambda wc: (
            "--svtypes " + ",".join(get_sv_benchmark_filter_value(wc, "svtypes", []))
            if get_sv_benchmark_filter_value(wc, "svtypes", [])
            else ""
        ),
    conda:
        f"{ENVS}/bcftools.yaml"
    log:
        f"{RESULTS}/sv_benchmarks/{{benchmark}}/truth/logs/filter_sv.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.vcf})
        python {SCRIPTS}/filter_sv_vcf.py \
            --input {input.vcf} \
            {params.sample_flag} \
            --min-size {params.min_size} \
            {params.svtypes_flag} \
            {params.pass_flag} 2>{log} \
            | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule filter_sv_benchmark_query:
    """Filter the query VCF for one generic SV benchmark."""
    input:
        vcf=get_sv_benchmark_query_vcf,
    output:
        vcf=temp(f"{RESULTS}/sv_benchmarks/{{benchmark}}/query/filtered/sv.vcf.gz"),
        tbi=temp(f"{RESULTS}/sv_benchmarks/{{benchmark}}/query/filtered/sv.vcf.gz.tbi"),
    params:
        sample=lambda wc: get_sv_benchmark_sample(wc, "query"),
        sample_flag=lambda wc: (
            f"--sample {get_sv_benchmark_sample(wc, 'query')}"
            if get_sv_benchmark_sample(wc, "query")
            else ""
        ),
        pass_flag=lambda wc: (
            "--pass-only"
            if get_sv_benchmark_filter_value(
                wc,
                "query_pass_only",
                get_sv_benchmark_filter_value(wc, "pass_only", True),
            )
            else ""
        ),
        min_size=lambda wc: get_sv_benchmark_filter_value(wc, "min_size", config["truvari"]["sizemin"]),
        svtypes=lambda wc: ",".join(get_sv_benchmark_filter_value(wc, "svtypes", [])),
        svtypes_flag=lambda wc: (
            "--svtypes " + ",".join(get_sv_benchmark_filter_value(wc, "svtypes", []))
            if get_sv_benchmark_filter_value(wc, "svtypes", [])
            else ""
        ),
    conda:
        f"{ENVS}/bcftools.yaml"
    log:
        f"{RESULTS}/sv_benchmarks/{{benchmark}}/query/logs/filter_sv.log",
    shell:
        """
        mkdir -p $(dirname {log}) $(dirname {output.vcf})
        python {SCRIPTS}/filter_sv_vcf.py \
            --input {input.vcf} \
            {params.sample_flag} \
            --min-size {params.min_size} \
            {params.svtypes_flag} \
            {params.pass_flag} 2>{log} \
            | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule truvari_bench_sv_benchmark:
    """Benchmark one configured truth/query SV pair."""
    input:
        base=f"{RESULTS}/sv_benchmarks/{{benchmark}}/truth/filtered/sv.vcf.gz",
        base_tbi=f"{RESULTS}/sv_benchmarks/{{benchmark}}/truth/filtered/sv.vcf.gz.tbi",
        comp=f"{RESULTS}/sv_benchmarks/{{benchmark}}/query/filtered/sv.vcf.gz",
        comp_tbi=f"{RESULTS}/sv_benchmarks/{{benchmark}}/query/filtered/sv.vcf.gz.tbi",
    output:
        summary=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/summary.json",
        tp_base=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/tp-base.vcf.gz",
        tp_comp=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/tp-comp.vcf.gz",
        fp=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/fp.vcf.gz",
        fn=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari/fn.vcf.gz",
    params:
        outdir=f"{RESULTS}/sv_benchmarks/{{benchmark}}/sv/truvari",
        refdist=lambda wc: get_sv_benchmark_truvari_value(wc, "refdist", config["truvari"]["refdist"]),
        pctsize=lambda wc: get_sv_benchmark_truvari_value(wc, "pctsize", config["truvari"]["pctsize"]),
        pctseq=lambda wc: get_sv_benchmark_truvari_value(wc, "pctseq", config["truvari"]["pctseq"]),
        pctovl=lambda wc: get_sv_benchmark_truvari_value(wc, "pctovl", None),
        sizemin=lambda wc: get_sv_benchmark_truvari_value(wc, "sizemin", config["truvari"]["sizemin"]),
        ref_flag=lambda wc: (
            f"-f {config['reference']}" if config.get("reference") else ""
        ),
        pctovl_flag=lambda wc: (
            f"--pctovl {get_sv_benchmark_truvari_value(wc, 'pctovl', None)}"
            if get_sv_benchmark_truvari_value(wc, "pctovl", None) is not None
            else ""
        ),
    conda:
        f"{ENVS}/truvari.yaml"
    log:
        f"{RESULTS}/sv_benchmarks/{{benchmark}}/logs/truvari_bench.log",
    shell:
        """
        mkdir -p $(dirname {log})
        rm -rf {params.outdir}
        truvari bench \
            -b {input.base} \
            -c {input.comp} \
            {params.ref_flag} \
            -o {params.outdir} \
            --refdist {params.refdist} \
            --pctsize {params.pctsize} \
            --pctseq  {params.pctseq} \
            {params.pctovl_flag} \
            --sizemin {params.sizemin} 2>{log}
        """
