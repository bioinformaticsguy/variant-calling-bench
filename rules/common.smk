# ─────────────────────────────────────────────────────────────────────────────
#  common.smk — shared constants and input-resolution helpers
#
#  Included before all other rule files so RESULTS, SNV_PIPELINES, etc.
#  are available everywhere.
# ─────────────────────────────────────────────────────────────────────────────

import os

RESULTS = config.get("results_dir", "results")

# Derive pipeline lists automatically from whichever pipelines have an snv:/sv:
# key defined — no need to list them separately in the config.
SNV_PIPELINES = [name for name, cfg in config["pipelines"].items() if "snv" in cfg]
SV_PIPELINES  = [name for name, cfg in config["pipelines"].items() if "sv"  in cfg]

GGTYPED_CFG = config.get("ggtyped", {})
GGTYPED_ENABLED = bool(GGTYPED_CFG.get("enabled", False))
GGTYPED_THRESHOLDS = [str(t) for t in GGTYPED_CFG.get("certainty_thresholds", [])]

SV_TRUTH_CONFIGS = config.get("sv_truths", {})
if not SV_TRUTH_CONFIGS and "sv" in config.get("truth", {}):
    SV_TRUTH_CONFIGS = {
        "default": {
            "name": config["truth"].get("name", "truth"),
            "vcf": config["truth"]["sv"]["vcf"],
            "bed": config["truth"]["sv"].get("bed", ""),
        }
    }
SV_TRUTHSETS = list(SV_TRUTH_CONFIGS.keys())
STANDARD_SV_TRUTH_COMPARISONS = [
    {"truthset": truthset, "query": pipeline, "kind": "pipeline", "threshold": ""}
    for truthset in SV_TRUTHSETS
    for pipeline in SV_PIPELINES
]
GGTYPED_SV_TRUTH_COMPARISONS = [
    {"truthset": truthset, "query": "ggtyped", "kind": "ggtyped", "threshold": threshold}
    for truthset in SV_TRUTHSETS
    for threshold in (GGTYPED_THRESHOLDS if GGTYPED_ENABLED else [])
]
ALL_SV_TRUTH_COMPARISONS = STANDARD_SV_TRUTH_COMPARISONS + GGTYPED_SV_TRUTH_COMPARISONS

SV_BENCHMARK_CONFIGS = config.get("sv_benchmarks", {})
SV_BENCHMARKS = list(SV_BENCHMARK_CONFIGS.keys())

# Absolute path to the envs/ directory — conda: directives in included .smk
# files resolve relative to the .smk file's own directory, not the project
# root. Using workflow.basedir (= directory of the main Snakefile) avoids the
# resulting rules/envs/ lookup error.
ENVS    = os.path.join(workflow.basedir, "envs")
SCRIPTS = os.path.join(workflow.basedir, "scripts")


# ── Input helpers (used as `input:` lambdas in rules) ────────────────────────

def get_snv_vcf(wc):
    """Return the SNV VCF for wildcard callset ('truth' or a pipeline name)."""
    if wc.callset == "truth":
        return config["truth"]["snv"]["vcf"]
    return config["pipelines"][wc.callset]["snv"]


def get_sv_vcf(wc):
    """Return the SV VCF for wildcard callset ('truth' or a pipeline name)."""
    if wc.callset == "truth":
        return config["truth"]["sv"]["vcf"]
    return config["pipelines"][wc.callset]["sv"]


def get_sv_truth_vcf(wc):
    """Return the SV truth VCF for a truthset wildcard."""
    return SV_TRUTH_CONFIGS[wc.truthset]["vcf"]


def get_sv_truth_bed(wc):
    """Return the SV truth BED for a truthset wildcard, or an empty string."""
    return SV_TRUTH_CONFIGS[wc.truthset].get("bed", "")


def get_sv_truth_label(wc):
    """Return the display label for a truthset wildcard."""
    return SV_TRUTH_CONFIGS[wc.truthset].get("name", wc.truthset)


def get_pipeline_label(wc):
    """Return the display label for a pipeline — the pipeline key name."""
    return wc.pipeline


def get_sv_benchmark_truth_vcf(wc):
    return SV_BENCHMARK_CONFIGS[wc.benchmark]["truth"]["vcf"]


def get_sv_benchmark_query_vcf(wc):
    return SV_BENCHMARK_CONFIGS[wc.benchmark]["query"]["vcf"]


def get_sv_benchmark_label(wc):
    return SV_BENCHMARK_CONFIGS[wc.benchmark].get("name", wc.benchmark)


def get_sv_benchmark_truth_label(wc):
    truth = SV_BENCHMARK_CONFIGS[wc.benchmark]["truth"]
    return truth.get("name", "truth")


def get_sv_benchmark_query_label(wc):
    query = SV_BENCHMARK_CONFIGS[wc.benchmark]["query"]
    return query.get("name", "query")


def get_sv_benchmark_sample(wc, role):
    return SV_BENCHMARK_CONFIGS[wc.benchmark][role].get("sample", "")


def get_sv_benchmark_filter_value(wc, key, default=None):
    filters = SV_BENCHMARK_CONFIGS[wc.benchmark].get("filters", {})
    return filters.get(key, default)


def get_sv_benchmark_truvari_value(wc, key, default):
    params = dict(config.get("truvari", {}))
    params.update(SV_BENCHMARK_CONFIGS[wc.benchmark].get("truvari", {}))
    return params.get(key, default)
