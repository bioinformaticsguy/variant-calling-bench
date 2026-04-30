# ─────────────────────────────────────────────────────────────────────────────
#  common.smk — shared constants and input-resolution helpers
#
#  Included before all other rule files so RESULTS, SNV_PIPELINES, etc.
#  are available everywhere.
# ─────────────────────────────────────────────────────────────────────────────

import os

RESULTS       = config.get("results_dir", "results")
SNV_PIPELINES = config["snv_pipelines"]
SV_PIPELINES  = config["sv_pipelines"]

# Absolute path to the envs/ directory — conda: directives in included .smk
# files resolve relative to the .smk file's own directory, not the project
# root. Using workflow.basedir (= directory of the main Snakefile) avoids the
# resulting rules/envs/ lookup error.
ENVS = os.path.join(workflow.basedir, "envs")


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


def get_pipeline_label(wc):
    """Return the human-readable label for a pipeline wildcard."""
    return config["pipelines"][wc.pipeline].get("label", wc.pipeline)
