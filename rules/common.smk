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


def get_pipeline_label(wc):
    """Return the display label for a pipeline — the pipeline key name."""
    return wc.pipeline
