# ─────────────────────────────────────────────────────────────────────────────
#  common.smk — shared constants and input-resolution helpers
#
#  Included before all other rule files so RESULTS, SNV_PIPELINES, etc.
#  are available everywhere.
# ─────────────────────────────────────────────────────────────────────────────

RESULTS       = config.get("results_dir", "results")
SNV_PIPELINES = config["snv_pipelines"]
SV_PIPELINES  = config["sv_pipelines"]


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
