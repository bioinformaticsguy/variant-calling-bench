# Variant Calling Benchmark Pipeline

Snakemake workflow for benchmarking HG002 variant calls against GIAB truth sets.
It supports SNV benchmarking, SV benchmarking with Truvari, GGtyped certainty
threshold sweeps, and side-by-side SV truth-set comparison.

## Current Comparisons

The active config runs these comparisons:

| Query callset | Variant type | Truth set(s) | Output area |
| --- | --- | --- | --- |
| `GGtyped` HG002 | SV | Active `truth.sv` CMRG block | `results/ggtyped/` |
| `GGtyped` HG002 | SV | CMRG and SV_Tier1 | `results/sv_truths/` |

The default `pipelines:` block is empty because the current `input_data/` tree
contains GGtyped outputs and truth sets, but no VarientPiper VCFs. Additional
SNV/SV callsets can still be added under `pipelines:` or supplied with
`submit_job.sh --pipeline ... --snv ... --sv ...`.

The side-by-side SV truth-set comparison is controlled by `sv_truths` in
`config/config.yaml`:

```yaml
sv_truths:
  cmrg:
    name: "HG002 CMRG SV v1.00"
    vcf: "input_data/truth/HG002/SV/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz"
    bed: "input_data/truth/HG002/SV/HG002_GRCh38_CMRG_SV_v1.00.bed"
  sv_tier1:
    name: "HG002 SV Tier1"
    vcf: "input_data/truth/HG002/SV_Tier1/GRCh38_HG2-T2TQ100-V1.1.vcf.gz"
    bed: "input_data/truth/HG002/SV_Tier1/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"
```

Every configured SV pipeline and every enabled GGtyped certainty threshold are
benchmarked against every truth set in this block.

## Required Inputs

Download the standard HG002 truth sets:

```bash
bash install_truth_set.sh --all
```

At minimum, the multi-truth SV comparison expects:

```text
input_data/truth/HG002/SV/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz
input_data/truth/HG002/SV/HG002_GRCh38_CMRG_SV_v1.00.bed
input_data/truth/HG002/SV_Tier1/GRCh38_HG2-T2TQ100-V1.1.vcf.gz
input_data/truth/HG002/SV_Tier1/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
input_data/GGTyped/out_ggtyper_annotated.vcf.gz
```

For optional non-GGtyped pipeline benchmarking, add that callset under
`pipelines:` or pass it to `submit_job.sh`.

## Running

```bash
# Dry run
snakemake -n --use-conda --cores 4

# Full run: GGtyped single-truth outputs, multi-truth outputs, all plots
snakemake --use-conda --cores 4

# Combined CMRG vs SV_Tier1 plot and CSV only
snakemake --use-conda --cores 4 \
  results/sv_truths/combined_sv_truthset_metrics.png \
  results/sv_truths/combined_sv_truthset_metrics.csv

# One SV_Tier1 benchmark for an optional configured pipeline
snakemake --use-conda --cores 4 \
  results/sv_truths/sv_tier1/{pipeline}/sv/truvari/summary.json

# One SV_Tier1 GGtyped threshold benchmark
snakemake --use-conda --cores 4 \
  results/sv_truths/sv_tier1/ggtyped/certainty_thresholds/0.9/sv/truvari/summary.json
```

## Output Layout

Standard SNV and legacy single-truth SV outputs:

```text
results/{pipeline}/metrics/snv_metrics.json
results/{pipeline}/plots/snv_concordance_summary.png
results/{pipeline}/plots/snv_chromosome_concordance.png
results/{pipeline}/sv/truvari/summary.json
results/{pipeline}/plots/sv_concordance_summary.png
results/{pipeline}/plots/sv_concordance_by_type.png
```

GGtyped single-truth certainty sweep:

```text
results/ggtyped/certainty_thresholds/{threshold}/sv/truvari/summary.json
results/ggtyped/certainty_thresholds/{threshold}/plots/sv_concordance_summary.png
results/ggtyped/certainty_thresholds/combined_sv_metrics.csv
results/ggtyped/certainty_thresholds/combined_sv_metrics.png
```

New multi-truth SV outputs:

```text
results/sv_truths/{truthset}/{pipeline}/sv/truvari/summary.json
results/sv_truths/{truthset}/{pipeline}/plots/sv_concordance_summary.png
results/sv_truths/{truthset}/{pipeline}/plots/sv_concordance_by_type.png

results/sv_truths/{truthset}/ggtyped/certainty_thresholds/{threshold}/sv/truvari/summary.json
results/sv_truths/{truthset}/ggtyped/certainty_thresholds/{threshold}/plots/sv_concordance_summary.png

results/sv_truths/combined_sv_truthset_metrics.csv
results/sv_truths/combined_sv_truthset_metrics.png
```

`truthset` is currently `cmrg` or `sv_tier1`.

## Workflow Summary

### SNV Benchmark

```text
truth/query VCF
  -> filter_snvs
  -> normalize_snvs
  -> bcftools_isec
  -> calculate_metrics
  -> plot_summary
  -> plot_chromosome
```

### Standard SV Benchmark

```text
truth/query VCF
  -> filter_svs
  -> truvari_bench
  -> plot_sv_summary
  -> plot_sv_by_type
```

### Multi-Truth SV Benchmark

```text
each configured sv_truths truth VCF
  -> filter_sv_truthset

each SV query VCF
  -> filter_svs
  -> truvari_bench_sv_truthset for each truth set
  -> plot_sv_truthset_summary
  -> plot_sv_truthset_by_type

each GGtyped threshold VCF
  -> filter_ggtyped_certainty
  -> truvari_bench_ggtyped_sv_truthset for each truth set
  -> plot_ggtyped_sv_truthset_summary

all multi-truth summaries
  -> plot_combined_sv_truthset_comparison
```

## GGtyped Filtering

The GGtyped branch filters `out_ggtyper_annotated.vcf.gz` by:

- `GGT_MAX_CERT >= threshold`
- `GGT_GT` in `REF/VAR` or `VAR/VAR`

Thresholds are configured here:

```yaml
ggtyped:
  enabled: true
  certainty_thresholds: [0.5, 0.7, 0.8, 0.9, 0.95]
  keep_genotypes: ["REF/VAR", "VAR/VAR"]
```

## Metrics

| Metric | Formula | Meaning |
| --- | --- | --- |
| Sensitivity / recall | `TP / (TP + FN)` | Fraction of truth variants recovered by the query |
| Precision | `TP / (TP + FP)` | Fraction of query variants confirmed by truth |
| F1 | `2 * precision * recall / (precision + recall)` | Balance of precision and recall |

Truvari does not define true negatives for variant records, so specificity is not
reported. The combined plots focus on sensitivity, precision, F1, and TP/FP/FN
counts.

## Notes

- CMRG covers medically relevant regions only, so SV recall can look low for
  genome-wide callers because many query calls fall outside the benchmark BED.
- SV_Tier1 gives broader SV truth coverage and is usually the better target for a
  genome-wide HG002 SV comparison.
- `reference` is empty in the active config, so `bcftools norm -f` left-alignment
  and Truvari sequence similarity matching are disabled.
- `truvari.pctseq` is set to `0.0`; set a reference FASTA and raise `pctseq` only
  if sequence-level matching is required.
