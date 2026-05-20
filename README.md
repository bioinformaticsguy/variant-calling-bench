# Variant Calling Benchmark Pipeline

Snakemake workflow for benchmarking HG002 variant calls against GIAB truth sets.
It supports SNV benchmarking, SV benchmarking with Truvari, GGtyped certainty
threshold sweeps, side-by-side HG002 SV truth-set comparison, and generic
truth/query SV benchmarks such as S021-reseq vs S021.

## Current Comparisons

The active config runs these comparisons:

| Query callset | Variant type | Truth set(s) | Output area |
| --- | --- | --- | --- |
| `GGtyped` HG002 | SV | Active `truth.sv` CMRG block | `results/ggtyped/` |
| `GGtyped` HG002 | SV | CMRG and SV_Tier1 | `results/sv_truths/` |
| configured `sv_benchmarks` | SV | benchmark-specific truth VCF | `results/sv_benchmarks/` |

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

## Generic SV Benchmarks

Use `sv_benchmarks` when you want to compare arbitrary SV callsets, for example
treating a resequenced sample as truth and benchmarking the original sample
against it. The current S021 example is:

```yaml
sv_benchmarks:
  s021_reseq_vs_s021:
    name: "S021 vs S021-reseq"
    truth:
      name: "S021-reseq"
      vcf: "insp/S021-and-S021-reseq/input/sv/S021_reseq.truth.vcf.gz"
      sample: "SFB166S-S021-2_X92-25-X83-25_WGS-Rv-0310-WGS-Rv-0286"
    query:
      name: "S021"
      vcf: "insp/S021-and-S021-reseq/input/sv/S021.vcf.gz"
      sample: "SFB166S-S021_X83-25_WGS-Rv-0286"
    filters:
      pass_only: true
      min_size: 50
      svtypes: ["DEL", "INS", "DUP", "INV", "BND"]
    truvari:
      refdist: 500
      pctsize: 0.7
      pctseq: 0.0
      pctovl: 0.0
      sizemin: 50
```

`truth.sample` and `query.sample` are optional for single-sample VCFs. For
multi-sample VCFs, set them to the exact sample names in the VCF header. The
pipeline keeps only records where the selected sample has a non-reference
genotype.

To add multiple benchmarks in the same run, add more entries:

```yaml
sv_benchmarks:
  s021_reseq_vs_s021:
    ...
  sample2_reseq_vs_sample2:
    name: "sample2 vs sample2-reseq"
    truth:
      name: "sample2-reseq"
      vcf: "input_data/sample2/reseq.truth.vcf.gz"
      sample: "sample2_reseq"
    query:
      name: "sample2"
      vcf: "input_data/sample2/sample2.vcf.gz"
      sample: "sample2"
    filters:
      pass_only: true
      min_size: 50
      svtypes: ["DEL", "INS", "DUP", "INV", "BND"]
```

All configured benchmarks run together through `rule all`.

## Running

```bash
# Dry run
snakemake -n --use-conda --cores 4

# Full run: GGtyped single-truth outputs, multi-truth outputs, all plots
snakemake --use-conda --cores 4

# One generic SV benchmark only
snakemake --use-conda --cores 4 \
  results/sv_benchmarks/s021_reseq_vs_s021/sv/truvari/summary.json

# Combined plot for all configured generic SV benchmarks
snakemake --use-conda --cores 4 \
  results/sv_benchmarks/combined_sv_benchmark_metrics.png \
  results/sv_benchmarks/combined_sv_benchmark_metrics.csv

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

Generic benchmark outputs:

```text
results/sv_benchmarks/{benchmark}/sv/truvari/summary.json
results/sv_benchmarks/{benchmark}/sv/truvari/tp-base.vcf.gz
results/sv_benchmarks/{benchmark}/sv/truvari/tp-comp.vcf.gz
results/sv_benchmarks/{benchmark}/sv/truvari/fp.vcf.gz
results/sv_benchmarks/{benchmark}/sv/truvari/fn.vcf.gz
results/sv_benchmarks/{benchmark}/plots/sv_concordance_summary.png
results/sv_benchmarks/{benchmark}/plots/sv_concordance_by_type.png
results/sv_benchmarks/combined_sv_benchmark_metrics.csv
results/sv_benchmarks/combined_sv_benchmark_metrics.png
```

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

### Generic SV Benchmark

```text
each configured sv_benchmarks pair
  -> filter_sv_benchmark_truth
  -> filter_sv_benchmark_query
  -> truvari_bench_sv_benchmark
  -> plot_sv_benchmark_summary
  -> plot_sv_benchmark_by_type

all generic benchmark summaries
  -> plot_combined_sv_benchmark_comparison
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
