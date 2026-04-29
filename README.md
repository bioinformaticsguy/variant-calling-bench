# Variant Concordance Pipeline

A Snakemake pipeline that compares variant calls between two short-read pipelines
(**rare_disease** and **varient_piper**) and evaluates them against GIAB gold-standard
benchmarks. Covers both SNVs and structural variants (SVs).

---

## Purpose

Two independent variant-calling pipelines produce VCF files for the same samples.
This pipeline answers:

- **How well do the two pipelines agree with each other?** (pipeline-vs-pipeline, sample DNA_01)
- **How well does varient_piper perform against a known ground truth?** (pipeline-vs-benchmark, sample HG002)

---

## Samples and comparisons

| Sample | SNV comparison | SV comparison |
|--------|---------------|--------------|
| **DNA_01** | rare_disease vs varient_piper (no ground truth available) | rare_disease vs varient_piper |
| **HG002** | varient_piper vs GIAB v4.2.1 benchmark | varient_piper vs CMRG SV v1.00 benchmark |

HG002 (NA24385, the GIAB Ashkenazi son) is a control sample with well-characterised
variant call sets, making it ideal for measuring how accurate varient_piper is.

---

## Project layout

```
var_val/
├── Snakefile                  # Pipeline definition (10 rules)
├── config.yaml                # All tunable parameters (edit here to add samples/change paths)
├── input_data/
│   ├── rare_disease/
│   │   └── DNA_01/
│   │       ├── SNV/case_A4842_DNA_01_snv.vcf.gz
│   │       └── SV/case_A4842_DNA_01_sv.vcf.gz
│   ├── varient_piper/
│   │   ├── DNA_01/
│   │   │   ├── SNV/A4842.vcf.gz
│   │   │   └── SV/A4842.sv.vcf.gz
│   │   └── HG002/
│   │       ├── SNV/HG002.vcf.gz
│   │       └── SV/HG002.manta.sv.vcf.gz
│   └── benchmarks/
│       └── HG002/
│           ├── SNV/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz   # GIAB SNV truth
│           ├── SNV/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
│           ├── SV/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz            # CMRG SV truth
│           └── SV/HG002_GRCh38_CMRG_SV_v1.00.bed
├── scripts/
│   ├── calculate_metrics.py   # Count TP/FP/FN from isec output → JSON + CSV
│   ├── plot_concordance.py    # SNV summary plot (counts + metrics bars)
│   ├── plot_chromosome.py     # Per-chromosome sensitivity plot
│   ├── plot_sv_summary.py     # SV overall metrics plot
│   └── plot_sv_by_type.py     # Per-SVTYPE (DEL/INS/DUP/INV/BND) plot
├── envs/
│   ├── plotting.yaml          # matplotlib/numpy environment for SNV scripts
│   └── truvari.yaml           # truvari + matplotlib for SV scripts
└── results/
    └── {sample}/
        ├── filtered/          # Step 1 output – PASS SNVs only
        ├── normalized/        # Step 2 output – split multiallelics
        ├── isec/snv/          # Step 3 output – bcftools isec split files
        ├── metrics/           # Step 4 output – snv_metrics.json + .csv
        ├── plots/             # Steps 5,6,9,10 – all PNG and CSV plots
        └── sv/
            ├── filtered/      # Step 7 output – PASS SVs only
            └── truvari/       # Step 8 output – truvari bench results
```

---

## Running the pipeline

```bash
# Dry run – see what would execute without running anything
snakemake -n --use-conda

# Full run (4 CPU cores)
snakemake --use-conda --cores 4

# Re-run a specific sample's SNV metrics only
snakemake --use-conda --cores 4 results/DNA_01/metrics/snv_metrics.json

# Force re-run even if outputs already exist
snakemake --use-conda --cores 4 --forceall
```

---

## SNV workflow (steps 1–6)

```
input VCFs
    │
    ▼ Step 1: filter_snvs
    │   Keep SNPs only; apply PASS filter; restrict to high-confidence BED (if set)
    │
    ▼ Step 2: normalize
    │   Split multiallelic sites into biallelic records; left-align to reference
    │
    ▼ Step 3: bcftools_isec
    │   Intersect the two normalised VCFs
    │   Outputs:
    │     0000.vcf  – only in truth (rare_disease / GIAB)   → FN
    │     0001.vcf  – only in query (varient_piper)          → FP
    │     0002.vcf  – in both, truth representation          → TP
    │     0003.vcf  – in both, query representation
    │
    ▼ Step 4: calculate_metrics
    │   sensitivity = TP / (TP + FN)   "how much of truth does query recover?"
    │   precision   = TP / (TP + FP)   "how much of query is actually in truth?"
    │   F1          = harmonic mean of sensitivity and precision
    │   → results/{sample}/metrics/snv_metrics.json + .csv
    │
    ▼ Step 5: plot_summary
    │   Two-panel figure: variant count breakdown + metrics bars
    │   → results/{sample}/plots/snv_concordance_summary.png + .csv
    │
    ▼ Step 6: plot_chromosome
        Per-chromosome stacked bar (FN | TP | FP) + sensitivity bar chart
        → results/{sample}/plots/snv_chromosome_concordance.png + .csv
```

## SV workflow (steps 7–10)

```
input SV VCFs
    │
    ▼ Step 7: filter_svs
    │   Apply PASS filter
    │
    ▼ Step 8: truvari_bench
    │   Fuzzy breakpoint matching – SVs don't need to be at exactly the same position
    │   Outputs:
    │     tp-base.vcf.gz  – matched SVs (truth side)   → TP
    │     tp-comp.vcf.gz  – matched SVs (query side)   → TP
    │     fp.vcf.gz       – query-only SVs              → FP
    │     fn.vcf.gz       – truth-only SVs              → FN
    │     summary.json    – overall precision/recall/F1
    │
    ▼ Step 9: plot_sv_summary
    │   Overall SV concordance bars (sensitivity, precision, F1, counts)
    │   → results/{sample}/plots/sv_concordance_summary.png + .csv
    │
    ▼ Step 10: plot_sv_by_type
        Per-SVTYPE (DEL/INS/DUP/INV/BND) breakdown
        → results/{sample}/plots/sv_concordance_by_type.png + .csv
```

---

## Filtering steps explained

### SNV filtering

| Filter | Tool flag | Applied to | Why |
|--------|-----------|-----------|-----|
| SNPs only | `bcftools view -v snps` | Both pipelines | Indels have different error profiles and are handled separately |
| PASS only | `-f PASS,.` | Both pipelines | Removes variants flagged as low-quality by the caller |
| High-confidence BED | `-R {bed}` | HG002 only | Restricts comparison to regions the GIAB benchmark covers with high confidence; avoids penalising callers for uncertain regions |
| Split multiallelics | `bcftools norm -m-any` | Both pipelines | bcftools isec requires biallelic records; two callers may represent the same multiallelic site differently, causing false mismatches |
| Left-align | `bcftools norm -f {ref}` | Both pipelines | Normalises indel representation; currently skipped (reference not set) |

### SV filtering

| Filter | Tool flag | Applied to | Why |
|--------|-----------|-----------|-----|
| PASS only | `bcftools view -f PASS,.` | Both pipelines | Removes low-confidence SV calls |
| Minimum size 50 bp | `--sizemin 50` | Truvari | Standard SV definition; variants < 50 bp are small indels, not structural variants, and are already captured in the SNV comparison |
| Maximum size 50,000 bp | `--sizemax 50000` (Truvari default) | Truvari | Very large SVs (> 50 kb) are harder to validate and may represent copy-number segments rather than discrete events |
| Breakpoint distance 500 bp | `--refdist 500` | Truvari | SV callers report slightly different breakpoint positions for the same event; this allows 500 bp of slop |
| Size similarity 70% | `--pctsize 0.7` | Truvari | Two callers must agree on SV size within 30%; prevents matching a 100 bp deletion to a 10 kb deletion |
| High-confidence BED | `--includebed {bed}` | HG002 only | Restricts SV comparison to the 545 CMRG gene loci (see below) |

---

## GIAB benchmark sets

### SNV benchmark — GIAB v4.2.1

- **File:** `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
- **Coverage:** Whole genome (chromosomes 1–22), ~3.3 billion bp
- **High-confidence BED:** `_noinconsistent.bed` — excludes regions where the GIAB
  consortium's multiple technologies gave inconsistent results
- **Source:** Short-read + long-read + assembly-based integration by GIAB
  (Zook et al., Nature Biotechnology 2019)
- **Use here:** Acts as the truth set for SNV comparison on HG002; variants in the
  BED-excluded regions are not counted as FP or FN

### SV benchmark — CMRG v1.00

- **File:** `HG002_GRCh38_CMRG_SV_v1.00.vcf.gz`
- **CMRG =** Challenging Medically Relevant Genes
- **Coverage:** 545 gene loci (~11 Mb, < 0.4% of the genome)
- **High-confidence BED:** The 545 CMRG regions themselves
- **Source:** Haplotype-resolved assembly of HG002 using PacBio HiFi + ONT data,
  assembled with hifiasm, aligned back to GRCh38
  (Wagner et al., Nature Biotechnology 2022)
- **Why sequence-resolved:** The ALT alleles contain full sequences (not symbolic
  `<DEL>` tags), enabling precise sequence-level matching
- **Important:** Because this BED covers only specific gene windows, most
  genome-wide SV calls from varient_piper fall outside these regions and are
  silently excluded by Truvari's `--includebed`. This is expected and correct —
  the CMRG benchmark is designed only to evaluate calls in those specific genes.

---

## Understanding the metrics

| Metric | Formula | Interpretation |
|--------|---------|---------------|
| **Sensitivity** (Recall) | TP / (TP + FN) | Fraction of truth variants recovered by the query pipeline |
| **Precision** | TP / (TP + FP) | Fraction of query variants confirmed in the truth set |
| **F1** | 2 × P × R / (P + R) | Single number balancing both; closer to 1.0 is better |

For pipeline-vs-pipeline (DNA_01), neither set is a definitive ground truth.
Sensitivity and precision are symmetric — a high F1 means the pipelines largely agree.

For benchmark comparisons (HG002), the GIAB/CMRG set is the ground truth:
- Low sensitivity → varient_piper is missing real variants (false negatives)
- Low precision → varient_piper is calling extra variants not in the benchmark (false positives)

---

## Output files per sample

```
results/{sample}/
├── metrics/
│   ├── snv_metrics.json          # Full metrics with all counts
│   └── snv_metrics.csv           # Flat CSV for cohort aggregation
├── plots/
│   ├── snv_concordance_summary.png / .csv    # Count breakdown + metrics bars
│   ├── snv_chromosome_concordance.png / .csv # Per-chromosome sensitivity
│   ├── sv_concordance_summary.png / .csv     # Overall SV metrics
│   └── sv_concordance_by_type.png / .csv     # Per-SVTYPE breakdown
└── sv/truvari/
    ├── summary.json              # Truvari overall stats
    ├── tp-base.vcf.gz            # TP variants (truth side)
    ├── tp-comp.vcf.gz            # TP variants (query side)
    ├── fp.vcf.gz                 # FP variants (query only)
    └── fn.vcf.gz                 # FN variants (truth only)
```

The CSV files alongside every plot make it easy to:
- Inspect exact numbers behind a plot
- Aggregate results across multiple samples for cohort-level analysis
- Load into Excel / pandas for custom visualisation

---

## Adding a new sample

1. Add the sample name to `snv_samples` (and `sv_samples` if SV data is available) in `config.yaml`
2. Add VCF paths under `input_vcfs.rare_disease`, `input_vcfs.varient_piper` (and the `input_svs` equivalents)
3. Set `high_confidence_bed` to `""` if no BED is available, or provide a BED path
4. Run `snakemake --use-conda --cores 4`

---

## Known limitations

- **HG002 SV counts are very low** because the CMRG BED covers only ~11 Mb of the genome.
  Most varient_piper SV calls are outside those regions and are correctly excluded.
- **varient_piper HG002 SV VCF contains many small indels** (< 50 bp) that are excluded
  by the `sizemin: 50` Truvari filter. Only ~42 calls in that file qualify as true SVs
  by the standard ≥ 50 bp definition.
- **No reference FASTA configured** — bcftools left-alignment and Truvari sequence-similarity
  scoring (`--pctseq`) are both disabled. Set `reference:` in config.yaml to enable them.
- **Sequence similarity matching off** (`pctseq: 0.0`) — Truvari matches SVs by position
  and size only. With a reference FASTA, setting `pctseq: 0.7` would add sequence-level
  confirmation, reducing false matches.
