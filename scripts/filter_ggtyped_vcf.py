#!/usr/bin/env python3
"""Filter a GGtyped annotated VCF by certainty and genotyped call state."""

import argparse
import gzip
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_info(info):
    fields = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            fields[key] = value
        else:
            fields[item] = True
    return fields


def get_certainty(info_fields):
    value = info_fields.get("GGT_MAX_CERT")
    if value in (None, ".", ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def get_ggtyped_genotype(format_field, sample_field):
    if not format_field or not sample_field:
        return None
    keys = format_field.split(":")
    values = sample_field.split(":")
    if "GGT_GT" not in keys:
        return None
    idx = keys.index("GGT_GT")
    if idx >= len(values):
        return None
    return values[idx]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--threshold", type=float, required=True)
    parser.add_argument("--keep-genotypes", default="REF/VAR,VAR/VAR")
    args = parser.parse_args()

    keep_genotypes = {g.strip() for g in args.keep_genotypes.split(",") if g.strip()}
    total = kept = missing_certainty = missing_genotype = below_threshold = filtered_genotype = 0

    with open_text(args.input) as fh:
        for line in fh:
            if line.startswith("#"):
                sys.stdout.write(line)
                continue

            total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            info_fields = parse_info(parts[7])
            certainty = get_certainty(info_fields)
            if certainty is None:
                missing_certainty += 1
                continue
            if certainty < args.threshold:
                below_threshold += 1
                continue

            genotype = get_ggtyped_genotype(parts[8], parts[9])
            if genotype is None:
                missing_genotype += 1
                continue
            if genotype not in keep_genotypes:
                filtered_genotype += 1
                continue

            kept += 1
            sys.stdout.write(line)

    sys.stderr.write(
        "\n".join(
            [
                f"input_records={total}",
                f"kept_records={kept}",
                f"threshold={args.threshold:g}",
                f"keep_genotypes={','.join(sorted(keep_genotypes))}",
                f"missing_certainty={missing_certainty}",
                f"missing_ggtyped_genotype={missing_genotype}",
                f"below_threshold={below_threshold}",
                f"filtered_genotype={filtered_genotype}",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
