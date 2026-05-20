#!/usr/bin/env python3
"""Filter an SV VCF for benchmark-ready records."""

import argparse
import gzip
import signal
import sys


signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_info(info_text):
    info = {}
    if info_text == ".":
        return info
    for item in info_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = True
    return info


def has_alt_gt(gt):
    if not gt or gt == ".":
        return False
    alleles = gt.replace("|", "/").split("/")
    if any(allele == "." for allele in alleles):
        return False
    return any(allele != "0" for allele in alleles)


def first_int(value):
    try:
        return int(str(value).split(",")[0])
    except (TypeError, ValueError):
        return None


def sv_size(info, chrom_pos, ref, alt):
    for key in ("SVLEN", "INSLEN"):
        value = first_int(info.get(key))
        if value is not None:
            return abs(value)
    end = first_int(info.get("END"))
    if end is not None:
        try:
            return abs(end - int(chrom_pos))
        except ValueError:
            return None
    if alt.startswith("<") and alt.endswith(">"):
        return None
    return abs(len(alt) - len(ref))


def svtype(info, ref, alt):
    raw = info.get("SVTYPE")
    if raw and raw != ".":
        return str(raw).split(",")[0]
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1].split(":")[0]
    diff = len(alt) - len(ref)
    if diff <= -50:
        return "DEL"
    if diff >= 50:
        return "INS"
    return "OTHER"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample", default="")
    parser.add_argument("--min-size", type=int, default=50)
    parser.add_argument("--svtypes", default="")
    parser.add_argument("--pass-only", action="store_true")
    args = parser.parse_args()

    keep_svtypes = {x.strip() for x in args.svtypes.split(",") if x.strip()}
    sample_col = None

    with open_text(args.input) as handle:
        for line in handle:
            if line.startswith("##"):
                sys.stdout.write(line)
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                if args.sample:
                    if args.sample not in header[9:]:
                        raise SystemExit(f"Sample {args.sample} was not found in {args.input}")
                    sample_col = header.index(args.sample)
                sys.stdout.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            if args.pass_only and parts[6] not in ("PASS", "."):
                continue

            info = parse_info(parts[7])
            record_type = svtype(info, parts[3], parts[4])
            if keep_svtypes and record_type not in keep_svtypes:
                continue

            size = sv_size(info, parts[1], parts[3], parts[4])
            if record_type != "BND" and (size is None or size < args.min_size):
                continue

            if args.sample:
                if len(parts) <= sample_col:
                    continue
                values = dict(zip(parts[8].split(":"), parts[sample_col].split(":")))
                if not has_alt_gt(values.get("GT", ".")):
                    continue

            sys.stdout.write(line)


if __name__ == "__main__":
    main()
