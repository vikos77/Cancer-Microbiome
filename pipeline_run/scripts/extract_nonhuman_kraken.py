#!/usr/bin/env python3
"""
extract_nonhuman_kraken.py
--------------------------
Parse Kraken2 output and extract read pairs NOT classified as a target taxid
(default: 9606, Homo sapiens).

Usage:
    python extract_nonhuman_kraken.py \
        --kraken-output  sample_kraken2.out \
        --r1             sample_R1.fastq.gz \
        --r2             sample_R2.fastq.gz \
        --out-r1         sample_nonhuman_R1.fastq.gz \
        --out-r2         sample_nonhuman_R2.fastq.gz \
        --exclude-taxid  9606

Kraken2 output format (tab-separated):
    C/U  ReadID  TaxID  Length  LCA_kmer_map
    C    read1   9606   150     ...           <- classified as human
    U    read2   0      150     ...           <- unclassified

Logic:
    A pair is EXCLUDED (human) if EITHER mate is classified under --exclude-taxid
    or any of its descendants. For simplicity and speed we exclude any read
    where the taxid field exactly matches --exclude-taxid, plus unclassified
    reads are always KEPT (they are microbial candidates).
"""

import argparse
import gzip
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--kraken-output", required=True,
                   help="Path to Kraken2 output file (not the report)")
    p.add_argument("--r1", required=True,
                   help="Input R1 FASTQ (gzipped)")
    p.add_argument("--r2", required=True,
                   help="Input R2 FASTQ (gzipped)")
    p.add_argument("--out-r1", required=True,
                   help="Output R1 FASTQ (gzipped)")
    p.add_argument("--out-r2", required=True,
                   help="Output R2 FASTQ (gzipped)")
    p.add_argument("--exclude-taxid", type=int, default=9606,
                   help="Taxid to exclude (default: 9606 = Homo sapiens)")
    return p.parse_args()


def load_human_read_ids(kraken_output, exclude_taxid):
    """
    Return set of read IDs classified as exclude_taxid.
    Kraken2 output: col 0 = C/U, col 1 = read_id, col 2 = taxid
    Strip /1 /2 suffixes for consistent lookup.
    """
    human_ids = set()
    classified = 0
    unclassified = 0

    with open(kraken_output) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            status, read_id, taxid_str = parts[0], parts[1], parts[2]

            # Normalise read name: strip /1 /2 suffix
            base_id = read_id.rstrip("/12").rstrip()
            if base_id.endswith(("/1", "/2")):
                base_id = base_id[:-2]

            try:
                taxid = int(taxid_str)
            except ValueError:
                continue

            if status == "C":
                classified += 1
                if taxid == exclude_taxid:
                    human_ids.add(base_id)
            else:
                unclassified += 1

    total = classified + unclassified
    human_frac = len(human_ids) / max(total, 1) * 100
    print(f"  Kraken2 reads parsed:    {total:,}", file=sys.stderr)
    print(f"  Classified:              {classified:,}", file=sys.stderr)
    print(f"  Unclassified:            {unclassified:,}", file=sys.stderr)
    print(f"  Human (taxid {exclude_taxid}): {len(human_ids):,} ({human_frac:.2f}%)",
          file=sys.stderr)
    return human_ids


def fastq_records(fh):
    """Yield (header, seq, plus, qual) tuples from a FASTQ filehandle."""
    while True:
        header = fh.readline()
        if not header:
            break
        seq  = fh.readline()
        plus = fh.readline()
        qual = fh.readline()
        yield header, seq, plus, qual


def read_name_from_header(header):
    """Extract base read name from FASTQ @header line, strip /1 /2."""
    name = header.lstrip("@").split()[0]
    if name.endswith(("/1", "/2")):
        name = name[:-2]
    return name


def filter_pairs(r1_in, r2_in, r1_out, r2_out, human_ids):
    """Stream R1/R2 FASTQ pairs, write those NOT in human_ids."""
    kept = 0
    removed = 0

    open_r1_in  = gzip.open(r1_in,  "rt")
    open_r2_in  = gzip.open(r2_in,  "rt")
    open_r1_out = gzip.open(r1_out, "wt")
    open_r2_out = gzip.open(r2_out, "wt")

    try:
        for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(
                fastq_records(open_r1_in),
                fastq_records(open_r2_in)):

            name = read_name_from_header(h1)

            if name in human_ids:
                removed += 1
            else:
                open_r1_out.write(h1 + s1 + p1 + q1)
                open_r2_out.write(h2 + s2 + p2 + q2)
                kept += 1
    finally:
        for fh in (open_r1_in, open_r2_in, open_r1_out, open_r2_out):
            fh.close()

    return kept, removed


def main():
    args = parse_args()

    print(f"Loading Kraken2 output: {args.kraken_output}", file=sys.stderr)
    human_ids = load_human_read_ids(args.kraken_output, args.exclude_taxid)

    print(f"\nFiltering FASTQ pairs...", file=sys.stderr)
    kept, removed = filter_pairs(
        args.r1, args.r2,
        args.out_r1, args.out_r2,
        human_ids
    )

    total = kept + removed
    print(f"\n  Pairs written (non-human): {kept:,} ({kept/max(total,1)*100:.2f}%)",
          file=sys.stderr)
    print(f"  Pairs removed (human):     {removed:,} ({removed/max(total,1)*100:.2f}%)",
          file=sys.stderr)
    print(f"  Output R1: {args.out_r1}", file=sys.stderr)
    print(f"  Output R2: {args.out_r2}", file=sys.stderr)


if __name__ == "__main__":
    main()
