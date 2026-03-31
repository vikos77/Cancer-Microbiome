#!/usr/bin/env python3
"""
decontam_contigs.py — Contig decontamination filter for CLAUDE pipeline Phase 6.

Removes contigs that:
  1. Hit T2T-CHM13 with ≥90% identity AND alignment ≥ 500 bp (human genomic)
  2. Hit UniVec with ≥95% identity AND alignment ≥ 100 bp (vector contamination)

Short UniVec hits (primer traces at contig ends) are logged but NOT removed.

Usage:
    python decontam_contigs.py \
        --contigs all_contigs.fasta \
        --t2t-blast blast_vs_t2t.tsv \
        --univec-blast blast_vs_univec.tsv \
        --out clean_contigs.fasta \
        --removed removed_contigs.txt \
        --t2t-min-len 500 \
        --t2t-min-pident 90.0 \
        --univec-min-len 100 \
        --univec-min-pident 95.0

Output columns in --removed:
    contig_id  reason  pident  aln_len  qlen  fraction_covered
"""

import argparse
import sys
from collections import defaultdict


def parse_blast6(fname, min_pident, min_alen, label):
    """Parse BLAST -outfmt 6 with qlen as 13th column. Return set of flagged contig IDs."""
    flagged = {}  # contig_id → (reason, pident, alen, qlen)
    if not fname:
        return flagged
    try:
        with open(fname) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 12:
                    continue
                qseqid = parts[0]
                pident = float(parts[2])
                alen = int(parts[3])
                qlen = int(parts[12]) if len(parts) >= 13 else 0
                if pident >= min_pident and alen >= min_alen:
                    frac = alen / qlen if qlen > 0 else 0
                    if qseqid not in flagged:
                        flagged[qseqid] = (label, pident, alen, qlen, frac)
    except FileNotFoundError:
        print(f"WARNING: {fname} not found, skipping {label} filter.", file=sys.stderr)
    return flagged


def parse_fasta(fname):
    """Yield (header, seq) pairs from FASTA file."""
    header = None
    seqs = []
    with open(fname) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seqs)
                header = line[1:]
                seqs = []
            else:
                seqs.append(line)
    if header is not None:
        yield header, ''.join(seqs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contigs',          required=True)
    parser.add_argument('--t2t-blast',        default=None)
    parser.add_argument('--univec-blast',     default=None)
    parser.add_argument('--out',              required=True)
    parser.add_argument('--removed',          required=True)
    parser.add_argument('--t2t-min-len',      type=int,   default=500)
    parser.add_argument('--t2t-min-pident',   type=float, default=90.0)
    parser.add_argument('--univec-min-len',   type=int,   default=100)
    parser.add_argument('--univec-min-pident',type=float, default=95.0)
    args = parser.parse_args()

    # Collect contigs to remove
    t2t_hits   = parse_blast6(args.t2t_blast,   args.t2t_min_pident,   args.t2t_min_len,   'T2T_human')
    uvec_hits  = parse_blast6(args.univec_blast, args.univec_min_pident, args.univec_min_len,'UniVec')
    all_remove = {**t2t_hits, **uvec_hits}

    # Write removed log
    n_removed = 0
    with open(args.removed, 'w') as rf:
        rf.write("contig_id\treason\tpident\taln_len\tqlen\tfrac_covered\n")
        for cid, (reason, pident, alen, qlen, frac) in sorted(all_remove.items()):
            rf.write(f"{cid}\t{reason}\t{pident:.1f}\t{alen}\t{qlen}\t{frac:.3f}\n")
            n_removed += 1

    # Write clean FASTA
    n_total = n_kept = 0
    with open(args.out, 'w') as out_fh:
        for header, seq in parse_fasta(args.contigs):
            contig_id = header.split()[0]
            n_total += 1
            if contig_id in all_remove:
                continue
            out_fh.write(f">{header}\n{seq}\n")
            n_kept += 1

    print(f"Total contigs:   {n_total:>8,}")
    print(f"T2T hits (human):{len(t2t_hits):>8,}")
    print(f"UniVec hits:     {len(uvec_hits):>8,}")
    print(f"Removed total:   {n_removed:>8,}")
    print(f"Kept (clean):    {n_kept:>8,}")
    print(f"Output:          {args.out}")
    print(f"Removed log:     {args.removed}")


if __name__ == '__main__':
    main()
