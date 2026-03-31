#!/usr/bin/env python3
"""
viral_consensus.py — Phase 7 consensus viral contig caller

Applies the CLAUDE pipeline ≥2/3 tool consensus rule for viral contig classification.
Three independent tools with orthogonal algorithms must agree before a contig is
called viral.

Tools and thresholds:
    geNomad  : virus_score >= 0.7
    DVF      : score >= 0.9 AND pvalue <= 0.05
    VirSorter2: max_score >= 0.5

Usage:
    python viral_consensus.py \
        --genomad   <genomad_virus_summary.tsv> \
        --dvf       <dvfpred.txt> \
        --vs2       <final-viral-score.tsv>     [optional] \
        --contigs   <clean_contigs.fasta> \
        --out-fasta <consensus_viral_contigs.fasta> \
        --out-table <viral_classification_summary.tsv>
"""

import argparse
import sys
import os

def load_genomad(tsv_path, min_score=0.7):
    """
    Load geNomad virus summary TSV.
    Returns dict: contig_name -> {'genomad_score': float, 'genomad_hit': bool, 'taxonomy': str}
    """
    hits = {}
    if not os.path.exists(tsv_path):
        print(f"[WARN] geNomad file not found: {tsv_path}", file=sys.stderr)
        return hits

    with open(tsv_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        # Expected columns: seq_name, length, topology, coordinates, n_genes,
        #   genetic_code, virus_score, fdr, n_hallmarks, marker_enrichment, taxonomy
        try:
            score_idx = header.index('virus_score')
            tax_idx   = header.index('taxonomy')
        except ValueError:
            # Fallback: score at col 6, taxonomy at col 10
            score_idx = 6
            tax_idx   = 10

        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= score_idx:
                continue
            name  = parts[0]
            score = float(parts[score_idx])
            tax   = parts[tax_idx] if len(parts) > tax_idx else 'Unclassified'
            hits[name] = {
                'genomad_score': score,
                'genomad_hit':   score >= min_score,
                'genomad_tax':   tax,
            }
    return hits


def load_dvf(txt_path, min_score=0.9, max_pvalue=0.05):
    """
    Load DeepVirFinder output.
    Returns dict: contig_name -> {'dvf_score': float, 'dvf_pvalue': float, 'dvf_hit': bool}
    """
    hits = {}
    if not os.path.exists(txt_path):
        print(f"[WARN] DVF file not found: {txt_path}", file=sys.stderr)
        return hits

    with open(txt_path) as f:
        f.readline()  # skip header: name\tlen\tscore\tpvalue
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            name   = parts[0]
            score  = float(parts[2])
            pvalue = float(parts[3])
            hits[name] = {
                'dvf_score':  score,
                'dvf_pvalue': pvalue,
                'dvf_hit':    (score >= min_score and pvalue <= max_pvalue),
            }
    return hits


def load_vs2(tsv_path, min_score=0.5):
    """
    Load VirSorter2 final-viral-score.tsv.
    Returns dict: contig_name -> {'vs2_score': float, 'vs2_hit': bool, 'vs2_group': str}
    Column names vary by VS2 version; handles both.
    """
    hits = {}
    if not os.path.exists(tsv_path):
        print(f"[WARN] VS2 file not found: {tsv_path}", file=sys.stderr)
        return hits

    with open(tsv_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        # Typical columns: seqname, dsDNAphage, NCLDV, RNA, ssDNA, lavidaviridae,
        #   max_score, max_score_group, length, hallmark, viral, cellular
        try:
            score_idx = header.index('max_score')
            group_idx = header.index('max_score_group')
        except ValueError:
            score_idx = 6
            group_idx = 7

        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= score_idx:
                continue
            # VS2 appends ||full or ||provirus to seqname — strip it for matching
            raw_name = parts[0]
            name     = raw_name.split('||')[0]
            score    = float(parts[score_idx])
            group    = parts[group_idx] if len(parts) > group_idx else 'unknown'
            # Keep highest score if contig appears multiple times (provirus segments)
            if name not in hits or score > hits[name]['vs2_score']:
                hits[name] = {
                    'vs2_score': score,
                    'vs2_hit':   score >= min_score,
                    'vs2_group': group,
                }
    return hits


def read_fasta_names(fasta_path):
    """Return list of (header_line, sequence_lines) tuples, preserving order."""
    records = []
    header  = None
    seqlines = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    records.append((header, seqlines))
                header   = line
                seqlines = []
            else:
                seqlines.append(line)
    if header is not None:
        records.append((header, seqlines))
    return records


def main():
    parser = argparse.ArgumentParser(
        description='CLAUDE Phase 7 — viral consensus caller (≥2/3 tool agreement)'
    )
    parser.add_argument('--genomad',   required=True,  help='geNomad virus_summary.tsv')
    parser.add_argument('--dvf',       required=True,  help='DVF dvfpred.txt')
    parser.add_argument('--vs2',       default=None,   help='VirSorter2 final-viral-score.tsv (optional)')
    parser.add_argument('--contigs',   required=True,  help='Input FASTA (clean_contigs.fasta)')
    parser.add_argument('--out-fasta', required=True,  help='Output FASTA of consensus viral contigs')
    parser.add_argument('--out-table', required=True,  help='Output TSV summary table')
    parser.add_argument('--genomad-min-score', type=float, default=0.7)
    parser.add_argument('--dvf-min-score',     type=float, default=0.9)
    parser.add_argument('--dvf-max-pvalue',    type=float, default=0.05)
    parser.add_argument('--vs2-min-score',     type=float, default=0.5)
    parser.add_argument('--min-votes',         type=int,   default=2,
                        help='Minimum tool votes to call a contig viral (default: 2)')
    args = parser.parse_args()

    # --- Load tool results ---
    print('[INFO] Loading geNomad results...', file=sys.stderr)
    genomad = load_genomad(args.genomad, min_score=args.genomad_min_score)

    print('[INFO] Loading DVF results...', file=sys.stderr)
    dvf = load_dvf(args.dvf, min_score=args.dvf_min_score, max_pvalue=args.dvf_max_pvalue)

    vs2 = {}
    if args.vs2:
        print('[INFO] Loading VirSorter2 results...', file=sys.stderr)
        vs2 = load_vs2(args.vs2, min_score=args.vs2_min_score)
    else:
        print('[INFO] VS2 file not provided — running with 2-tool mode (geNomad + DVF)', file=sys.stderr)

    n_tools = 2 + (1 if args.vs2 else 0)
    print(f'[INFO] Tools loaded: {n_tools}  |  Min votes for consensus: {args.min_votes}',
          file=sys.stderr)

    # --- Read contigs ---
    print('[INFO] Reading contigs...', file=sys.stderr)
    records = read_fasta_names(args.contigs)
    print(f'[INFO] Total contigs: {len(records)}', file=sys.stderr)

    # --- Build per-contig vote table ---
    table_rows = []
    consensus_names = set()

    for header, seqlines in records:
        name = header[1:].split()[0]  # strip '>' and take first token

        g_info = genomad.get(name, {'genomad_score': 'NA', 'genomad_hit': False, 'genomad_tax': 'NA'})
        d_info = dvf.get(name,     {'dvf_score': 'NA',  'dvf_pvalue': 'NA', 'dvf_hit': False})
        v_info = vs2.get(name,     {'vs2_score': 'NA',  'vs2_hit': False,   'vs2_group': 'NA'})

        votes = sum([
            int(g_info['genomad_hit']),
            int(d_info['dvf_hit']),
            int(v_info['vs2_hit']),
        ])

        is_viral = votes >= args.min_votes

        if is_viral:
            consensus_names.add(name)

        table_rows.append({
            'contig':          name,
            'genomad_score':   g_info['genomad_score'],
            'genomad_hit':     int(g_info['genomad_hit']),
            'genomad_tax':     g_info['genomad_tax'],
            'dvf_score':       d_info['dvf_score'],
            'dvf_pvalue':      d_info['dvf_pvalue'],
            'dvf_hit':         int(d_info['dvf_hit']),
            'vs2_score':       v_info['vs2_score'],
            'vs2_hit':         int(v_info['vs2_hit']),
            'vs2_group':       v_info['vs2_group'],
            'total_votes':     votes,
            'n_tools':         n_tools,
            'consensus_viral': int(is_viral),
        })

    # --- Write consensus FASTA ---
    os.makedirs(os.path.dirname(os.path.abspath(args.out_fasta)), exist_ok=True)
    n_written = 0
    with open(args.out_fasta, 'w') as fw:
        for header, seqlines in records:
            name = header[1:].split()[0]
            if name in consensus_names:
                fw.write(header + '\n')
                fw.write('\n'.join(seqlines) + '\n')
                n_written += 1

    # --- Write summary table ---
    os.makedirs(os.path.dirname(os.path.abspath(args.out_table)), exist_ok=True)
    cols = ['contig', 'genomad_score', 'genomad_hit', 'genomad_tax',
            'dvf_score', 'dvf_pvalue', 'dvf_hit',
            'vs2_score', 'vs2_hit', 'vs2_group',
            'total_votes', 'n_tools', 'consensus_viral']
    with open(args.out_table, 'w') as fw:
        fw.write('\t'.join(cols) + '\n')
        for row in table_rows:
            fw.write('\t'.join(str(row[c]) for c in cols) + '\n')

    # --- Summary stats ---
    votes_hist = {0: 0, 1: 0, 2: 0, 3: 0}
    for row in table_rows:
        votes_hist[row['total_votes']] = votes_hist.get(row['total_votes'], 0) + 1

    print(f'\n[RESULT] Consensus viral contigs (≥{args.min_votes}/{n_tools} tools): {n_written}', file=sys.stderr)
    print(f'[RESULT] Vote distribution:', file=sys.stderr)
    for v in sorted(votes_hist):
        print(f'         {v}/{n_tools} votes: {votes_hist[v]} contigs', file=sys.stderr)

    genomad_hits = sum(1 for r in table_rows if r['genomad_hit'])
    dvf_hits     = sum(1 for r in table_rows if r['dvf_hit'])
    vs2_hits     = sum(1 for r in table_rows if r['vs2_hit'])
    print(f'[RESULT] Per-tool hits: geNomad={genomad_hits}  DVF={dvf_hits}  VS2={vs2_hits}',
          file=sys.stderr)
    print(f'[RESULT] Output FASTA : {args.out_fasta}', file=sys.stderr)
    print(f'[RESULT] Output table : {args.out_table}', file=sys.stderr)


if __name__ == '__main__':
    main()
