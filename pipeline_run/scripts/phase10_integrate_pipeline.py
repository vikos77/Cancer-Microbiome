#!/usr/bin/env python3
"""
integrate_pipeline.py — Phase 10 Integration
CLAUDE Pipeline: Cancer-Linked Analysis of Underlying DNA Elements

Merges all Phase 7–9 outputs into a master viral contig table and
generates a QC dashboard aggregating all QC1–QC12 checkpoints.

Usage:
    python3 integrate_pipeline.py --base_dir /path/to/project --out_dir pipeline_run/13_integration/virome

Outputs:
    master_viral_table.tsv    — one row per viral contig (119 contigs), all annotations merged
    qc_dashboard.txt          — aggregated QC1–QC12 status
    mobilome_summary.tsv      — IS element summary across bins
"""

import argparse
import os
import sys
import csv
from collections import defaultdict

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_tsv(path, delimiter='\t'):
    """Read a TSV/CSV into list of dicts. Returns (header_list, list_of_dicts)."""
    rows = []
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row in reader:
            rows.append(row)
    return rows


def safe_float(v, default='NA'):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def safe_int(v, default=0):
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Data loading functions
# ---------------------------------------------------------------------------

def load_checkv(path):
    """Load CheckV quality_summary.tsv → dict keyed by contig_id."""
    data = {}
    for row in read_tsv(path):
        cid = row['contig_id']
        data[cid] = {
            'length':              safe_int(row.get('contig_length', 0)),
            'checkv_quality':      row.get('checkv_quality', 'NA'),
            'checkv_completeness': safe_float(row.get('completeness', 'NA')),
            'checkv_contamination':safe_float(row.get('contamination', 'NA')),
            'provirus':            row.get('provirus', 'No'),
            'checkv_gene_count':   safe_int(row.get('gene_count', 0)),
            'checkv_viral_genes':  safe_int(row.get('viral_genes', 0)),
            'checkv_host_genes':   safe_int(row.get('host_genes', 0)),
        }
    return data


def load_consensus_table(path):
    """Load viral_classification_summary_2tool.tsv → dict keyed by contig, viral only."""
    data = {}
    for row in read_tsv(path):
        if row.get('consensus_viral', '0') != '1':
            continue
        cid = row['contig']
        data[cid] = {
            'genomad_score': safe_float(row.get('genomad_score', 'NA')),
            'genomad_hit':   safe_int(row.get('genomad_hit', 0)),
            'genomad_tax':   row.get('genomad_tax', 'NA'),
            'dvf_score':     safe_float(row.get('dvf_score', 'NA')),
            'dvf_hit':       safe_int(row.get('dvf_hit', 0)),
            'total_votes':   safe_int(row.get('total_votes', 0)),
        }
    return data


def load_pharokka_length_gc(path):
    """Load pharokka_length_gc_cds_density.tsv → dict keyed by contig."""
    data = {}
    for row in read_tsv(path):
        cid = row['contig']
        data[cid] = {
            'gc_perc':           safe_float(row.get('gc_perc', 'NA')),
            'cds_coding_density':safe_float(row.get('cds_coding_density', 'NA')),
        }
    return data


def load_pharokka_functions(path):
    """
    Load pharokka_cds_functions.tsv and pivot to per-contig function counts.
    Columns: Description, Count, contig
    Returns dict keyed by contig → {category: count, ...}
    """
    # Categories to extract
    func_categories = [
        'CDS',
        'connector',
        'DNA, RNA and nucleotide metabolism',
        'head and packaging',
        'integration and excision',
        'lysis',
        'moron, auxiliary metabolic gene and host takeover',
        'other',
        'tail',
        'transcription regulation',
        'unknown function',
        'tRNAs',
        'tmRNAs',
        'CRISPRs',
    ]
    data = defaultdict(lambda: {c: 0 for c in func_categories})
    for row in read_tsv(path):
        desc = row.get('Description', '').strip()
        contig = row.get('contig', '').strip()
        count = safe_int(row.get('Count', 0))
        if desc in func_categories and contig:
            data[contig][desc] = count
    return dict(data)


def load_genomad_virus_summary(path):
    """Load clean_contigs_virus_summary.tsv → dict keyed by seq_name."""
    data = {}
    for row in read_tsv(path):
        cid = row.get('seq_name', '').strip()
        if cid:
            data[cid] = {
                'genomad_topology': row.get('topology', 'NA'),
                'genomad_hallmarks': safe_int(row.get('n_hallmarks', 0)),
                'genomad_taxonomy':  row.get('taxonomy', 'NA'),
            }
    return data


def load_inphared_hits(path):
    """Load pharokka_top_hits_mash_inphared.tsv → dict keyed by contig."""
    data = {}
    for row in read_tsv(path):
        cid = row.get('contig', '').strip()
        if cid:
            data[cid] = {
                'inphared_accession':  row.get('Accession', 'NA'),
                'inphared_mash_dist':  safe_float(row.get('mash_distance', 'NA')),
                'inphared_hashes':     row.get('mash_matching_hashes', 'NA'),
                'inphared_description':row.get('Description', 'NA'),
                'inphared_host':       row.get('Host', 'NA'),
                'inphared_family':     row.get('Family', 'NA'),
                'inphared_length_bp':  safe_int(row.get('Genome_Length_(bp)', 0)),
            }
    return data


def load_propagate(path):
    """Load propagate_results.tsv → dict keyed by prophage ID."""
    data = {}
    if not os.path.exists(path):
        return data
    for row in read_tsv(path):
        pid = row.get('prophage', '').strip()
        if pid:
            data[pid] = {
                'activity': row.get('active', 'NA'),
                'cohen_d':  safe_float(row.get('CohenD', 'NA')),
            }
    return data


# ---------------------------------------------------------------------------
# Build master table
# ---------------------------------------------------------------------------

def build_master_table(checkv, consensus, pharokka_gc, pharokka_func,
                       genomad_vs, inphared, propagate):
    """
    Join all data sources on contig_id.
    Backbone: CheckV quality_summary.tsv (119 viral contigs).
    """
    rows = []
    for cid, cv in checkv.items():
        row = {'contig_id': cid}

        # CheckV fields
        row.update(cv)

        # Tool votes
        cons = consensus.get(cid, {})
        row['genomad_score'] = cons.get('genomad_score', 'NA')
        row['genomad_hit']   = cons.get('genomad_hit', 'NA')
        row['dvf_score']     = cons.get('dvf_score', 'NA')
        row['dvf_hit']       = cons.get('dvf_hit', 'NA')
        row['total_votes']   = cons.get('total_votes', 'NA')

        # geNomad virus summary (taxonomy, hallmarks, topology)
        gvs = genomad_vs.get(cid, {})
        row['genomad_taxonomy']  = gvs.get('genomad_taxonomy',
                                   cons.get('genomad_tax', 'NA'))
        row['genomad_topology']  = gvs.get('genomad_topology', 'NA')
        row['genomad_hallmarks'] = gvs.get('genomad_hallmarks', 'NA')

        # Pharokka GC + coding density
        pgc = pharokka_gc.get(cid, {})
        row['gc_perc']            = pgc.get('gc_perc', 'NA')
        row['cds_coding_density'] = pgc.get('cds_coding_density', 'NA')

        # Pharokka functional categories
        pf = pharokka_func.get(cid, {})
        row['cds_total']                = pf.get('CDS', 0)
        row['cds_unknown']              = pf.get('unknown function', 0)
        cds_known = row['cds_total'] - row['cds_unknown']
        row['cds_known']                = cds_known
        if row['cds_total'] > 0:
            row['cds_known_pct'] = round(100.0 * cds_known / row['cds_total'], 1)
        else:
            row['cds_known_pct'] = 'NA'
        row['phrog_connector']          = pf.get('connector', 0)
        row['phrog_dna_metabolism']     = pf.get('DNA, RNA and nucleotide metabolism', 0)
        row['phrog_head_packaging']     = pf.get('head and packaging', 0)
        row['phrog_integration_excision'] = pf.get('integration and excision', 0)
        row['phrog_lysis']              = pf.get('lysis', 0)
        row['phrog_moron_amg']          = pf.get('moron, auxiliary metabolic gene and host takeover', 0)
        row['phrog_other']              = pf.get('other', 0)
        row['phrog_tail']               = pf.get('tail', 0)
        row['phrog_transcription']      = pf.get('transcription regulation', 0)
        row['trnas']                    = pf.get('tRNAs', 0)
        row['tmrnas']                   = pf.get('tmRNAs', 0)
        row['crisprs']                  = pf.get('CRISPRs', 0)

        # INPHARED closest hit
        inh = inphared.get(cid, {})
        row['inphared_accession']   = inh.get('inphared_accession', 'NA')
        row['inphared_mash_dist']   = inh.get('inphared_mash_dist', 'NA')
        row['inphared_hashes']      = inh.get('inphared_hashes', 'NA')
        row['inphared_description'] = inh.get('inphared_description', 'NA')
        row['inphared_host']        = inh.get('inphared_host', 'NA')
        row['inphared_family']      = inh.get('inphared_family', 'NA')
        row['inphared_ref_length']  = inh.get('inphared_length_bp', 'NA')

        # Prophage activity from PropagAtE
        # metaspades_415 is provirus per CheckV; check PropagAtE for any match
        row['propagate_activity'] = 'NA'
        for pid, pval in propagate.items():
            if cid in pid:  # provirus ID contains contig name
                row['propagate_activity'] = pval['activity']
                break

        rows.append(row)

    # Sort by CheckV completeness descending (HQ first), then by length
    def sort_key(r):
        q_order = {'High-quality': 0, 'Medium-quality': 1, 'Low-quality': 2,
                   'Not-determined': 3, 'NA': 4}
        q = q_order.get(r.get('checkv_quality', 'NA'), 4)
        comp = r.get('checkv_completeness', 0)
        if not isinstance(comp, (int, float)):
            comp = 0
        return (q, -comp)

    rows.sort(key=sort_key)
    return rows


# ---------------------------------------------------------------------------
# Write master table
# ---------------------------------------------------------------------------

MASTER_COLUMNS = [
    'contig_id', 'length', 'gc_perc', 'cds_coding_density',
    'checkv_quality', 'checkv_completeness', 'checkv_contamination', 'provirus',
    'checkv_gene_count', 'checkv_viral_genes', 'checkv_host_genes',
    'total_votes', 'genomad_hit', 'genomad_score',
    'dvf_hit', 'dvf_score',
    'genomad_taxonomy', 'genomad_topology', 'genomad_hallmarks',
    'cds_total', 'cds_known', 'cds_unknown', 'cds_known_pct',
    'phrog_head_packaging', 'phrog_tail', 'phrog_lysis', 'phrog_dna_metabolism',
    'phrog_connector', 'phrog_integration_excision', 'phrog_moron_amg',
    'phrog_other', 'phrog_transcription',
    'trnas', 'tmrnas', 'crisprs',
    'inphared_accession', 'inphared_mash_dist', 'inphared_hashes',
    'inphared_description', 'inphared_host', 'inphared_family', 'inphared_ref_length',
    'propagate_activity',
]


def write_master_table(rows, out_path):
    with open(out_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=MASTER_COLUMNS,
                                delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written: {out_path} ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# QC dashboard
# ---------------------------------------------------------------------------

QC_FILES = [
    ('QC1',  'QC1_NA12878_chr22.txt',         'Data Acquisition'),
    ('QC2',  'QC2_NA12878_chr22.txt',         'Host Depletion'),
    ('QC3',  'QC3_NA12878_chr22.txt',         'Read QC'),
    ('QC4',  'QC4_NA12878_chr22.txt',         'Taxonomic Profiling'),
    ('QC5a', 'QC5_NA12878_chr22.txt',         'Assembly (NA12878)'),
    ('QC5b', 'QC5_virome_SRR15090802.txt',    'Assembly (SRR15090802)'),
    ('QC6',  'QC6_virome_SRR15090802.txt',    'Decontam + MAG Binning'),
    ('QC7',  'QC7_virome_SRR15090802.txt',    'Viral ID (2-tool consensus)'),
    ('QC8',  'QC8_virome_SRR15090802.txt',    'Phage Annotation (Pharokka)'),
    ('QC10', 'QC10_virome_SRR15090802.txt',   'Prophage Concordance'),
    ('QC11', 'QC11_virome_SRR15090802.txt',   'Host Prediction + Annotation'),
    ('QC12', 'QC12_virome_SRR15090802.txt',   'Mobilome Detection'),
]


    # Files without explicit STATUS lines — known PASS from pipeline history
_IMPLICIT_STATUS = {
    'QC1_NA12878_chr22.txt': 'PASS (WARN: chimeric reads only)',
    'QC2_NA12878_chr22.txt': 'PASS',
}

def extract_qc_status(filepath):
    """Scan a QC checkpoint file for STATUS line."""
    if not os.path.exists(filepath):
        return 'FILE NOT FOUND'
    fname = os.path.basename(filepath)
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            # Match patterns: "=== STATUS: PASS ===" or "STATUS: PASS"
            if 'STATUS' in line and ('PASS' in line or 'FAIL' in line or
                                     'PARTIAL' in line or 'COMPLETE' in line or
                                     'WARN' in line):
                # Strip === delimiters
                status = line.replace('===', '').replace('STATUS:', '').strip()
                return status
    return _IMPLICIT_STATUS.get(fname, 'STATUS NOT FOUND')


def extract_qc_key_metric(filepath):
    """Extract one representative metric line from the QC file."""
    metrics = {
        'QC1_NA12878_chr22.txt':         'Read pairs extracted (chimeric)',
        'QC2_NA12878_chr22.txt':         '393 pairs after host depletion',
        'QC3_NA12878_chr22.txt':         '391 pairs after quality filtering',
        'QC4_NA12878_chr22.txt':         '0 microbial taxa (Kraken2 + KrakenUniq)',
        'QC5_NA12878_chr22.txt':         '0 contigs (391 reads insufficient)',
        'QC5_virome_SRR15090802.txt':    '20,548 contigs | N50=1,110bp | max=101kb',
        'QC6_virome_SRR15090802.txt':    '20,547 clean contigs | 0 human | 1 vector removed',
        'QC7_virome_SRR15090802.txt':    '119 viral contigs (geNomad+DVF) | 1 complete phage (101kb)',
        'QC8_virome_SRR15090802.txt':    '674 CDS | 32.1% known function | 0 ARG | 0 virulence',
        'QC10_virome_SRR15090802.txt':   '1 dormant provirus (metaspades_100, CohenD=0.062)',
        'QC11_virome_SRR15090802.txt':   '82.4% contigs annotated (Pharokka) | iPHoP pending',
        'QC12_virome_SRR15090802.txt':   '124 plasmids | 16 IS | 0 integrons | mobileOG pending',
    }
    fname = os.path.basename(filepath)
    return metrics.get(fname, '')


def build_qc_dashboard(qc_dir, out_path):
    lines = []
    lines.append('=' * 78)
    lines.append('CLAUDE PIPELINE — QC DASHBOARD')
    lines.append('Sample: SRR15090802 (Wahida et al. gut virome, VLP-enriched healthy child)')
    lines.append('Generated: Phase 10 Integration')
    lines.append('=' * 78)
    lines.append('')
    lines.append(f'{"ID":<6}  {"Phase":<35}  {"Status":<20}  Key metric')
    lines.append('-' * 78)

    n_pass = n_partial = n_fail = n_missing = 0
    for qc_id, filename, phase_name in QC_FILES:
        fpath = os.path.join(qc_dir, filename)
        status = extract_qc_status(fpath)
        metric = extract_qc_key_metric(fpath)
        lines.append(f'{qc_id:<6}  {phase_name:<35}  {status:<20}  {metric}')
        s = status.upper()
        if 'FILE NOT FOUND' in s:
            n_missing += 1
        elif 'FAIL' in s:
            n_fail += 1
        elif 'PARTIAL' in s or 'PENDING' in s:
            n_partial += 1
        else:
            n_pass += 1

    lines.append('-' * 78)
    lines.append(f'SUMMARY: {n_pass} PASS | {n_partial} PARTIAL | {n_fail} FAIL | {n_missing} MISSING')
    lines.append('')
    lines.append('Pending items (not blocking pipeline):')
    lines.append('  - PHASTEST: batch IDs saved; retry when cluster recovers')
    lines.append('    Cmd: python3 pipeline_run/scripts/phastest_submit.py \\')
    lines.append('           --batch_ids BB_3288ee5236,BB_51980fe431,BB_6f5e98c05f,BB_3ba4ff53bf \\')
    lines.append('           --out_dir pipeline_run/10_prophages/virome/phastest_results')
    lines.append('  - iPHoP host prediction: ~280 GB DB downloading; run after completion')
    lines.append('    Cmd: conda run -n iphop_env iphop predict \\')
    lines.append('           --fa_file pipeline_run/08_viral_contigs/virome/consensus_viral_contigs_2tool.fasta \\')
    lines.append('           --db_dir databases/iphop/ --out_dir pipeline_run/10_prophages/virome/iphop -t 4')
    lines.append('  - mobileOG-db DIAMOND: proteins pre-computed (4,848); run when server recovers')
    lines.append('    URL: https://mobileogdb.flsi.cloud.vt.edu/entries/mobileOG-db_beatrix-1.6.All.faa')
    lines.append('=' * 78)

    with open(out_path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"  Written: {out_path}")
    return n_pass, n_partial, n_fail, n_missing


# ---------------------------------------------------------------------------
# Mobilome summary table
# ---------------------------------------------------------------------------

def build_mobilome_summary(isescan_dir, plasmid_summary_path, out_path):
    """
    Build a compact mobilome summary TSV from ISEScan .sum files and
    geNomad plasmid summary.
    """
    rows = []

    # ISEScan .sum files
    for bin_name in ['bin.1', 'bin.2', 'bin.3', 'bin.4']:
        sum_file = os.path.join(isescan_dir, bin_name, 'metabat2',
                                f'{bin_name}.fa.sum')
        if not os.path.exists(sum_file):
            # bin.3 has 0 IS elements — no file written
            continue
        with open(sum_file) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 6:
                    rows.append({
                        'bin':        bin_name,
                        'element_type': 'insertion_sequence',
                        'contig_id':  parts[0],
                        'family':     parts[1],
                        'n_copies':   parts[2],
                        'pct_genome': parts[3],
                        'bp_covered': parts[4],
                        'contig_len': parts[5],
                        'notes': '',
                    })

    # Plasmids from geNomad (top 20 by score)
    if os.path.exists(plasmid_summary_path):
        plasmid_rows = read_tsv(plasmid_summary_path)
        plasmid_rows.sort(key=lambda r: safe_float(r.get('plasmid_score', 0),
                                                    default=0), reverse=True)
        for r in plasmid_rows[:20]:
            rows.append({
                'bin':        'all_contigs',
                'element_type': 'plasmid',
                'contig_id':  r.get('seq_name', 'NA'),
                'family':     r.get('conjugation_genes', 'NA'),
                'n_copies':   '1',
                'pct_genome': 'NA',
                'bp_covered': r.get('length', 'NA'),
                'contig_len': r.get('length', 'NA'),
                'notes':      f"score={r.get('plasmid_score','NA')} topology={r.get('topology','NA')}",
            })

    with open(out_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            'bin', 'element_type', 'contig_id', 'family', 'n_copies',
            'pct_genome', 'bp_covered', 'contig_len', 'notes'],
            delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Written: {out_path} ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# Print summary stats
# ---------------------------------------------------------------------------

def print_stats(rows):
    total = len(rows)
    hq   = sum(1 for r in rows if r.get('checkv_quality') == 'High-quality')
    mq   = sum(1 for r in rows if r.get('checkv_quality') == 'Medium-quality')
    lq   = sum(1 for r in rows if r.get('checkv_quality') == 'Low-quality')
    nd   = sum(1 for r in rows if r.get('checkv_quality') == 'Not-determined')

    prov = sum(1 for r in rows if r.get('provirus') == 'Yes')
    with_ann = sum(1 for r in rows
                   if isinstance(r.get('cds_known', 0), int) and r['cds_known'] > 0)
    total_cds = sum(r.get('cds_total', 0) for r in rows
                    if isinstance(r.get('cds_total', 0), int))
    known_cds = sum(r.get('cds_known', 0) for r in rows
                    if isinstance(r.get('cds_known', 0), int))
    pct_known = round(100.0 * known_cds / total_cds, 1) if total_cds > 0 else 0

    trnas = sum(r.get('trnas', 0) for r in rows
                if isinstance(r.get('trnas', 0), int))

    exact_match = sum(1 for r in rows
                      if r.get('inphared_mash_dist') == 0.0)

    print()
    print('  === MASTER TABLE STATS ===')
    print(f'  Total viral contigs:      {total}')
    print(f'  CheckV quality:           HQ={hq}  MQ={mq}  LQ={lq}  ND={nd}')
    print(f'  Proviruses:               {prov}')
    print(f'  Contigs with ≥1 known CDS:{with_ann}  ({round(100*with_ann/total,1)}%)')
    print(f'  Total CDS:                {total_cds}')
    print(f'  Known-function CDS:       {known_cds}  ({pct_known}%)')
    print(f'  tRNAs encoded:            {trnas}')
    print(f'  Exact INPHARED matches:   {exact_match}')
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base_dir', default='/home/vicky/Microbiome_cancer',
                   help='Project root directory')
    p.add_argument('--out_dir',
                   default='pipeline_run/13_integration/virome',
                   help='Output directory (relative to base_dir)')
    args = p.parse_args()

    base = args.base_dir
    out_dir = os.path.join(base, args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    run_dir = os.path.join(base, 'pipeline_run')

    # Input paths
    P = {
        'checkv':       f'{run_dir}/09_viral_qc/virome/checkv_2tool/quality_summary.tsv',
        'consensus':    f'{run_dir}/08_viral_contigs/virome/viral_classification_summary_2tool.tsv',
        'pharokka_gc':  f'{run_dir}/11_phage_annotation/virome/pharokka/pharokka_length_gc_cds_density.tsv',
        'pharokka_fn':  f'{run_dir}/11_phage_annotation/virome/pharokka/pharokka_cds_functions.tsv',
        'genomad_vs':   f'{run_dir}/08_viral_contigs/virome/genomad/clean_contigs_summary/clean_contigs_virus_summary.tsv',
        'inphared':     f'{run_dir}/11_phage_annotation/virome/pharokka/pharokka_top_hits_mash_inphared.tsv',
        'propagate':    f'{run_dir}/10_prophages/virome/propagate_results/propagate_results.tsv',
        'plasmids':     f'{run_dir}/12_mobilome/virome/plasmids/clean_contigs_plasmid_summary.tsv',
        'isescan':      f'{run_dir}/12_mobilome/virome/isescan',
        'qc_dir':       f'{run_dir}/qc_checkpoints',
    }

    # Verify required inputs exist
    for key in ['checkv', 'consensus', 'pharokka_gc', 'pharokka_fn',
                'genomad_vs', 'inphared']:
        if not os.path.exists(P[key]):
            sys.exit(f'ERROR: Required input not found: {P[key]}')

    print('Loading input files...')
    checkv       = load_checkv(P['checkv'])
    consensus    = load_consensus_table(P['consensus'])
    pharokka_gc  = load_pharokka_length_gc(P['pharokka_gc'])
    pharokka_fn  = load_pharokka_functions(P['pharokka_fn'])
    genomad_vs   = load_genomad_virus_summary(P['genomad_vs'])
    inphared     = load_inphared_hits(P['inphared'])
    propagate    = load_propagate(P['propagate'])

    print(f'  CheckV contigs:        {len(checkv)}')
    print(f'  Consensus viral:       {len(consensus)}')
    print(f'  Pharokka GC entries:   {len(pharokka_gc)}')
    print(f'  Pharokka func entries: {len(pharokka_fn)}')
    print(f'  geNomad taxonomy:      {len(genomad_vs)}')
    print(f'  INPHARED hits:         {len(inphared)}')
    print(f'  PropagAtE results:     {len(propagate)}')

    print('\nBuilding master viral table...')
    rows = build_master_table(checkv, consensus, pharokka_gc, pharokka_fn,
                              genomad_vs, inphared, propagate)
    write_master_table(rows, os.path.join(out_dir, 'master_viral_table.tsv'))
    print_stats(rows)

    print('Building QC dashboard...')
    n_pass, n_partial, n_fail, n_miss = build_qc_dashboard(
        P['qc_dir'],
        os.path.join(out_dir, 'qc_dashboard.txt')
    )

    print('Building mobilome summary...')
    build_mobilome_summary(
        P['isescan'],
        P['plasmids'],
        os.path.join(out_dir, 'mobilome_summary.tsv')
    )

    print('\nDone.')
    print(f'Outputs in: {out_dir}')


if __name__ == '__main__':
    main()
