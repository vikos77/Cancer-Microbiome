#!/usr/bin/env python3
"""
integrate_wirbel_crc.py — Phase 10 Integration (Wirbel CRC Production Run)
CLAUDE Pipeline: Cancer-Linked Analysis of Underlying DNA Elements

Merges Phase 7–11 outputs for CRC_ERR2726414 and control_ERR2726517
into master tables and a QC dashboard.

Usage:
    python3 integrate_wirbel_crc.py [--base_dir /path/to/project]
                                    [--out_dir pipeline_run/13_integration/Wirbel_CRC]

Outputs (per sample):
    master_viral_table_{sample}.tsv   — one row per viral contig
    defence_systems_summary_{sample}.tsv — per-bin defence counts (DF+PADLOC+CCF)
    mobilome_summary_{sample}.tsv     — per-bin plasmid/IS/mobileOG counts
    qc_dashboard_Wirbel_CRC.txt       — QC2–QC14 status for both samples
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

SAMPLES = ['CRC_ERR2726414', 'control_ERR2726517']
SAMPLE_NBINS = {'CRC_ERR2726414': 30, 'control_ERR2726517': 92}

# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def read_tsv(path, delimiter='\t'):
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
        return int(float(v))
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# Phase 7 — Viral ID
# ---------------------------------------------------------------------------

def load_checkv(path):
    """Load CheckV quality_summary.tsv → dict keyed by contig_id."""
    data = {}
    for row in read_tsv(path):
        cid = row['contig_id']
        data[cid] = {
            'length':               safe_int(row.get('contig_length', 0)),
            'checkv_quality':       row.get('checkv_quality', 'NA'),
            'checkv_completeness':  safe_float(row.get('completeness', 'NA')),
            'checkv_contamination': safe_float(row.get('contamination', 'NA')),
            'provirus':             row.get('provirus', 'No'),
            'checkv_gene_count':    safe_int(row.get('gene_count', 0)),
            'checkv_viral_genes':   safe_int(row.get('viral_genes', 0)),
            'checkv_host_genes':    safe_int(row.get('host_genes', 0)),
        }
    return data


def load_consensus(path):
    """Load consensus_viral_table.tsv → dict keyed by contig (viral only)."""
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


def load_genomad_virus_summary(path):
    """Load clean_contigs_virus_summary.tsv → dict keyed by seq_name."""
    data = {}
    if not os.path.exists(path):
        return data
    for row in read_tsv(path):
        cid = row.get('seq_name', '').strip()
        if cid:
            data[cid] = {
                'genomad_topology':  row.get('topology', 'NA'),
                'genomad_hallmarks': safe_int(row.get('n_hallmarks', 0)),
                'genomad_taxonomy':  row.get('taxonomy', 'NA'),
            }
    return data


# ---------------------------------------------------------------------------
# Phase 8 — Phage Annotation
# ---------------------------------------------------------------------------

def load_pharokka_length_gc(path):
    data = {}
    for row in read_tsv(path):
        cid = row['contig']
        data[cid] = {
            'gc_perc':            safe_float(row.get('gc_perc', 'NA')),
            'cds_coding_density': safe_float(row.get('cds_coding_density', 'NA')),
        }
    return data


def load_pharokka_functions(path):
    FUNC_CATS = [
        'CDS', 'connector', 'DNA, RNA and nucleotide metabolism',
        'head and packaging', 'integration and excision', 'lysis',
        'moron, auxiliary metabolic gene and host takeover', 'other',
        'tail', 'transcription regulation', 'unknown function',
        'tRNAs', 'tmRNAs', 'CRISPRs',
    ]
    data = defaultdict(lambda: {c: 0 for c in FUNC_CATS})
    for row in read_tsv(path):
        desc   = row.get('Description', '').strip()
        contig = row.get('contig', '').strip()
        count  = safe_int(row.get('Count', 0))
        if desc in FUNC_CATS and contig:
            data[contig][desc] = count
    return dict(data)


def load_inphared_hits(path):
    data = {}
    if not os.path.exists(path):
        return data
    for row in read_tsv(path):
        cid = row.get('contig', '').strip()
        if cid:
            data[cid] = {
                'inphared_accession':   row.get('Accession', 'NA'),
                'inphared_mash_dist':   safe_float(row.get('mash_distance', 'NA')),
                'inphared_hashes':      row.get('mash_matching_hashes', 'NA'),
                'inphared_description': row.get('Description', 'NA'),
                'inphared_host':        row.get('Host', 'NA'),
                'inphared_family':      row.get('Family', 'NA'),
                'inphared_ref_length':  safe_int(row.get('Genome_Length_(bp)', 0)),
            }
    return data


def load_iphop(path):
    """Load iPHoP Host_prediction_to_genus_m90.csv → dict keyed by virus name."""
    data = {}
    if not os.path.exists(path):
        return data
    for row in read_tsv(path, delimiter=','):
        vid = row.get('Virus', '').strip()
        if vid:
            host = row.get('Host genus', 'NA')
            # Shorten full lineage to genus only
            if ';' in host:
                host = host.split(';')[-1].strip().lstrip('g__')
            data[vid] = {
                'iphop_host_genus':  host,
                'iphop_confidence':  safe_float(row.get('Confidence score', 'NA')),
                'iphop_methods':     row.get('List of methods', 'NA'),
            }
    return data


def load_propagate(path):
    """Load PropagAtE results. prophage column = 'contig|provirus_start_end'."""
    data = {}
    if not os.path.exists(path):
        return data
    for row in read_tsv(path):
        pid = row.get('prophage', '').strip()
        if pid:
            # Extract just the contig name (before the | separator)
            contig_key = pid.split('|')[0]
            data[contig_key] = {
                'propagate_activity': row.get('active', 'NA'),
                'propagate_cohen_d':  safe_float(row.get('CohenD', 'NA')),
                'propagate_ratio':    safe_float(row.get('prophage-host_ratio', 'NA')),
                'propagate_full_id':  pid,
            }
    return data


# ---------------------------------------------------------------------------
# Build master viral table
# ---------------------------------------------------------------------------

MASTER_COLUMNS = [
    'contig_id', 'length', 'gc_perc', 'cds_coding_density',
    'checkv_quality', 'checkv_completeness', 'checkv_contamination', 'provirus',
    'checkv_gene_count', 'checkv_viral_genes', 'checkv_host_genes',
    'total_votes', 'genomad_hit', 'genomad_score', 'dvf_hit', 'dvf_score',
    'genomad_taxonomy', 'genomad_topology', 'genomad_hallmarks',
    'cds_total', 'cds_known', 'cds_unknown', 'cds_known_pct',
    'phrog_head_packaging', 'phrog_tail', 'phrog_lysis', 'phrog_dna_metabolism',
    'phrog_connector', 'phrog_integration_excision', 'phrog_moron_amg',
    'phrog_other', 'phrog_transcription', 'trnas', 'tmrnas', 'crisprs',
    'inphared_accession', 'inphared_mash_dist', 'inphared_hashes',
    'inphared_description', 'inphared_host', 'inphared_family', 'inphared_ref_length',
    'iphop_host_genus', 'iphop_confidence', 'iphop_methods',
    'propagate_activity', 'propagate_cohen_d', 'propagate_ratio',
]


def build_master_table(checkv, consensus, pharokka_gc, pharokka_fn,
                       genomad_vs, inphared, iphop, propagate):
    rows = []
    for cid, cv in checkv.items():
        row = {'contig_id': cid}
        row.update(cv)

        # Consensus votes
        cons = consensus.get(cid, {})
        row['genomad_score'] = cons.get('genomad_score', 'NA')
        row['genomad_hit']   = cons.get('genomad_hit', 'NA')
        row['dvf_score']     = cons.get('dvf_score', 'NA')
        row['dvf_hit']       = cons.get('dvf_hit', 'NA')
        row['total_votes']   = cons.get('total_votes', 'NA')

        # geNomad summary
        gvs = genomad_vs.get(cid, {})
        row['genomad_taxonomy']  = gvs.get('genomad_taxonomy',
                                   cons.get('genomad_tax', 'NA'))
        row['genomad_topology']  = gvs.get('genomad_topology', 'NA')
        row['genomad_hallmarks'] = gvs.get('genomad_hallmarks', 'NA')

        # Pharokka
        pgc = pharokka_gc.get(cid, {})
        row['gc_perc']            = pgc.get('gc_perc', 'NA')
        row['cds_coding_density'] = pgc.get('cds_coding_density', 'NA')

        pf = pharokka_fn.get(cid, {})
        row['cds_total']    = pf.get('CDS', 0)
        row['cds_unknown']  = pf.get('unknown function', 0)
        cds_known = row['cds_total'] - row['cds_unknown']
        row['cds_known']    = cds_known
        row['cds_known_pct'] = (round(100.0 * cds_known / row['cds_total'], 1)
                                if row['cds_total'] > 0 else 'NA')
        row['phrog_connector']            = pf.get('connector', 0)
        row['phrog_dna_metabolism']       = pf.get('DNA, RNA and nucleotide metabolism', 0)
        row['phrog_head_packaging']       = pf.get('head and packaging', 0)
        row['phrog_integration_excision'] = pf.get('integration and excision', 0)
        row['phrog_lysis']                = pf.get('lysis', 0)
        row['phrog_moron_amg']            = pf.get('moron, auxiliary metabolic gene and host takeover', 0)
        row['phrog_other']                = pf.get('other', 0)
        row['phrog_tail']                 = pf.get('tail', 0)
        row['phrog_transcription']        = pf.get('transcription regulation', 0)
        row['trnas']   = pf.get('tRNAs', 0)
        row['tmrnas']  = pf.get('tmRNAs', 0)
        row['crisprs'] = pf.get('CRISPRs', 0)

        # INPHARED
        inh = inphared.get(cid, {})
        for k in ('inphared_accession', 'inphared_mash_dist', 'inphared_hashes',
                  'inphared_description', 'inphared_host', 'inphared_family',
                  'inphared_ref_length'):
            row[k] = inh.get(k, 'NA')

        # iPHoP — key is the virus name from CheckV/consensus (may differ from iPHoP key)
        iph = iphop.get(cid, {})
        # Also try stripping suffixes added by CheckV (e.g. _1, _2 for provirus regions)
        if not iph:
            base = cid.rsplit('_', 1)[0]
            iph = iphop.get(base, {})
        row['iphop_host_genus']  = iph.get('iphop_host_genus', 'NA')
        row['iphop_confidence']  = iph.get('iphop_confidence', 'NA')
        row['iphop_methods']     = iph.get('iphop_methods', 'NA')

        # PropagAtE — keyed by contig name (before | separator)
        pval = propagate.get(cid, {})
        row['propagate_activity'] = pval.get('propagate_activity', 'NA')
        row['propagate_cohen_d']  = pval.get('propagate_cohen_d', 'NA')
        row['propagate_ratio']    = pval.get('propagate_ratio', 'NA')

        rows.append(row)

    # Sort: HQ first, then by completeness desc
    q_order = {'High-quality': 0, 'Medium-quality': 1, 'Low-quality': 2,
               'Not-determined': 3, 'NA': 4}
    rows.sort(key=lambda r: (
        q_order.get(r.get('checkv_quality', 'NA'), 4),
        -(r.get('checkv_completeness', 0) if isinstance(r.get('checkv_completeness'), float) else 0)
    ))
    return rows


def write_master_table(rows, path):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=MASTER_COLUMNS, delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    print(f"  Written: {path} ({len(rows)} rows)")


def print_viral_stats(rows, sample):
    total = len(rows)
    hq  = sum(1 for r in rows if r.get('checkv_quality') == 'High-quality')
    mq  = sum(1 for r in rows if r.get('checkv_quality') == 'Medium-quality')
    lq  = sum(1 for r in rows if r.get('checkv_quality') == 'Low-quality')
    nd  = sum(1 for r in rows if r.get('checkv_quality') == 'Not-determined')
    prov = sum(1 for r in rows if r.get('provirus') == 'Yes')
    total_cds = sum(r.get('cds_total', 0) for r in rows
                    if isinstance(r.get('cds_total', 0), int))
    known_cds = sum(r.get('cds_known', 0) for r in rows
                    if isinstance(r.get('cds_known', 0), int))
    pct_known = round(100.0 * known_cds / total_cds, 1) if total_cds else 0
    active = sum(1 for r in rows if r.get('propagate_activity') == 'active')
    iphop_pred = sum(1 for r in rows if r.get('iphop_host_genus', 'NA') != 'NA')
    print(f"\n  [{sample}] Viral contigs: {total}  HQ={hq} MQ={mq} LQ={lq} ND={nd}")
    print(f"    Proviruses: {prov}  Active prophages: {active}")
    print(f"    CDS: {total_cds} total, {known_cds} known ({pct_known}%)")
    print(f"    iPHoP host predictions: {iphop_pred}")


# ---------------------------------------------------------------------------
# Phase 9 — Mobilome
# ---------------------------------------------------------------------------

def load_isescan_per_bin(isescan_dir, nbins):
    """Return dict bin → {total_IS, families: {family: count}}."""
    data = {}
    for i in range(1, nbins + 1):
        bin_name = f'bin.{i}'
        sum_file = os.path.join(isescan_dir, bin_name, 'metabat2',
                                f'{bin_name}.fa.sum')
        total = 0
        families = defaultdict(int)
        if os.path.exists(sum_file):
            with open(sum_file) as fh:
                for line in fh:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) >= 3:
                        fam = parts[1]
                        n   = safe_int(parts[2])
                        families[fam] += n
                        total += n
        data[bin_name] = {'total_IS': total, 'families': dict(families)}
    return data


def load_plasmid_counts_per_bin(plasmid_summary_path):
    """Count plasmids per bin from geNomad plasmid_summary.tsv.
    seq_name contains the contig name; map back to bin via contig→bin if needed.
    For now, return total count only (plasmids are from clean_contigs, not per-bin).
    """
    if not os.path.exists(plasmid_summary_path):
        return 0, []
    rows = read_tsv(plasmid_summary_path)
    return len(rows), rows


def load_mobileog_per_bin(hits_path, nbins):
    """Count mobileOG hits per bin.
    Format: query_id (=bin.N_proteins_{contig}_{N}), hit_id|gene|...|class|...|
    """
    data = {f'bin.{i}': defaultdict(int) for i in range(1, nbins + 1)}
    if not os.path.exists(hits_path):
        return data
    # Keep best hit per query (first occurrence, already sorted by score in DIAMOND)
    seen = set()
    with open(hits_path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            if not row:
                continue
            query = row[0]
            # Extract bin name from query: format is {bin.N}_{protein_id}
            # e.g. "megahit_k119_79903_4" — need to map back
            # The FAA files are named bin.N_proteins.faa so queries come from those
            # But the query name = contig_N from prodigal, not "bin.N_..."
            # We need to use a bin-to-contig mapping or use the per-bin FAA indexing
            # Since we can't easily map without the index, just count total
            if query not in seen:
                seen.add(query)
                hit = row[1] if len(row) > 1 else ''
                # Extract class from hit: format mobileOG_ID|gene|prot|class|...
                parts = hit.split('|')
                mob_class = parts[3] if len(parts) > 3 else 'unknown'
                data['all'] = data.get('all', defaultdict(int))
                data['all'][mob_class] += 1
    return data


def count_mobileog_total(hits_path):
    """Count unique best-hit mobileOG queries."""
    if not os.path.exists(hits_path):
        return 0, {}
    seen = set()
    classes = defaultdict(int)
    with open(hits_path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            if not row:
                continue
            query = row[0]
            if query not in seen:
                seen.add(query)
                hit = row[1] if len(row) > 1 else ''
                parts = hit.split('|')
                mob_class = parts[3] if len(parts) > 3 else 'unknown'
                classes[mob_class] += 1
    return len(seen), dict(classes)


def build_mobilome_summary(sample, run_dir, nbins, out_path):
    isescan_dir   = f'{run_dir}/12_mobilome/{sample}/isescan'
    plasmid_path  = f'{run_dir}/12_mobilome/{sample}/plasmids/clean_contigs_plasmid_summary.tsv'
    mobileog_path = f'{run_dir}/12_mobilome/{sample}/mobileog/mobileog_all_hits.tsv'

    isescan   = load_isescan_per_bin(isescan_dir, nbins)
    n_plasmid, plasmid_rows = load_plasmid_counts_per_bin(plasmid_path)
    n_mobileog, mobileog_classes = count_mobileog_total(mobileog_path)

    rows = []
    for i in range(1, nbins + 1):
        bin_name = f'bin.{i}'
        is_data  = isescan[bin_name]
        top_fam  = (sorted(is_data['families'].items(), key=lambda x: -x[1])[0][0]
                    if is_data['families'] else 'none')
        rows.append({
            'bin':         bin_name,
            'total_IS':    is_data['total_IS'],
            'IS_top_fam':  top_fam,
            'IS_families': ';'.join(f"{k}:{v}" for k, v in
                           sorted(is_data['families'].items(), key=lambda x: -x[1])),
        })

    # Add summary rows for plasmids and mobileOG
    rows.append({
        'bin': 'ALL_CONTIGS',
        'total_IS': '',
        'IS_top_fam': f'plasmids={n_plasmid}',
        'IS_families': f'mobileOG_hits={n_mobileog}; classes={mobileog_classes}',
    })

    with open(out_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['bin', 'total_IS', 'IS_top_fam', 'IS_families'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    total_IS = sum(isescan[f'bin.{i}']['total_IS'] for i in range(1, nbins + 1))
    print(f"  Written: {out_path} | IS={total_IS}  Plasmids={n_plasmid}  mobileOG={n_mobileog}")
    return total_IS, n_plasmid, n_mobileog


# ---------------------------------------------------------------------------
# Phase 11 — Defence Systems
# ---------------------------------------------------------------------------

def load_defensefinder_per_bin(df_dir, nbins):
    """Return {bin: {total, types: {type: count}}}."""
    data = {}
    for i in range(1, nbins + 1):
        bin_name = f'bin.{i}'
        bin_dir  = os.path.join(df_dir, bin_name)
        total = 0
        types = defaultdict(int)
        if os.path.isdir(bin_dir):
            for fname in os.listdir(bin_dir):
                if fname.endswith('_systems.tsv'):
                    fpath = os.path.join(bin_dir, fname)
                    for row in read_tsv(fpath):
                        sys_type = row.get('type', '').strip()
                        if sys_type:
                            types[sys_type] += 1
                            total += 1
        data[bin_name] = {'total': total, 'types': dict(types)}
    return data


def load_padloc_per_bin(padloc_dir, nbins):
    """Return {bin: {total, types: {type: count}}}."""
    data = {}
    for i in range(1, nbins + 1):
        bin_name = f'bin.{i}'
        bin_dir  = os.path.join(padloc_dir, bin_name)
        total = 0
        types = defaultdict(int)
        seen  = set()
        if os.path.isdir(bin_dir):
            for fname in os.listdir(bin_dir):
                if fname.endswith('_padloc.csv'):
                    fpath = os.path.join(bin_dir, fname)
                    for row in read_tsv(fpath, delimiter=','):
                        sys_num  = row.get('system.number', '')
                        seqid    = row.get('seqid', '')
                        sys_type = row.get('system', '').strip()
                        key = (sys_num, seqid)
                        if key not in seen and sys_type:
                            seen.add(key)
                            types[sys_type] += 1
                            total += 1
        data[bin_name] = {'total': total, 'types': dict(types)}
    return data


def load_ccf_per_bin(ccf_dir, nbins):
    """Return {bin: {total_arrays}}."""
    data = {}
    for i in range(1, nbins + 1):
        bin_name = f'bin.{i}'
        tsv = os.path.join(ccf_dir, bin_name, 'TSV', 'Crisprs_REPORT.tsv')
        n_arrays = 0
        if os.path.exists(tsv):
            with open(tsv) as fh:
                for j, line in enumerate(fh):
                    if j > 0 and line.strip():
                        n_arrays += 1
        data[bin_name] = {'total_arrays': n_arrays}
    return data


def build_defence_summary(sample, run_dir, nbins, out_path):
    df_dir  = f'{run_dir}/14_defence_systems/{sample}/defensefinder'
    pl_dir  = f'{run_dir}/14_defence_systems/{sample}/padloc'
    ccf_dir = f'{run_dir}/14_defence_systems/{sample}/crisprcasfinder'

    df_data  = load_defensefinder_per_bin(df_dir, nbins)
    pl_data  = load_padloc_per_bin(pl_dir, nbins)
    ccf_data = load_ccf_per_bin(ccf_dir, nbins)

    rows = []
    for i in range(1, nbins + 1):
        bin_name = f'bin.{i}'
        df  = df_data[bin_name]
        pl  = pl_data[bin_name]
        ccf = ccf_data[bin_name]

        df_top  = (sorted(df['types'].items(), key=lambda x: -x[1])[0][0]
                   if df['types'] else 'none')
        pl_top  = (sorted(pl['types'].items(), key=lambda x: -x[1])[0][0]
                   if pl['types'] else 'none')

        rows.append({
            'bin':                bin_name,
            'df_total':           df['total'],
            'df_top_type':        df_top,
            'df_types':           ';'.join(f"{k}:{v}" for k, v in
                                  sorted(df['types'].items(), key=lambda x: -x[1])[:5]),
            'padloc_total':       pl['total'],
            'padloc_top_type':    pl_top,
            'padloc_types':       ';'.join(f"{k}:{v}" for k, v in
                                  sorted(pl['types'].items(), key=lambda x: -x[1])[:5]),
            'crispr_arrays':      ccf['total_arrays'],
        })

    total_df  = sum(df_data[f'bin.{i}']['total'] for i in range(1, nbins + 1))
    total_pl  = sum(pl_data[f'bin.{i}']['total'] for i in range(1, nbins + 1))
    total_ccf = sum(ccf_data[f'bin.{i}']['total_arrays'] for i in range(1, nbins + 1))

    with open(out_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=[
            'bin', 'df_total', 'df_top_type', 'df_types',
            'padloc_total', 'padloc_top_type', 'padloc_types', 'crispr_arrays'],
            delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    df_types_n  = len(set(t for d in df_data.values() for t in d['types']))
    pl_types_n  = len(set(t for d in pl_data.values() for t in d['types']))
    bins_pos_df = sum(1 for d in df_data.values() if d['total'] > 0)
    bins_pos_pl = sum(1 for d in pl_data.values() if d['total'] > 0)
    bins_pos_ccf = sum(1 for d in ccf_data.values() if d['total_arrays'] > 0)
    print(f"  Written: {out_path}")
    print(f"    DF:  {total_df} systems / {df_types_n} types / {bins_pos_df}/{nbins} bins")
    print(f"    PADLOC: {total_pl} systems / {pl_types_n} types / {bins_pos_pl}/{nbins} bins")
    print(f"    CCF: {total_ccf} arrays / {bins_pos_ccf}/{nbins} bins")
    return total_df, total_pl, total_ccf


# ---------------------------------------------------------------------------
# QC dashboard
# ---------------------------------------------------------------------------

QC_FILES_WIRBEL = [
    # (qc_id, filename, phase_name)
    ('QC2-CRC',  'QC2_CRC_ERR2726414.txt',     'Host Depletion (CRC)'),
    ('QC2-CTL',  'QC2_control_ERR2726517.txt',  'Host Depletion (Control)'),
    ('QC3-CRC',  'QC3_CRC_ERR2726414.txt',      'Read QC (CRC)'),
    ('QC3-CTL',  'QC3_control_ERR2726517.txt',  'Read QC (Control)'),
    ('QC4-CRC',  'QC4_CRC_ERR2726414.txt',      'Tax Profiling (CRC)'),
    ('QC4-CTL',  'QC4_control_ERR2726517.txt',  'Tax Profiling (Control)'),
    ('QC5-CRC',  'QC5_CRC_ERR2726414.txt',      'Assembly (CRC)'),
    ('QC5-CTL',  'QC5_control_ERR2726517.txt',  'Assembly (Control)'),
    ('QC6-CRC',  'QC6_CRC_ERR2726414.txt',      'Decontam+Binning (CRC)'),
    ('QC6-CTL',  'QC6_control_ERR2726517.txt',  'Decontam+Binning (Control)'),
    ('QC7-CRC',  'QC7_CRC_ERR2726414.txt',      'Viral ID (CRC)'),
    ('QC7-CTL',  'QC7_control_ERR2726517.txt',  'Viral ID (Control)'),
    ('QC8-CRC',  'QC8_CRC_ERR2726414.txt',      'Phage Annotation (CRC)'),
    ('QC8-CTL',  'QC8_control_ERR2726517.txt',  'Phage Annotation (Control)'),
    ('QC10-CRC', 'QC10_CRC_ERR2726414.txt',     'Prophage Activity (CRC)'),
    ('QC10-CTL', 'QC10_control_ERR2726517.txt', 'Prophage Activity (Control)'),
    ('QC11-CRC', 'QC11_CRC_ERR2726414.txt',     'Host Prediction (CRC)'),
    ('QC11-CTL', 'QC11_control_ERR2726517.txt', 'Host Prediction (Control)'),
    ('QC12-CRC', 'QC12_CRC_ERR2726414.txt',     'Mobilome (CRC)'),
    ('QC12-CTL', 'QC12_control_ERR2726517.txt', 'Mobilome (Control)'),
    ('QC14',     'QC14_defence_systems.txt',    'Defence Systems (both)'),
]


def extract_qc_status(filepath):
    if not os.path.exists(filepath):
        return 'FILE NOT FOUND'
    with open(filepath) as fh:
        for line in fh:
            s = line.strip()
            # Match lines with STATUS/VERDICT/OVERALL + outcome word
            if any(kw in s for kw in ('STATUS', 'VERDICT', 'OVERALL')) and \
               any(w in s for w in ('PASS', 'FAIL', 'PARTIAL', 'COMPLETE', 'WARN')):
                return (s.replace('===', '')
                          .replace('STATUS:', '').replace('VERDICT:', '')
                          .replace('OVERALL', '').strip(' :'))
    # Fallback: return first Status: line
    with open(filepath) as fh:
        for line in fh:
            s = line.strip()
            if s.lower().startswith('status:'):
                return s[7:].strip()
    return 'NOT FOUND'


def build_qc_dashboard(qc_dir, out_path, crc_stats, ctl_stats):
    lines = []
    lines.append('=' * 80)
    lines.append('CLAUDE PIPELINE — QC DASHBOARD (Wirbel CRC Production Run)')
    lines.append('CRC:     ERR2726414 (Stage III CRC, F.nucleatum 2.30%)')
    lines.append('Control: ERR2726517 (healthy, depth-matched)')
    lines.append('=' * 80)
    lines.append('')
    lines.append(f'{"ID":<10}  {"Phase":<30}  {"Status"}')
    lines.append('-' * 80)

    n_pass = n_partial = n_fail = n_missing = 0
    for qc_id, filename, phase_name in QC_FILES_WIRBEL:
        fpath  = os.path.join(qc_dir, filename)
        status = extract_qc_status(fpath)
        lines.append(f'{qc_id:<10}  {phase_name:<30}  {status}')
        s = status.upper()
        if 'FILE NOT FOUND' in s or 'NOT FOUND' in s:
            n_missing += 1
        elif 'FAIL' in s:
            n_fail += 1
        elif 'PARTIAL' in s or 'PENDING' in s or 'WARN' in s:
            n_partial += 1
        else:
            n_pass += 1

    lines.append('-' * 80)
    lines.append(f'SUMMARY: {n_pass} PASS | {n_partial} PARTIAL | {n_fail} FAIL | {n_missing} MISSING')
    lines.append('')
    lines.append('KEY FINDINGS SUMMARY')
    lines.append('-' * 80)
    for label, stats in [('CRC (ERR2726414)', crc_stats),
                          ('Control (ERR2726517)', ctl_stats)]:
        lines.append(f'{label}:')
        for k, v in stats.items():
            lines.append(f'  {k}: {v}')
        lines.append('')
    lines.append('=' * 80)

    with open(out_path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"  Written: {out_path}  ({n_pass} PASS | {n_partial} PARTIAL | {n_fail} FAIL)")
    return n_pass, n_partial, n_fail, n_missing


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base_dir', default='/home/vicky/Microbiome_cancer')
    p.add_argument('--out_dir',  default='pipeline_run/13_integration/Wirbel_CRC')
    args = p.parse_args()

    base    = args.base_dir
    out_dir = os.path.join(base, args.out_dir)
    run_dir = os.path.join(base, 'pipeline_run')
    qc_dir  = os.path.join(run_dir, 'qc_checkpoints')
    os.makedirs(out_dir, exist_ok=True)

    all_stats = {}

    for sample in SAMPLES:
        nbins = SAMPLE_NBINS[sample]
        print(f'\n{"="*60}')
        print(f'Processing: {sample}  ({nbins} bins)')
        print(f'{"="*60}')

        P = {
            'checkv':      f'{run_dir}/08_viral_contigs/{sample}/checkv/quality_summary.tsv',
            'consensus':   f'{run_dir}/08_viral_contigs/{sample}/consensus_viral_table.tsv',
            'genomad_vs':  f'{run_dir}/08_viral_contigs/{sample}/genomad/clean_contigs_summary/clean_contigs_virus_summary.tsv',
            'pharokka_gc': f'{run_dir}/11_phage_annotation/{sample}/pharokka/pharokka_length_gc_cds_density.tsv',
            'pharokka_fn': f'{run_dir}/11_phage_annotation/{sample}/pharokka/pharokka_cds_functions.tsv',
            'inphared':    f'{run_dir}/11_phage_annotation/{sample}/pharokka/pharokka_top_hits_mash_inphared.tsv',
            'iphop':       f'{run_dir}/10_prophages/{sample}/iphop/Host_prediction_to_genus_m90.csv',
            'propagate':   f'{run_dir}/10_prophages/{sample}/propagate_results/propagate_results.tsv',
        }

        # Verify required files
        for key in ('checkv', 'consensus', 'pharokka_gc', 'pharokka_fn'):
            if not os.path.exists(P[key]):
                sys.exit(f'ERROR: Required input missing: {P[key]}')

        print('Loading viral annotation data...')
        checkv      = load_checkv(P['checkv'])
        consensus   = load_consensus(P['consensus'])
        genomad_vs  = load_genomad_virus_summary(P['genomad_vs'])
        pharokka_gc = load_pharokka_length_gc(P['pharokka_gc'])
        pharokka_fn = load_pharokka_functions(P['pharokka_fn'])
        inphared    = load_inphared_hits(P['inphared'])
        iphop       = load_iphop(P['iphop'])
        propagate   = load_propagate(P['propagate'])

        print(f'  CheckV contigs:  {len(checkv)}')
        print(f'  Consensus viral: {len(consensus)}')
        print(f'  iPHoP hosts:     {len(iphop)}')
        print(f'  PropagAtE:       {len(propagate)}')

        print('Building master viral table...')
        rows = build_master_table(checkv, consensus, pharokka_gc, pharokka_fn,
                                  genomad_vs, inphared, iphop, propagate)
        viral_out = os.path.join(out_dir, f'master_viral_table_{sample}.tsv')
        write_master_table(rows, viral_out)
        print_viral_stats(rows, sample)

        print('Building defence systems summary...')
        def_out = os.path.join(out_dir, f'defence_systems_summary_{sample}.tsv')
        total_df, total_pl, total_ccf = build_defence_summary(
            sample, run_dir, nbins, def_out)

        print('Building mobilome summary...')
        mob_out = os.path.join(out_dir, f'mobilome_summary_{sample}.tsv')
        total_IS, n_plasmid, n_mobileog = build_mobilome_summary(
            sample, run_dir, nbins, mob_out)

        total = len(rows)
        hq  = sum(1 for r in rows if r.get('checkv_quality') == 'High-quality')
        mq  = sum(1 for r in rows if r.get('checkv_quality') == 'Medium-quality')
        active = sum(1 for r in rows if r.get('propagate_activity') == 'active')

        all_stats[sample] = {
            'Viral contigs (CheckV)':   total,
            'HQ/MQ contigs':            f'{hq}/{mq}',
            'Active prophages':         active,
            'DefenseFinder systems':    total_df,
            'PADLOC systems':           total_pl,
            'CRISPR arrays':            total_ccf,
            'Plasmids':                 n_plasmid,
            'IS elements':              total_IS,
            'mobileOG hits':            n_mobileog,
        }

    print('\nBuilding QC dashboard...')
    qc_out = os.path.join(out_dir, 'qc_dashboard_Wirbel_CRC.txt')
    build_qc_dashboard(qc_dir, qc_out,
                       all_stats.get('CRC_ERR2726414', {}),
                       all_stats.get('control_ERR2726517', {}))

    print(f'\nAll outputs in: {out_dir}')
    print('Done.')


if __name__ == '__main__':
    main()
