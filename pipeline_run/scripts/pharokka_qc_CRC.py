#!/usr/bin/env python3
"""
Pharokka QC aggregation for CRC ERR2726414.
Reads pharokka_cds_functions.tsv and pharokka_top_hits_mash_inphared.tsv
and reports per-category totals plus annotated contig fraction.
"""
import csv
from collections import defaultdict

BASE = "/home/vicky/Microbiome_cancer/pipeline_run/11_phage_annotation/CRC_ERR2726414/pharokka"

# ---- 1. Aggregate functional categories ----
cat_totals = defaultdict(int)
contigs_with_known = set()
contigs_all = set()

with open(f"{BASE}/pharokka_cds_functions.tsv") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        desc  = row["Description"].strip()
        count = int(row["Count"])
        contig = row["contig"].strip()
        contigs_all.add(contig)
        cat_totals[desc] += count
        # "known function" = any PHROG category except unknown and non-CDS rows
        known_cats = {
            "connector", "DNA, RNA and nucleotide metabolism", "head and packaging",
            "integration and excision", "lysis",
            "moron, auxiliary metabolic gene and host takeover",
            "other", "tail", "transcription regulation"
        }
        if desc in known_cats and count > 0:
            contigs_with_known.add(contig)

total_cds     = cat_totals["CDS"]
total_unknown = cat_totals["unknown function"]
total_known   = total_cds - total_unknown
total_trnas   = cat_totals["tRNAs"]
total_crisprs = cat_totals["CRISPRs"]
total_card    = cat_totals["CARD_AMR_Genes"]
total_vfdb    = cat_totals["VFDB_Virulence_Factors"]

annotated_n   = len(contigs_with_known)
total_contigs = len(contigs_all)
annotated_pct = 100.0 * annotated_n / total_contigs if total_contigs else 0

print("=" * 60)
print("PHAROKKA QC — CRC ERR2726414 (335 contigs)")
print("=" * 60)
print(f"\nGene Prediction (Pyrodigal-gv):")
print(f"  Total CDS predicted : {total_cds}")
print(f"  tRNAs detected      : {total_trnas}")
print(f"  CRISPRs (MinCED)    : {total_crisprs}")
print(f"  Avg CDS per contig  : {total_cds/total_contigs:.1f}")

print(f"\nPHROG Functional Annotation:")
for cat in [
    "DNA, RNA and nucleotide metabolism", "head and packaging", "tail",
    "connector", "lysis", "integration and excision",
    "moron, auxiliary metabolic gene and host takeover",
    "transcription regulation", "other", "unknown function"
]:
    n = cat_totals[cat]
    pct = 100.0 * n / total_cds if total_cds else 0
    print(f"  {cat:<48s} {n:>5d}  ({pct:.1f}%)")

print(f"\n  Known function total  : {total_known}  ({100.0*total_known/total_cds:.1f}%)")
print(f"  Annotated contigs    : {annotated_n}/{total_contigs} ({annotated_pct:.1f}%)")

print(f"\nBiosafety Screens:")
print(f"  CARD AMR genes      : {total_card}")
print(f"  VFDB virulence      : {total_vfdb}")

# ---- 2. Top Mash INPHARED hits ----
print(f"\nTop Mash INPHARED Hits (closest known phages, distance < 0.2):")
hits = []
with open(f"{BASE}/pharokka_top_hits_mash_inphared.tsv") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        acc = row.get("Accession", "").strip()
        if acc and acc != "no_inphared_mash_hit":
            hits.append({
                "contig": row["contig"],
                "acc": acc,
                "dist": row.get("mash_distance", ""),
                "hashes": row.get("mash_matching_hashes", ""),
                "desc": row.get("Description", ""),
                "host": row.get("Host", ""),
            })

if hits:
    for h in sorted(hits, key=lambda x: float(x["dist"]))[:10]:
        print(f"  Contig : {h['contig'][:55]}")
        print(f"  Hit    : {h['acc']}  dist={h['dist']}  hashes={h['hashes']}")
        print(f"  Desc   : {h['desc'][:80]}")
        print(f"  Host   : {h['host']}")
        print()
else:
    print("  (no hits within Mash distance 0.2)")

# ---- 3. CARD detail ----
print("CARD ARG hits (detail):")
with open(f"{BASE}/top_hits_card.tsv") as fh:
    for line in fh:
        print(f"  {line.rstrip()}")

print("\nVFDB hits (detail):")
with open(f"{BASE}/top_hits_vfdb.tsv") as fh:
    for line in fh:
        print(f"  {line.rstrip()}")
