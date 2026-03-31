#!/usr/bin/env python3
"""
PHASTEST web API submission script for CLAUDE pipeline Phase 8.
Submits FASTA files to phastest.ca, polls for completion, downloads results.

Usage:
    python3 phastest_submit.py --input_dir <dir_with_fasta_files> --out_dir <output_dir>
    python3 phastest_submit.py --batch_ids BB_abc123,BB_def456 --out_dir <output_dir>

Notes:
    - PHASTEST processes sequences as a batch (up to 10 per submission).
    - Batch IDs have format BB_<hex>; individual submission IDs are ZZ_<hex>.
    - Results downloadable from /batches/<batch_id>.zip when complete.
    - Backend cluster outages are common; retry --batch_ids to check completed jobs.

Requirements:
    pip install requests
"""

import os
import sys
import time
import json
import argparse
import requests
from pathlib import Path
from html.parser import HTMLParser
import re


BASE_URL = "https://phastest.ca"
POLL_INTERVAL = 120  # seconds between status checks
MAX_WAIT = 7200      # 2 hours max wait per batch


# ---------------------------------------------------------------------------
# HTML parsing helpers
# ---------------------------------------------------------------------------

class CSRFParser(HTMLParser):
    """Extract CSRF token from meta tag or hidden input."""
    def __init__(self):
        super().__init__()
        self.token = None

    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if tag == 'meta' and d.get('name') == 'csrf-token':
            self.token = d.get('content')
        elif tag == 'input' and d.get('name') == 'authenticity_token':
            self.token = d.get('value')


class BatchTableParser(HTMLParser):
    """Extract individual submission IDs from a batch page table."""
    def __init__(self):
        super().__init__()
        self.submission_ids = []

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            d = dict(attrs)
            href = d.get('href', '')
            m = re.match(r'/submissions/(ZZ_[a-f0-9]+)', href)
            if m:
                self.submission_ids.append(m.group(1))


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def get_csrf_token(session):
    """Fetch CSRF token from the submissions/new page."""
    resp = session.get(f"{BASE_URL}/submissions/new", timeout=30)
    resp.raise_for_status()
    p = CSRFParser()
    p.feed(resp.text)
    return p.token


def submit_fasta(session, fasta_path, csrf_token):
    """
    Submit a FASTA file to PHASTEST via the paste-text form.
    Returns batch_id string (BB_<hex>) or None on failure.

    PHASTEST limits each batch to 10 sequences; additional sequences
    are silently discarded. Submit individual large contigs if needed.
    """
    with open(fasta_path, 'r') as f:
        seq_content = f.read()

    data = {
        'utf8': '\u2713',
        'authenticity_token': csrf_token,
        'submission[category]': 'text',
        'sequence_text': seq_content,
        'submission[bacterial_sensitivity]': 'lite',
    }
    headers = {
        'X-CSRF-Token': csrf_token,
        'Referer': f"{BASE_URL}/submissions/new",
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    }
    resp = session.post(f"{BASE_URL}/submissions", data=data, headers=headers,
                        allow_redirects=True, timeout=120)

    # Successful submission redirects to /batches/BB_<hex>
    if '/batches/' in resp.url:
        batch_id = resp.url.rstrip('/').split('/batches/')[-1]
        return batch_id

    # Fallback: parse from response body
    m = re.search(r'/batches/(BB_[a-f0-9]+)', resp.text)
    if m:
        return m.group(1)

    print(f"  Response status: {resp.status_code}, URL: {resp.url}")
    return None


def get_batch_status(session, batch_id):
    """
    Check status of a PHASTEST batch.
    Returns (overall_status, {sub_id: status, ...}).

    Statuses: 'complete', 'running', 'queued', 'error', 'cluster_error', 'unknown'
    """
    headers = {'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'}
    resp = session.get(f"{BASE_URL}/batches/{batch_id}", timeout=30, headers=headers)

    if resp.status_code != 200 or resp.content[:2] == b'PK':
        # Empty ZIP means results not ready yet
        return 'running', {}

    text = resp.text

    # Extract individual submission IDs
    p = BatchTableParser()
    p.feed(text)
    sub_ids = p.submission_ids

    # Check each submission status
    sub_statuses = {}
    for sid in sub_ids:
        sub_resp = session.get(f"{BASE_URL}/submissions/{sid}", timeout=30, headers=headers)
        sub_text = sub_resp.text
        if 'backend computing cluster' in sub_text:
            sub_statuses[sid] = 'cluster_error'
        elif 'Complete' in sub_text and 'alert-danger' not in sub_text:
            sub_statuses[sid] = 'complete'
        elif 'alert-danger' in sub_text or ('failed' in sub_text.lower() and 'cluster' not in sub_text):
            sub_statuses[sid] = 'error'
        elif 'Queued' in sub_text:
            sub_statuses[sid] = 'queued'
        else:
            sub_statuses[sid] = 'running'

    # Aggregate
    statuses_set = set(sub_statuses.values()) if sub_statuses else {'unknown'}
    if all(s == 'complete' for s in statuses_set):
        overall = 'complete'
    elif 'running' in statuses_set or 'queued' in statuses_set:
        overall = 'running'
    elif 'cluster_error' in statuses_set:
        overall = 'cluster_error'
    elif all(s in ('error', 'cluster_error') for s in statuses_set):
        overall = 'error'
    else:
        overall = 'running'

    return overall, sub_statuses


def download_batch_results(session, batch_id, out_dir):
    """Download completed batch results ZIP from PHASTEST."""
    headers = {'Accept': 'application/zip,*/*;q=0.8'}
    zip_url = f"{BASE_URL}/batches/{batch_id}.zip"
    resp = session.get(zip_url, timeout=120, headers=headers)

    if resp.status_code == 200 and resp.content[:2] == b'PK' and len(resp.content) > 100:
        zip_path = out_dir / f"{batch_id}.zip"
        zip_path.write_bytes(resp.content)
        print(f"  Results ZIP saved: {zip_path} ({len(resp.content):,} bytes)")
        return zip_path
    elif len(resp.content) <= 22:
        print(f"  Batch {batch_id}: ZIP empty (results not yet available)")
    else:
        print(f"  Batch {batch_id}: download failed (status {resp.status_code})")
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Submit FASTAs to PHASTEST web server')
    parser.add_argument('--input_dir', '-i',
                        help='Directory containing FASTA files (.fa, .fasta, .fna)')
    parser.add_argument('--out_dir', '-o', required=True,
                        help='Output directory for results')
    parser.add_argument('--batch_ids', '-b',
                        help='Comma-separated batch IDs (BB_...) to retrieve (skip submission)')
    args = parser.parse_args()

    if not args.input_dir and not args.batch_ids:
        parser.error('Either --input_dir or --batch_ids is required')

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0'
    })

    # --- Mode 2: retrieve already-submitted batches ---
    if args.batch_ids:
        batch_ids = [b.strip() for b in args.batch_ids.split(',')]
        print(f"Checking {len(batch_ids)} batch(es)...")
        for batch_id in batch_ids:
            print(f"\nBatch {batch_id}:")
            status, sub_statuses = get_batch_status(session, batch_id)
            print(f"  Overall status: {status}")
            for sid, ss in sub_statuses.items():
                print(f"    {sid}: {ss}")
            if status == 'complete':
                download_batch_results(session, batch_id, out_dir)
            elif status == 'cluster_error':
                print(f"  PHASTEST cluster unavailable — retry later")
        return

    # --- Mode 1: submit FASTA files ---
    input_dir = Path(args.input_dir)
    fastas = sorted(list(input_dir.glob('*.fa')) +
                    list(input_dir.glob('*.fasta')) +
                    list(input_dir.glob('*.fna')))
    if not fastas:
        print(f"No FASTA files found in {input_dir}")
        sys.exit(1)

    print(f"Found {len(fastas)} FASTA file(s) to submit")
    print("Note: PHASTEST processes up to 10 sequences per batch")

    # Get CSRF token
    print("\nFetching CSRF token...")
    csrf_token = get_csrf_token(session)
    if not csrf_token:
        print("ERROR: Could not retrieve CSRF token. PHASTEST may be unavailable.")
        sys.exit(1)
    print(f"  CSRF token obtained")

    submitted_batches = {}
    for fasta in fastas:
        print(f"\nSubmitting {fasta.name}...")
        batch_id = submit_fasta(session, fasta, csrf_token)
        if batch_id:
            print(f"  Batch ID: {batch_id}")
            print(f"  Batch URL: {BASE_URL}/batches/{batch_id}")
            submitted_batches[fasta.stem] = batch_id
        else:
            print(f"  WARNING: Submission failed for {fasta.name}")
        time.sleep(5)  # polite delay between submissions

        # Refresh CSRF token for next submission
        csrf_token = get_csrf_token(session)

    # Save batch IDs
    batches_file = out_dir / 'phastest_batch_ids.json'
    with open(batches_file, 'w') as f:
        json.dump(submitted_batches, f, indent=2)
    print(f"\nBatch IDs saved to: {batches_file}")
    if submitted_batches:
        ids_str = ','.join(submitted_batches.values())
        print(f"Retrieve later with: --batch_ids {ids_str}")

    # Poll for completion
    if not submitted_batches:
        print("No batches to poll.")
        return

    print(f"\nPolling for results (every {POLL_INTERVAL}s, max {MAX_WAIT//60}min per batch)...")
    for name, batch_id in submitted_batches.items():
        elapsed = 0
        while elapsed < MAX_WAIT:
            status, sub_statuses = get_batch_status(session, batch_id)
            n_done = sum(1 for s in sub_statuses.values() if s == 'complete')
            n_total = len(sub_statuses)
            print(f"  [{elapsed//60:.0f}m] {name} (batch {batch_id}): {status} ({n_done}/{n_total} complete)")

            if status == 'complete':
                download_batch_results(session, batch_id, out_dir)
                break
            elif status == 'cluster_error':
                print(f"  PHASTEST cluster unavailable. Save batch ID {batch_id} and retry later.")
                break
            elif status == 'error':
                print(f"  Batch {batch_id} failed. Check {BASE_URL}/batches/{batch_id}")
                break

            time.sleep(POLL_INTERVAL)
            elapsed += POLL_INTERVAL


if __name__ == '__main__':
    main()
