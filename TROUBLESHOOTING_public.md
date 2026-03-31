# CLAUDE Pipeline — Troubleshooting Log

> **Purpose:** Comprehensive record of problems encountered during pipeline development
> and validation. Use this as a first reference when re-running or porting the pipeline.
>
> **Format:** Each issue includes Symptom, Root Cause, Fix, and Prevention where applicable.

---

## Table of Contents

1. [Conda / Environment Issues](#1-conda--environment-issues)
2. [Hostile (Host Depletion)](#2-hostile-host-depletion)
3. [KrakenUniq Build](#3-krakenuniq-build)
4. [Assembly](#4-assembly)
5. [iPHoP — Database Download](#5-iphop--database-download)
6. [iPHoP — Tool Execution](#6-iphop--tool-execution)
7. [iPHoP — Memory / System Crashes](#7-iphop--memory--system-crashes)
8. [ISEScan (Mobilome)](#8-isescan-mobilome)
9. [IntegronFinder (Mobilome)](#9-integronfinder-mobilome)
10. [mobileOG-db (Mobilome)](#10-mobileog-db-mobilome)
11. [PHASTEST (Prophage)](#11-phastest-prophage)
12. [VirSorter2](#12-virsorter2)
13. [PropagAtE](#13-propagate)
14. [Integration Script](#14-integration-script)
15. [Ubuntu System Issues](#15-ubuntu-system-issues)

---

## 1. Conda / Environment Issues

### 1.1 pip install blocked on Ubuntu 24.04

**Symptom:** `pip install` fails with PEP 668 error ("externally managed environment")

**Root cause:** Ubuntu 24.04 protects the system Python from pip modifications.

**Fix:** Always install into a named conda environment using mamba, or use the
explicit conda env pip binary:
```bash
$HOME/miniconda3/envs/<env_name>/bin/pip install <package>
```

**Prevention:** Never run bare `pip install` outside a conda environment.

---

### 1.2 `conda run` doesn't inherit environment variables

**Symptom:** Commands using `conda run -n env cmd` fail when the command uses env vars
set in the calling shell (e.g., `$DB_DIR`, `$PROJECT`).

**Root cause:** `conda run` spawns a clean shell; parent env vars are not inherited.

**Fix:** Pass variables explicitly inside a `bash -c "..."` wrapper:
```bash
conda run -n env bash -c "export MYVAR=/path && command"
```

**Prevention:** Use `bash -c "..."` for all multi-command conda run calls.

---

### 1.3 `conda run` + pipes

**Symptom:** `conda run -n env cmd1 | cmd2` pipes silently fail or output goes
to the wrong place.

**Fix:**
```bash
conda run -n env bash -c "cmd1 | cmd2"
```

---

### 1.4 System Python hijacks subprocess calls inside conda env

**Symptom:** iPHoP fails with `ModuleNotFoundError: No module named 'joblib'`
even though joblib is installed in the conda env.

**Root cause:** iPHoP's internal subprocess calls use bare `python3`, which
resolves to `/usr/bin/python3` (Ubuntu 24.04 system Python 3.12) instead of
the conda env's Python 3.8.

**Fix:** Prepend the conda env bin to PATH before calling iPHoP:
```bash
conda run -n iphop_env bash -c "
export PATH=$HOME/miniconda3/envs/iphop_env/bin:\$PATH
iphop predict ...
"
```

**Prevention:** Always set PATH explicitly when running tools that spawn
subprocesses that need to resolve to the env's Python.

---

### 1.5 `conda run -n claude_pipeline python3` hits system Python 3.12

**Symptom:** NumPy compatibility errors when running integration scripts via
`conda run -n claude_pipeline python3`.

**Root cause:** `python3` resolves to system Python 3.12; the claude_pipeline
env uses Python 3.10.

**Fix:** Use the explicit interpreter path:
```bash
$HOME/miniconda3/envs/claude_pipeline/bin/python3.10 script.py
```

---

## 2. Hostile (Host Depletion)

### 2.1 API changed between hostile 1.x and 2.x

**Symptom:** `hostile fetch --name X --destination /path` fails in hostile 2.x.

**Root cause:** hostile 2.0 changed the CLI completely.

**Fix (hostile 2.x):**
```bash
hostile index fetch --name human-t2t-hla-argos985 --bowtie2
```
Indexes are stored in `~/.local/share/hostile/` automatically. No `--destination` flag.

---

## 3. KrakenUniq Build

### 3.1 jellyfish version conflict

**Symptom:** `krakenuniq-build` fails with jellyfish error.

**Root cause:** KrakenUniq requires jellyfish 1.x; jellyfish 2.x (default in
many conda channels) is incompatible.

**Fix:**
```bash
mamba install -c bioconda jellyfish=1.1.12
export JELLYFISH_BIN=''   # must be empty string, not unset
krakenuniq-build --db ...
```

**Prevention:** Always install jellyfish=1.1.12 explicitly when setting up KrakenUniq.

---

### 3.2 Thermal shutdown during KrakenUniq build on laptop

**Symptom:** System powers off mid-build.

**Root cause:** KrakenUniq build is CPU-intensive and sustains 100% load across
all threads. On laptops without active thermal management this triggers automatic
shutdown.

**Fix:**
```bash
nice -n 19 krakenuniq-build --db ... --threads 1
```
Run at lowest scheduling priority with a single thread to stay within thermal limits.

---

## 4. Assembly

### 4.1 MEGAHIT fails if output directory exists

**Symptom:** MEGAHIT exits immediately with an error about existing output directory.

**Root cause:** MEGAHIT does not overwrite existing output directories.

**Fix:** Always remove before running:
```bash
rm -rf /path/to/megahit_out && megahit -1 R1.fq -2 R2.fq -o megahit_out
```

---

### 4.2 All assemblers produce 0 contigs from test reads

**Symptom:** A small test dataset produces no assembled contigs.

**Root cause:** This is expected. De novo assembly requires a minimum number of
reads to build a de Bruijn graph with sufficient coverage. A read set of only a
few hundred pairs is far below the assembly threshold for any microbial genome.

**Prevention:** Document this in the SOP; use a real metagenome dataset for
assembly validation.

---

## 5. iPHoP — Database Download

### 5.1 Jun25 database: wrong download flag

**Symptom:** Running `iphop download --full` does nothing useful or runs verification only.

**Root cause:** `--full` in iPHoP 1.4.2 is a verification flag, not a download flag.

**Fix:**
```bash
iphop download --db_dir databases/iphop/ --split --no_prompt
```

---

### 5.2 Jun25 database: disk space exhausted during chunk merge

**Symptom:** `cat iphop_db_chunk_*.tar.gz > iPHoP.latest_rw.tar.gz` runs out
of disk space mid-merge. Partially written archive wastes ~175 GB.

**Root cause:** 28 chunks × ~10 GB = 278 GB compressed. Uncompressed BLAST
binary databases expand to ~350–500 GB. The disk must have ≥500 GB free before
starting the merge.

**Fix:** Delete the partial archive to reclaim space:
```bash
rm databases/iphop/iPHoP.latest_rw.tar.gz
```
Then either move to a server with ≥500 GB free, or use the test database for
tool validation.

---

### 5.3 Test database: wrong name causes 404

**Symptom:** iPHoP test database download fails on md5 check.

**Root cause:** The test database is named `iPHoP_db_rw_1.4_for-test` (with `_rw_`
and `_1.4_`). Earlier attempts used `iPHoP_db_for-test` which returns 404 on the NERSC server.

**Fix:** Use the exact name:
```bash
iphop download --db_dir databases/iphop_test/ --split --no_prompt \
    --db_name iPHoP_db_rw_1.4_for-test
```

---

## 6. iPHoP — Tool Execution

### 6.1 `--no_overlap` flag does not exist

**Symptom:** `iphop: error: unrecognized arguments: --no_overlap`

**Fix:** Remove the flag. It does not exist in iPHoP 1.4.2.

---

### 6.2 PHP step fails: `ModuleNotFoundError: No module named 'joblib'`

**See §1.4 above.** Fix: `export PATH=$HOME/miniconda3/envs/iphop_env/bin:$PATH`

---

### 6.3 Skipping RaFAH (R Random Forest step) due to memory constraints

**Symptom:** RaFAH's R process loads a 16 GB model and the system runs out of
memory (OOM) or freezes completely.

**Root cause:** The RaFAH Random Forest model requires ~16.4 GB just to load,
plus additional working memory during computation. Total peak RAM usage exceeds
30 GB. On a machine with less than 30 GB free RAM, the process will be killed
by the OS memory manager.

**Fix (laptop workaround):** Pre-create a stub `rafahparsed.csv` before running
iPHoP. iPHoP checks for this file and skips RaFAH if it exists:

```bash
mkdir -p pipeline_run/10_prophages/sample/iphop/Wdir
echo 'Virus,Host_genus,RaFAH_score,RaFAH_rank,Translated_genus' \
    > pipeline_run/10_prophages/sample/iphop/Wdir/rafahparsed.csv

conda run -n iphop_env bash -c "
export PATH=$HOME/miniconda3/envs/iphop_env/bin:\$PATH
iphop predict \
    --fa_file pipeline_run/08_viral_contigs/sample/consensus_viral_contigs.fasta \
    --db_dir databases/iphop_test/Test_db_rw_v1.4 \
    --out_dir pipeline_run/10_prophages/sample/iphop \
    -t 8
"
```

iPHoP will log a warning about RaFAH having no data but will complete
successfully using 5/6 tools (blast genomes, CRISPR spacers, WIsH, VHM, PHP).

**For production (server):** RaFAH should run normally. Requires ≥64 GB RAM.
Do not use the stub on a properly-resourced machine.

**Wdir checkpoint recovery:** If the system crashes after steps 1-5 are
complete (blastparsed.csv, crisprparsed.csv, wishparsed.csv, vhmparsed.csv,
phpparsed.csv all exist in Wdir), simply add the stub and re-run. iPHoP
will skip all completed steps via its own file-existence checks.

---

### 6.4 iPHoP output: understanding the warnings

When running with the RaFAH stub, iPHoP prints:
```
!#!#!#! WARNING
RaFAH results were empty. This may be ok, but is still unusual...
We skipped blast because there was no useable hits
Note that we tried to use RaFAH, but the file .../rafahparsed.csv had no data
```
These warnings are expected when using the stub or with a small test database.
Output files are still valid.

---

### 6.5 systemd-oomd kills GNOME Terminal during RaFAH

**Symptom:** The terminal window closes completely when RaFAH starts loading
its R model, even with sufficient apparent free memory.

**Root cause:** Ubuntu 24.04 runs `systemd-oomd`, a userspace OOM daemon that
is more aggressive than the kernel OOM killer. It monitors memory pressure and
kills entire cgroups (including the terminal) when pressure exceeds 50% for
more than 20 seconds. Critically, stopping `systemd-oomd.service` is not
sufficient because `systemd-oomd.socket` automatically restarts the daemon.

**Fix:** Stop both the service and the socket:
```bash
sudo systemctl stop systemd-oomd.socket systemd-oomd
```
Then expand swap before running:
```bash
sudo fallocate -l 13G /swapfile2
sudo chmod 600 /swapfile2
sudo mkswap /swapfile2
sudo swapon /swapfile2
```
Remember to re-enable oomd after the run:
```bash
sudo systemctl start systemd-oomd.socket systemd-oomd
```

---

## 7. iPHoP — Memory / System Crashes

### 7.1 Ubuntu unattended-upgrades auto-reboot

**Symptom:** Machine reboots without warning, killing long-running background
jobs. Occurs 10-15 minutes after boot.

**Root cause:** Ubuntu's `unattended-upgrades` service applies kernel updates
on a daily timer and reboots automatically if `Automatic-Reboot` is enabled
(the default on Ubuntu 24.04).

**Diagnosis:**
```bash
journalctl -b -1 --no-pager | grep -E "PackageKit|unattended|logind"
# Look for: PackageKit daemon start -> "The system will reboot now!"
```

**Fix:** Disable auto-reboot permanently:
```bash
sudo sed -i \
    's|^//Unattended-Upgrade::Automatic-Reboot "false";|Unattended-Upgrade::Automatic-Reboot "false";|' \
    /etc/apt/apt.conf.d/50unattended-upgrades

grep "Automatic-Reboot" /etc/apt/apt.conf.d/50unattended-upgrades
# Should show (uncommented): Unattended-Upgrade::Automatic-Reboot "false";
```

---

### 7.2 System OOM freeze from RaFAH R model

**Symptom:** Machine becomes completely unresponsive; requires hard power-off.
Signs: kernel OOM kills GNOME shell, terminal, and browser in rapid succession.

**Root cause:** R process loads 16.4 GB model then requires additional working
memory during Random Forest computation. With 30 GB RAM + swap, the process
can attempt to use ~42 GB, saturating both and causing the kernel to freeze.

**Diagnosis:**
```bash
journalctl -b -1 --no-pager | grep "oom-kill"
```

**Fix:** Use the RaFAH stub workaround (see §6.3) on machines with less than
64 GB RAM. The RaFAH step cannot run on a 30 GB RAM laptop regardless of
swap size.

---

### 7.3 Adding swap space

If additional swap is needed as a buffer:
```bash
sudo swapoff /swap.img
sudo dd if=/dev/zero of=/swap.img bs=1G count=20 status=progress
sudo mkswap /swap.img
sudo swapon /swap.img
free -h
```
Note that expanded swap delays but does not prevent OOM if RaFAH is run on a
machine with insufficient RAM. The fix is ≥64 GB physical RAM.

---

## 8. ISEScan (Mobilome)

### 8.1 ISEScan not on PATH in conda env

**Symptom:** `isescan.py` not found even inside claude_pipeline env.

**Root cause:** ISEScan installs as a Python script, not a wrapper binary;
it must be called via the env's Python explicitly.

**Fix:**
```bash
LD_LIBRARY_PATH=$HOME/miniconda3/envs/claude_pipeline/lib:$LD_LIBRARY_PATH \
$HOME/miniconda3/envs/claude_pipeline/bin/python3.10 \
    $(which isescan.py) --seqfile input.fasta --output outdir --nthread 4
```

---

### 8.2 ISEScan produces no output files on zero-result runs

**Symptom:** ISEScan completes successfully but the output directory has no
`.sum`, `.tsv`, or `.gff` files.

**Root cause:** ISEScan only writes output files when IS elements are found.
Zero IS elements means no output files — this is not an error.

**Fix:** Check the log for `ISEScan ends at` to confirm completion:
```bash
grep "ISEScan ends at" isescan.log
```

---

## 9. IntegronFinder (Mobilome)

### 9.1 `--no-func-annot` flag does not exist

**Symptom:** `integron_finder: error: unrecognized arguments: --no-func-annot`

**Fix:** The correct flag for linear/metagenome contigs is `--linear`:
```bash
integron_finder --linear --local-max input.fasta
```

---

## 10. mobileOG-db (Mobilome)

### 10.1 Server timeouts (504)

**Symptom:** wget or browser returns 504 Gateway Timeout when accessing
the mobileOG-db server.

**Root cause:** External server outages — not a pipeline issue.

**Fix:** Retry when the server recovers. Run DIAMOND against the local database:
```bash
diamond blastp \
    --query pipeline_run/12_mobilome/sample/mobileog/all_bin_proteins.faa \
    --db databases/mobileOG-db \
    --out pipeline_run/12_mobilome/sample/mobileog/mobileog_hits.tsv \
    --outfmt 6 --id 90 --query-cover 80 --threads 8
```

---

## 11. PHASTEST (Prophage)

### 11.1 Wrong form field name

**Symptom:** PHASTEST API submission fails or returns an error.

**Root cause:** The form field is `sequence_text`, not `submission[seq_data]`
(which was documented in some older guides).

**Fix:** Use `sequence_text` in the POST request.

---

### 11.2 Batch limit

**Root cause:** PHASTEST accepts max 10 sequences per batch submission.

**Fix:** Split large inputs into batches of 10 sequences or fewer.

---

### 11.3 HPC cluster unavailability

**Symptom:** PHASTEST returns cluster errors and jobs do not run.

**Root cause:** PHASTEST uses an HPC backend that experiences periodic outages.

**Fix:** Batch IDs persist on the server — retry later:
```bash
python3 pipeline_run/scripts/phastest_submit.py \
    --batch_ids <id1>,<id2>,<id3> \
    --out_dir pipeline_run/10_prophages/sample/phastest_results
```

---

## 12. VirSorter2

### 12.1 Runtime impractical on laptop (24h+ for Viruses HMM)

**Symptom:** VirSorter2 runs indefinitely on the Viruses HMM scan step.

**Root cause:** The combined HMM database (8.3 GB) requires very long search
times at low thread counts. On a server with 32+ threads this step completes
in 2-4 hours.

**Decision:** VirSorter2 was excluded from the consensus call on the laptop
validation run; a 2-tool consensus (geNomad + DVF) was used instead. Production
runs should include VS2 as the third tool.

---

### 12.2 VirSorter2 patch required

**Symptom:** VirSorter2 fails on installation or early execution with a
Python compatibility error.

**Fix:** VirSorter2 requires a patch for the current bioconda build.
See `virsorter2_env` setup in phase_7_sop.md for the exact patch applied.

---

## 13. PropagAtE

### 13.1 Installation path

PropagAtE is not available on PyPI and must be installed from source:

```bash
git clone https://github.com/AnantharamanLab/PropagAtE.git \
    pipeline_run/scripts/PropagAtE/
$HOME/miniconda3/envs/claude_pipeline/bin/pip install \
    pipeline_run/scripts/PropagAtE/
```

---

### 13.2 MEGAHIT contig header parsing bug

**Symptom:** PropagAtE reports 0 proviruses found even though geNomad predicted
proviruses in MEGAHIT contigs.

**Root cause:** MEGAHIT contig headers include extra metadata after the contig ID
(`flag=0 multi=... len=...`). PropagAtE parses the full header string as the
contig name, causing a mismatch when looking up contigs in the lengths dict.

**Fix:** Write clean FASTAs with only the contig ID as the header before
passing to PropagAtE:
```bash
seqkit replace -p " .+" -r "" megahit_contigs.fasta > megahit_contigs_clean.fasta
```

---

## 14. Integration Script

### 14.1 QC checkpoint files with no STATUS line

**Symptom:** `integrate_pipeline.py` cannot find STATUS in QC1 and QC2 checkpoint
files, which use box-drawing borders rather than a `STATUS:` line.

**Fix:** Added `_IMPLICIT_STATUS` dict at module level in `integrate_pipeline.py`:
```python
_IMPLICIT_STATUS = {
    'QC1_NA12878_chr22.txt': 'PASS (WARN: chimeric reads only)',
    'QC2_NA12878_chr22.txt': 'PASS',
}
```
`extract_qc_status()` falls back to this dict when no STATUS line is found.

---

### 14.2 Pharokka functions TSV pivot format

**Symptom:** Integration script cannot parse Pharokka functional output.

**Root cause:** `pharokka_cds_functions.tsv` uses a long-format table with
`Description`, `Count`, and `contig` columns — not one row per contig.

**Fix:** Pivot the table by contig in the integration script:
```python
def load_pharokka_functions(path) -> dict:
    # pivot Description/Count/contig into per-contig dict
```

---

## 15. Ubuntu System Issues

### 15.1 Key diagnostic commands

```bash
# List all recent boots (check for unexpected reboots)
journalctl --list-boots --no-pager | tail -10

# Find reboot cause in previous boot
journalctl -b -1 --no-pager | grep -E \
    "(PackageKit|unattended|logind|oom-kill|Killed process|thermal)" | tail -30

# Check current memory and swap
free -h && swapon --show

# Check uptime (confirm no silent reboot)
uptime
```

### 15.2 Production server checklist

Before running memory-intensive steps (iPHoP RaFAH, CheckM2, metaSPAdes) on
any machine:

* ≥64 GB RAM (iPHoP RaFAH minimum)
* ≥500 GB free disk (iPHoP Jun25 DB extraction)
* Auto-reboot disabled (`Unattended-Upgrade::Automatic-Reboot "false"`)
* Swap set to ≥32 GB as a safety buffer
* Run long jobs in `screen` or `tmux` so reboots do not kill sessions
* Use `nohup` and a log file for all background jobs

---

*CLAUDE Pipeline — Troubleshooting Log*
