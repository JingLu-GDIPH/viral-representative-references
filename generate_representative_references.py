#!/usr/bin/env python3
"""Generate representative reference sequences for viral mapping.

For each genotype group (defined by user-provided metadata, e.g. RdRp|VP1 combo):
1. QC filter sequences (coverage, missing data, mixed sites)
2. MAFFT align all sequences within group
3. Hierarchical clustering at 99.3% identity
4. Select medoid (real sequence closest to cluster center) per cluster
5. Verify ≥30% distinguishable windows between medoid pairs; merge if below
6. Output FASTA + QC report + clustering report

The medoid approach ensures references stay close to real sample sequences,
unlike ASR/consensus which can drift away from circulating diversity.

Works for ANY virus — provide a FASTA of curated sequences and a metadata
TSV mapping sequence IDs to genotype/serotype labels.

Usage:
  python generate_representative_references.py \
    --input-fasta sequences.fasta \
    --metadata metadata.tsv \
    --group-label "RdRp_VP1" \
    --output reference_panel.fasta \
    --threads 8

Metadata TSV format (tab-delimited):
  accession    genotype    [other columns...]
  AB123456     H1N1pdm     ...
  AB123457     H3N2        ...

The --metadata-genotype-column specifies which column holds the genotype label.
Sequences are grouped by this label, then processed independently.
"""

__version__ = "2.1"
__date__ = "2026-08-06"

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from collections import defaultdict
from itertools import combinations

from Bio import SeqIO
import numpy as np


# ============================================================
# Parameters
# ============================================================
WINDOW_SIZE = 150          # Sliding window size for distinguishability check
WINDOW_STEP = 50           # Step size for sliding window
MIN_MISMATCHES_DIST = 2    # ≥2 mismatches per window = distinguishable
MIN_DIST_WINDOW_PCT = 0.30 # ≥30% distinguishable windows to keep refs separate
CLUSTER_IDENTITY_THRESHOLD = 0.993  # Cluster at 99.3% identity

# QC parameters (Nextclade-style: coverage, missing data, mixed sites)
QC_COV_GOOD = 0.90      # ≥90% ACGT coverage = good
QC_COV_BAD = 0.80       # <80% = bad
QC_MISS_GOOD = 0.02     # <2% N = good
QC_MISS_BAD = 0.10      # ≥10% N = bad
QC_MIXED_THRESHOLD = 10 # >10 non-ACGTN ambiguous bases = bad


# ============================================================
# QC: Sequence quality filtering (3 checks)
# ============================================================
def run_sequence_qc(seq):
    """Run 3 Nextclade-style QC checks on a single ungapped sequence.

    Checks:
      1. Coverage — ACGT fraction of full sequence
      2. Missing data — N fraction
      3. Mixed sites — non-ACGTN ambiguous bases count

    Stop codon and frameshift checks are intentionally NOT included because
    many viruses have overlapping ORFs (e.g. norovirus ORF1/ORF2/ORF3) where
    single-frame translation produces spurious stops. For ORF-level QC, use
    the seqqc skill on individual ORF alignments.

    Returns dict with per-check status and overall qc_status.
    """
    seq = seq.upper()
    n = len(seq)

    # 1. Coverage: ACGT fraction
    acgt = sum(1 for c in seq if c in "ACGT")
    cov = acgt / n if n else 0.0

    # 2. Missing data: N fraction
    missing = seq.count("N") / n if n else 1.0

    # 3. Mixed sites: non-ACGTN ambiguous bases
    mixed = sum(1 for c in seq if c not in "ACGTN-")

    # Classify
    cov_st = "good" if cov >= QC_COV_GOOD else ("mediocre" if cov >= QC_COV_BAD else "bad")
    miss_st = "good" if missing < QC_MISS_GOOD else ("mediocre" if missing < QC_MISS_BAD else "bad")
    mix_st = "good" if mixed <= QC_MIXED_THRESHOLD else "bad"

    order = {"bad": 0, "mediocre": 1, "good": 2}
    checks = [cov_st, miss_st, mix_st]
    overall = min(checks, key=lambda x: order[x])

    return {
        "qc_status": overall,
        "coverage": cov_st,
        "missing": miss_st,
        "mixed": mix_st,
        "cov_pct": round(cov * 100, 2),
        "missing_pct": round(missing * 100, 2),
        "mixed_sites": mixed,
    }


def qc_filter_sequences(seqs, exclude_mediocre=False):
    """Filter sequences by QC status.

    Args:
        seqs: {seq_id: sequence} (ungapped)
        exclude_mediocre: if True, also exclude mediocre sequences

    Returns:
        (passed_seqs, qc_report_rows)
    """
    exclude_status = {"bad", "mediocre"} if exclude_mediocre else {"bad"}
    passed = {}
    report_rows = []

    for sid, seq in seqs.items():
        qc = run_sequence_qc(seq)
        qc["seq_id"] = sid
        report_rows.append(qc)

        if qc["qc_status"] not in exclude_status:
            passed[sid] = seq

    return passed, report_rows


# ============================================================
# Sequence utilities
# ============================================================
def run_mafft(input_fa, output_fa, fast=False):
    args = ["mafft", "--auto", "--quiet"]
    if fast:
        args.append("--retree")
        args.append("1")
    args.append(str(input_fa))
    subprocess.run(args, stdout=open(output_fa, "w"),
                   stderr=subprocess.DEVNULL, check=True)


def aln_identity(a, b):
    """Identity on aligned sequences (excluding gaps and N)."""
    same = eff = 0
    for i in range(min(len(a), len(b))):
        x, y = a[i], b[i]
        if x in "ACGT" and y in "ACGT":
            eff += 1
            if x == y:
                same += 1
    return same / eff if eff else 0.0


def count_distinguishable_windows(aln_a, aln_b, window=WINDOW_SIZE,
                                   step=WINDOW_STEP, min_mm=MIN_MISMATCHES_DIST):
    """Count fraction of windows with ≥min_mm mismatches between two aligned seqs."""
    n_win = n_dist = 0
    for start in range(0, len(aln_a) - window + 1, step):
        wa = aln_a[start:start + window]
        wb = aln_b[start:start + window]
        mm = sum(1 for i in range(window)
                 if wa[i] in "ACGT" and wb[i] in "ACGT" and wa[i] != wb[i])
        n_win += 1
        if mm >= min_mm:
            n_dist += 1
    return n_dist / n_win if n_win else 0.0


def ungapped_seq(aln_seq):
    """Remove gaps from aligned sequence."""
    return aln_seq.replace("-", "")


# ============================================================
# Clustering
# ============================================================
def cluster_sequences(aln_seqs, threshold=CLUSTER_IDENTITY_THRESHOLD):
    """Greedy single-linkage clustering at given identity threshold.

    Returns list of clusters (each = list of sequence IDs).
    """
    ids = list(aln_seqs.keys())
    if len(ids) <= 1:
        return [ids]

    n = len(ids)
    clusters = []
    assigned = set()

    # Sort by decreasing average identity to others (denser sequences first)
    avg_sim = {}
    for i in ids:
        sims = [aln_identity(aln_seqs[i], aln_seqs[j]) for j in ids if j != i]
        avg_sim[i] = np.mean(sims) if sims else 0

    seed_order = sorted(ids, key=lambda x: -avg_sim[x])

    for seed in seed_order:
        if seed in assigned:
            continue
        cluster = [seed]
        assigned.add(seed)
        for other in ids:
            if other in assigned:
                continue
            pid = aln_identity(aln_seqs[seed], aln_seqs[other])
            if pid >= threshold:
                cluster.append(other)
                assigned.add(other)
        clusters.append(cluster)

    return clusters


MIN_MEDOID_LENGTH = 6500  # Minimum ungapped length for a medoid candidate


def select_medoid(aln_seqs, cluster_ids, ungapped_lengths=None,
                  min_length=MIN_MEDOID_LENGTH):
    """Select the medoid: sequence with minimum average distance to all others.

    If ``ungapped_lengths`` is provided, sequences shorter than ``min_length``
    are only chosen as a last resort (when no longer sequence exists in the
    cluster). Among sequences meeting the length threshold, the one with the
    lowest average distance to all cluster members is selected.
    """
    if len(cluster_ids) == 1:
        return cluster_ids[0]

    # Partition cluster into length-qualified and short candidates
    if ungapped_lengths:
        qualified = [sid for sid in cluster_ids
                     if ungapped_lengths.get(sid, 0) >= min_length]
        short = [sid for sid in cluster_ids
                 if ungapped_lengths.get(sid, 0) < min_length]
    else:
        qualified = list(cluster_ids)
        short = []

    # Prefer qualified candidates; fall back to short only if none qualified
    candidates = qualified if qualified else short

    best_id = candidates[0]
    best_avg_dist = float("inf")

    for sid in candidates:
        dists = [1 - aln_identity(aln_seqs[sid], aln_seqs[other])
                 for other in cluster_ids if other != sid]
        avg_dist = np.mean(dists) if dists else 0
        if avg_dist < best_avg_dist:
            best_avg_dist = avg_dist
            best_id = sid

    return best_id


def merge_close_references(refs, aln_seqs, min_dist_pct=MIN_DIST_WINDOW_PCT):
    """Check all medoid pairs; merge if distinguishable windows < threshold."""
    if len(refs) <= 1:
        return refs, []

    merges = []
    changed = True
    while changed:
        changed = False
        min_pct = 1.0
        merge_pair = None
        for i, j in combinations(range(len(refs)), 2):
            pct = count_distinguishable_windows(
                aln_seqs[refs[i]], aln_seqs[refs[j]])
            if pct < min_dist_pct and pct < min_pct:
                min_pct = pct
                merge_pair = (i, j)

        if merge_pair:
            i, j = merge_pair
            merges.append((refs[j], refs[i], min_pct))
            refs.pop(j)
            changed = True

    return refs, merges


# ============================================================
# Metadata parsing
# ============================================================
def parse_metadata(tsv_path, genotype_column="genotype"):
    """Parse metadata TSV to get genotype labels per sequence.

    Returns: {seq_id: genotype_label}
    """
    genotypes = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seq_id = row.get("accession", row.get("Accession", row.get("sequence_id", "")))
            genotype = row.get(genotype_column, "")
            if seq_id and genotype:
                genotypes[seq_id] = genotype
    return genotypes


def infer_genotype_from_header(header, pattern=r"((?:GI|GII)\.P(?:NA)?\d+)[-_ ]?(G?I(?:I)?\.(?:NA)?\d+)"):
    """Try to extract genotype from FASTA header using a regex pattern.

    Default pattern matches norovirus RdRp_VP1 format (e.g. GII.P17_GII.17).
    Override with --header-pattern for other viruses.
    """
    m = re.search(pattern, header)
    if m:
        return "_".join(m.groups()) if len(m.groups()) > 1 else m.group(1)
    return "unknown"


# ============================================================
# Main pipeline
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Generate representative reference sequences for viral mapping.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Example:
  python generate_representative_references.py \\
    --input-fasta ncbi_sequences.fasta \\
    --metadata metadata.tsv \\
    --metadata-genotype-column "rdrp_vp1_combo" \\
    --output reference_panel.fasta \\
    --threads 8
""")
    parser.add_argument("--input-fasta", required=True, help="Input sequences FASTA")
    parser.add_argument("--metadata", default="", help="Metadata TSV with genotype info")
    parser.add_argument("--metadata-genotype-column", default="genotype",
                        help="Column name in metadata TSV containing genotype labels")
    parser.add_argument("--header-pattern", default="",
                        help="Regex to infer genotype from FASTA header (fallback when no metadata)")
    parser.add_argument("--output", required=True, help="Output FASTA")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-dist-pct", type=float, default=MIN_DIST_WINDOW_PCT,
                        help="Minimum distinguishable window %% to keep refs separate (default: 30)")
    parser.add_argument("--cluster-threshold", type=float, default=CLUSTER_IDENTITY_THRESHOLD,
                        help="Identity threshold for initial clustering (default: 0.993)")
    parser.add_argument("--skip-qc", action="store_true",
                        help="Skip QC filtering (coverage, missing data, mixed sites)")
    parser.add_argument("--exclude-mediocre", action="store_true",
                        help="Also exclude mediocre-quality sequences (default: only exclude bad)")
    parser.add_argument("--output-prefix-label", default="MAPREF",
                        help="Label prefix for output FASTA headers (default: MAPREF)")
    parser.add_argument("--min-medoid-length", type=int, default=MIN_MEDOID_LENGTH,
                        help=f"Minimum ungapped length for medoid selection (default: {MIN_MEDOID_LENGTH}bp). "
                             "If the best medoid is shorter, the closest longer sequence is chosen instead.")
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path = output_path.with_suffix(".report.tsv")

    print(f"Loading sequences from {args.input_fasta}...")
    all_seqs = {}
    for rec in SeqIO.parse(args.input_fasta, "fasta"):
        all_seqs[rec.id] = str(rec.seq).upper()
    print(f"  Loaded {len(all_seqs)} sequences")

    # Determine genotypes
    geno_map = {}
    if args.metadata and os.path.exists(args.metadata):
        geno_map = parse_metadata(args.metadata, args.metadata_genotype_column)
        print(f"  Loaded {len(geno_map)} genotype assignments from metadata")

    # Fallback: infer from header
    if args.header_pattern:
        for sid in all_seqs:
            if sid not in geno_map:
                geno_map[sid] = infer_genotype_from_header(sid, args.header_pattern)

    # Group by genotype
    group_seqs = defaultdict(list)
    for sid, genotype in geno_map.items():
        if sid in all_seqs:
            group_seqs[genotype].append(sid)

    # Also include sequences without genotype assignment
    for sid in all_seqs:
        if sid not in geno_map:
            group_seqs["unknown"].append(sid)

    print(f"\nFound {len(group_seqs)} genotype groups:")
    for group in sorted(group_seqs, key=lambda c: -len(group_seqs[c])):
        if len(group_seqs[group]) >= 1:
            print(f"  {group}: {len(group_seqs[group])} sequences")

    # Process each group
    all_refs = []  # (group, medoid_id, medoid_seq, n_cluster)
    all_reports = []
    all_qc_rows = []

    for group in sorted(group_seqs):
        seq_ids = group_seqs[group]
        if not seq_ids:
            continue

        print(f"\n--- Processing {group} ({len(seq_ids)} seqs) ---")

        seqs_subset = {sid: all_seqs[sid] for sid in seq_ids if sid in all_seqs}
        if len(seqs_subset) == 0:
            continue

        # ---- QC filtering (before alignment) ----
        if not args.skip_qc:
            seqs_subset, qc_rows = qc_filter_sequences(
                seqs_subset,
                exclude_mediocre=args.exclude_mediocre,
            )
            for row in qc_rows:
                row["group"] = group
                all_qc_rows.append(row)
            n_excluded = len(seq_ids) - len(seqs_subset)
            if n_excluded > 0:
                print(f"  QC: excluded {n_excluded}/{len(seq_ids)} sequences "
                      f"({len(seqs_subset)} passed)")
            if len(seqs_subset) == 0:
                print(f"  WARNING: all sequences excluded by QC for {group}, skipping")
                continue
        # ---- QC end ----

        if len(seqs_subset) == 1:
            # Single sequence, use directly
            sid = list(seqs_subset.keys())[0]
            all_refs.append((group, sid, seqs_subset[sid], 1))
            all_reports.append({
                "group": group, "cluster_id": 1, "n_seqs": 1,
                "medoid": sid, "method": "single",
                "dist_check": "N/A (single ref)"
            })
            continue

        # MAFFT align within group
        tmp_in = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
        for sid, seq in seqs_subset.items():
            tmp_in.write(f">{sid}\n{seq}\n")
        tmp_in.close()
        tmp_aln = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)
        run_mafft(tmp_in.name, tmp_aln.name, fast=len(seqs_subset) > 50)
        os.unlink(tmp_in.name)

        aln_seqs = {}
        for rec in SeqIO.parse(tmp_aln.name, "fasta"):
            aln_seqs[rec.id] = str(rec.seq).upper()
        os.unlink(tmp_aln.name)

        # Cluster
        clusters = cluster_sequences(aln_seqs, args.cluster_threshold)
        print(f"  Clustered into {len(clusters)} groups at {args.cluster_threshold*100:.1f}% identity")

        # Build ungapped length map for medoid selection (prefer ≥6500bp)
        ungapped_lengths = {sid: len(seq) for sid, seq in seqs_subset.items()}

        # Select medoid per cluster
        medoids = []
        for ci, cluster in enumerate(clusters, 1):
            medoid_id = select_medoid(aln_seqs, cluster, ungapped_lengths)
            medoid_aln = aln_seqs[medoid_id]
            medoid_seq = ungapped_seq(medoid_aln)
            medoids.append((medoid_id, medoid_seq, len(cluster), aln_seqs[medoid_id]))

            all_reports.append({
                "group": group, "cluster_id": ci, "n_seqs": len(cluster),
                "medoid": medoid_id, "method": "medoid",
                "dist_check": "pending"
            })

        # Verify distinguishable windows between medoid pairs
        medoid_alns = {m[0]: m[3] for m in medoids}
        medoid_ids = [m[0] for m in medoids]
        kept_ids, merges = merge_close_references(list(medoid_ids), medoid_alns, args.min_dist_pct)

        if merges:
            print(f"  Merged {len(merges)} pairs with < {args.min_dist_pct*100:.0f}% distinguishable windows:")
            for merged_id, kept_id, pct in merges:
                print(f"    {merged_id[:30]} → merged into {kept_id[:30]} ({pct*100:.1f}% windows)")

        # Record final refs
        for mid in kept_ids:
            for m_id, m_seq, m_n, m_aln in medoids:
                if m_id == mid:
                    all_refs.append((group, mid, m_seq, m_n))
                    break

    # Write output FASTA
    with open(output_path, "w") as f:
        group_counter = defaultdict(int)
        for group, medoid_id, seq, n_cluster in all_refs:
            group_counter[group] += 1
            header = f"{group}_{args.output_prefix_label}_{medoid_id[:50]}"
            f.write(f">{header}\n{seq}\n")

    print(f"\n{'='*60}")
    print(f"Output: {output_path}")
    print(f"Total references: {len(all_refs)}")

    # Write clustering report
    with open(report_path, "w") as f:
        f.write("group\tcluster_id\tn_seqs\tmedoid\tmethod\tdist_check\n")
        for r in all_reports:
            f.write(f"{r['group']}\t{r['cluster_id']}\t{r['n_seqs']}\t{r['medoid']}\t{r['method']}\t{r['dist_check']}\n")
    print(f"Clustering report: {report_path}")

    # Write QC report
    if all_qc_rows:
        qc_path = output_path.with_suffix(".qc_report.tsv")
        qc_fields = ["seq_id", "group", "qc_status", "coverage", "missing", "mixed",
                     "cov_pct", "missing_pct", "mixed_sites"]
        with open(qc_path, "w") as f:
            f.write("\t".join(qc_fields) + "\n")
            for row in all_qc_rows:
                f.write("\t".join(str(row.get(k, "")) for k in qc_fields) + "\n")
        n_bad = sum(1 for r in all_qc_rows if r["qc_status"] == "bad")
        n_med = sum(1 for r in all_qc_rows if r["qc_status"] == "mediocre")
        n_good = sum(1 for r in all_qc_rows if r["qc_status"] == "good")
        print(f"QC report: {qc_path}")
        print(f"QC summary: {n_good} good, {n_med} mediocre, {n_bad} bad "
              f"({n_bad} excluded)")


if __name__ == "__main__":
    main()
