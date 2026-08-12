# Viral Representative Reference Panel Generation — Method Documentation

## 1. Overview

This document describes the complete methodology for generating a compact
panel of representative reference sequences from public database sequences
(e.g., NCBI GenBank). The approach was developed for norovirus GI/GII
genotyping from sewage metagenomic data but is applicable to any virus with
genotype/serotype classification.

**Core idea**: From thousands of public sequences per genotype, select a small
set (~100-300) of **real NCBI sequences** (medoids) that:
1. Represent within-genotype diversity (cluster coverage)
2. Are distinguishable from each other by short reads (≥30% of 150bp windows
   have ≥2 mismatches between any reference pair)
3. Pass sequence quality control (coverage, missing data, mixed sites)

This ensures that EM-based read assignment during downstream consensus
assembly has sufficient signal to correctly assign reads to references.

---

## 2. Input Data

### 2.1 Sequence source

- **Database**: NCBI GenBank nucleotide database
- **Filter**: Sequences ≥5800 bp (covering ORF1 + ORF2/VP1)
- **Date range**: Collected from year 2000 onward
- **Example**: 3,775 GII sequences, 575 GI sequences

### 2.2 Metadata

A tab-delimited TSV file with genotype assignments:

| Column | Description |
|--------|-------------|
| `accession` | GenBank accession matching FASTA headers |
| `vp1_genotype` | VP1 capsid genotype (e.g., GII.4, GII.17) |
| `vp1_assignment_status` | Confidence of genotype assignment |
| `vp1_percent_identity` | % identity to best reference genotype |

VP1 genotyping was performed using BLAST against a curated reference set
for each known genotype. Sequences with `no_significant_vp1_hit` status
(7 out of 3,775 GII) are excluded.

---

## 3. Pipeline Steps

### Step 1: Sequence Quality Control (QC)

**Purpose**: Remove low-quality sequences before they contaminate the
alignment and clustering.

**Timing**: Performed on ungapped sequences **before** MAFFT alignment,
within each genotype group.

**Three checks** (adapted from Nextclade/ncov QC framework):

| Check | Metric | Good | Mediocre | Bad (excluded) |
|-------|--------|------|----------|-----------------|
| **Coverage** | ACGT bases / total length | ≥90% | 80–90% | <80% |
| **Missing data** | N count / total length | <2% | 2–10% | ≥10% |
| **Mixed sites** | Non-ACGTN ambiguous bases (R, Y, K, M, S, W, B, D, H, V) | ≤10 | — | >10 |

**Overall QC status** = worst of all 3 checks.
- `bad` sequences are **excluded** from downstream processing
- `mediocre` sequences are **retained** (can be excluded with `--exclude-mediocre`)
- `good` sequences are fully retained

**Rationale for 3 checks (not 5)**:

The full Nextclade QC framework includes 5 checks (adding stop codons and
frameshifts). We intentionally omit these two because:

1. **Overlapping ORFs**: Norovirus (and many other viruses) has 3 overlapping
   ORFs (ORF1/ORF2/ORF3) in different reading frames. Translating the full
   genome in any single frame produces spurious stop codons at ORF boundaries.
2. **Variable 5' start positions**: NCBI sequences start at varying positions
   (some include 5'UTR, others start mid-ORF1), making frame determination
   unreliable without per-sequence annotation.
3. **Low marginal value**: For reference panel generation, the 3 retained
   checks (coverage, missing data, mixed sites) already capture the main
   quality issues. Frame-level QC is more appropriate for phylogenetic
   analysis of individual ORFs.

**Validation on norovirus GII** (3,775 sequences):

| QC Status | Count | Percentage |
|-----------|:-----:|:----------:|
| Good | 3,646 | 96.6% |
| Mediocre | 69 | 1.8% |
| **Bad (excluded)** | **60** | **1.6%** |

Main exclusion reasons: low coverage (<80% ACGT), excessive missing data
(≥10% N), high mixed site count (>10 ambiguous bases).

### Step 2: Multiple Sequence Alignment (MAFFT)

**Purpose**: Align all QC-passed sequences within each genotype group to
enable pairwise identity calculation and clustering.

**Method**: MAFFT with `--auto` (default FFT-NS-2 iterative refinement).
For groups with >50 sequences, `--retree 1` is used for faster tree
construction (L-INS-i is too slow for 1000+ sequences).

**Scope**: Alignment is performed **within each genotype group independently**
(e.g., all GII.4 sequences aligned together, all GII.17 together). This
ensures:
- Alignment quality is high (sequences are similar within a genotype)
- Computational cost is manageable (not aligning 3,775 sequences at once)
- Clustering captures within-genotype micro-diversity

### Step 3: Hierarchical Clustering

**Purpose**: Group highly similar sequences (>99.3% identity) into clusters,
so each cluster can be represented by a single medoid.

**Method**: Greedy single-linkage clustering:

1. Compute pairwise identity for all sequences in the alignment
   (excluding gaps and N positions)
2. Sort sequences by decreasing average similarity (denser sequences first)
3. For each unassigned sequence (seed), assign it and all unassigned
   sequences with ≥99.3% identity to a new cluster
4. Repeat until all sequences are assigned

**Threshold rationale (99.3%)**:

The 99.3% identity threshold corresponds to ~52 SNPs across a 7,500bp
norovirus genome. This captures epidemiologically meaningful clusters
(transmission chains / variant lineages) without splitting at the level
of individual sequencing errors (~1-2 SNPs per genome).

For 150bp short reads, 99.3% genome identity means most 150bp windows
have 0-1 mismatches — reads from these sequences would be indistinguishable
by mapping. Therefore, one representative per cluster is sufficient.

### Step 4: Medoid Selection

**Purpose**: For each cluster, select the single real sequence that best
represents all sequences in the cluster.

**Method**: The medoid is the sequence with the **minimum average distance**
to all other sequences in the cluster:

```
medoid = argmin_i  mean(1 - identity(i, j))  for all j in cluster
```

**Why medoid (not consensus)?**:

| Approach | Example | Issue |
|----------|---------|-------|
| Consensus | Majority-vote at each position | May not match any real sequence |
| ASR | Ancestral state reconstruction | Even more distant from real samples |
| **Medoid** | **Real NCBI sequence** | **Closest to circulating diversity** |

A medoid is guaranteed to be a real observed sequence. When sample reads
are mapped against medoid references, the reads come from circulating
viruses that are epidemiologically close to at least one medoid.

### Step 5: Distinguishability Validation

**Purpose**: Ensure that selected medoid pairs are sufficiently different
for short-read mapping to distinguish them.

**Method**: Sliding window analysis on the alignment:

1. For each pair of medoids, slide a 150bp window across the alignment
   (step size = 50bp)
2. Count windows with **≥2 mismatches** (both positions must be ACGT)
3. If **<30%** of windows are distinguishable, merge the pair (remove one
   medoid, keeping the one from the larger cluster)

**Parameters**:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Window size | 150 bp | Matches typical Illumina short-read length |
| Step size | 50 bp | Balance between resolution and compute |
| Min mismatches | 2 | Q20 150bp read can reliably detect ≥2 differences |
| Min dist windows | 30% | Ensures EM has sufficient signal across the genome |

**Scientific basis**: At Q20 (1% error rate per base), a 150bp read has
an expected 1.5 sequencing errors. Two mismatches between references is
above this noise floor, meaning a read can be correctly assigned to one
reference over the other with >95% confidence in that window. Requiring
30% of windows to meet this criterion ensures the entire genome has
enough distinguishing positions for EM convergence.

### Step 6: Output

**FASTA format** with informative headers:
```
>{genotype}_{LABEL}_{accession}
{sequence}
```

Example: `>GII.4_MAPREF_MK764013_11-Jan-2016`

**Clustering report** (`*.report.tsv`):
- Group, cluster ID, cluster size, medoid accession, method

**QC report** (`*.qc_report.tsv`):
- Per-sequence QC status (coverage, missing, mixed) for all input sequences

---

## 4. Parameters Summary

| Parameter | Default | Command-line flag | Description |
|-----------|---------|-------------------|-------------|
| Cluster identity threshold | 0.993 | `--cluster-threshold` | Sequences ≥99.3% identical are clustered |
| Min distinguishable windows | 30% | `--min-dist-pct` | Min % of 150bp windows with ≥2 mismatches |
| QC coverage good/bad | 90%/80% | `--cov-good` / `--cov-bad` | ACGT coverage thresholds |
| QC missing good/bad | 2%/10% | `--miss-good` / `--miss-bad` | N fraction thresholds |
| QC mixed threshold | 10 | `--mixed-threshold` | Max allowed ambiguous bases |
| Exclude mediocre | off | `--exclude-mediocre` | Also exclude mediocre sequences |

---

## 5. Validation Results

### 5.1 Norovirus GII

| Metric | Value |
|--------|-------|
| Input sequences | 3,775 |
| QC passed | 3,715 (98.4%) |
| QC excluded (bad) | 60 (1.6%) |
| Genotype groups | 16+ |
| Clusters (at 99.3%) | ~200 |
| **Final reference panel** | **~288 sequences** |
| **Reduction ratio** | **92.4%** (3,775 → ~288) |

### 5.2 Cluster coverage

The medoid panel achieves 100% cluster coverage: every sequence in the
input is within 99.3% identity of at least one panel member.

### 5.3 Distinguishability

All medoid pairs have ≥30% distinguishable 150bp windows (≥2 mismatches),
ensuring EM-based read assignment can correctly distinguish references.

---

## 6. Usage

```bash
# Standard run with QC filtering
python generate_representative_references.py \
  --input-fasta data/ncbi_norovirus_gii_gt5800_collected_2000_onward.fasta \
  --metadata data/ncbi_norovirus_gii_gt5800_metadata_vp1_typed.tsv \
  --metadata-genotype-column vp1_genotype \
  --output output/gii_mapping_references_v2.fasta \
  --threads 8

# Strict mode: exclude both bad and mediocre
python generate_representative_references.py \
  ... \
  --exclude-mediocre

# Skip QC (use all sequences)
python generate_representative_references.py \
  ... \
  --skip-qc
```

---

## 7. Limitations

1. **Genotype labeling dependency**: The tool requires pre-assigned genotype
   labels. Misclassified sequences will contaminate genotype groups.
2. **No recombination detection**: Intra-genotype recombinants are not
   specifically flagged (though they may form separate clusters).
3. **MAFFT scalability**: Groups with >1,000 sequences (e.g., GII.4 with
   1,739 sequences) take significant alignment time (~30-60 min).
4. **Reference currency**: The panel reflects NCBI deposits at the time of
   download. Emerging variants may require periodic panel updates.
