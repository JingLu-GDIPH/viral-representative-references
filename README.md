# Viral Representative References

Generate a panel of representative reference sequences from public database sequences (e.g. NCBI GenBank) for use in viral genotyping and consensus assembly pipelines.

## Overview

Given a FASTA of curated viral genome sequences and metadata assigning each sequence to a genotype/serotype, this tool:

1. **QC filters** sequences (coverage, missing data, mixed sites)
2. **MAFFT aligns** sequences within each genotype group
3. **Hierarchical clusters** at 99.3% identity
4. **Selects medoid** (real sequence closest to cluster center) per cluster
5. **Validates distinguishability** — verifies ≥30% of 150bp windows have ≥2 mismatches between medoid pairs; merges indistinguishable references
6. **Outputs** a compact reference panel FASTA + QC/clustering reports

The medoid approach ensures references stay close to real circulating sequences, unlike ASR/consensus sequences which can drift from observed diversity.

## Why Medoid References?

| Approach | Pros | Cons |
|----------|------|------|
| **Single reference** per genotype | Simple | Misses within-genotype diversity |
| **All NCBI sequences** | Maximum diversity | Too many (1000s), slow mapping, redundant |
| **ASR/consensus** | Synthetic "average" | May not match any real sequence |
| **Medoid (this tool)** | Real sequence, representative, compact | Requires clustering parameter tuning |

## Quick Start

```bash
python generate_representative_references.py \
  --input-fasta ncbi_sequences.fasta \
  --metadata metadata.tsv \
  --metadata-genotype-column rdrp_vp1_combo \
  --output reference_panel.fasta \
  --threads 8
```

## Input Format

### Sequences FASTA
Standard FASTA format. Headers can contain genotype info as fallback if metadata is incomplete.

### Metadata TSV
Tab-delimited with at least two columns:

| Column | Description |
|--------|-------------|
| `accession` | Sequence ID matching FASTA headers |
| `<genotype_column>` | Genotype/serotype label (specify via `--metadata-genotype-column`) |

Example:
```
accession	rdrp_vp1_combo	collection_date
LC036467	GII.P17_GII.17	2017-01-15
MK764013	GII.P16_GII.4	2016-01-11
```

## QC Checks

Three Nextclade-style checks are performed on each sequence **before** alignment:

| Check | Good | Mediocre | Bad (excluded) |
|-------|------|----------|----------------|
| **Coverage** | ≥90% ACGT | 80–90% | <80% |
| **Missing data** | <2% N | 2–10% | ≥10% |
| **Mixed sites** | ≤10 | — | >10 |

**Note:** Stop codon and frameshift checks are intentionally omitted because many viruses (e.g. norovirus, influenza) have overlapping ORFs in different reading frames, making single-frame translation unreliable. For ORF-level QC, use [seqqc](https://github.com/...) on individual ORF alignments.

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--cluster-threshold` | 0.993 | Identity threshold for clustering |
| `--min-dist-pct` | 0.30 | Min % distinguishable windows between refs |
| `--skip-qc` | off | Skip QC filtering |
| `--exclude-mediocre` | off | Also exclude mediocre sequences |
| `--output-prefix-label` | MAPREF | Label prefix for output headers |

## Output Files

| File | Content |
|------|---------|
| `reference_panel.fasta` | Selected medoid sequences |
| `reference_panel.report.tsv` | Clustering: group, cluster, medoid, n_seqs |
| `reference_panel.qc_report.tsv` | Per-sequence QC results |

## Clustering Science

The 99.3% clustering threshold and 30% distinguishable window requirement are designed for **150bp short-read mapping**:

- 150bp window with ≥2 mismatches ≈ 98.7% pairwise identity
- This is the minimum difference a Q20 150bp read can reliably distinguish
- ≥30% of windows meeting this criterion ensures EM-based read assignment has sufficient signal

## Dependencies

- Python ≥3.8
- Biopython (`pip install biopython`)
- NumPy (`pip install numpy`)
- MAFFT (`conda install -c bioconda mafft`)

## License

MIT
