# WhatsGNU-ATB
A custom reimplementation of [WhatsGNU](https://github.com/ahmedmagds/WhatsGNU) optimised for the scale of AllTheBacteria. It uses LMDB-backed sharded storage (8 shards) with numpy for hashing. The query tool is also custom-built for this database format.
Protein allele frequency analysis at the scale of [AllTheBacteria](https://allthebacteria.org/) (2.4M+ bacterial genomes).

WhatsGNU-ATB builds a sharded LMDB database from Bakta protein annotations and lets you query any bacterial genome to find out, for each protein, how many of the 2,438,285 genomes carry an identical copy — along with which species they belong to and which genomes are most similar.

A pre-built database covering all AllTheBacteria genomes is available on [OSF](https://osf.io/6jr4u/). If you just want to query genomes, skip to [Quick Start (Query)](#quick-start-query).

## Features

- **GNU scores**: for every protein in a query genome, reports the exact number of genomes (out of 2.4M+) containing an identical allele
- **Species breakdown**: top-K species contributing to each allele, with counts (other metadata like MLST contributions are coming soon)
- **Genome similarity**: ranks all 2.4M+ genomes by shared protein alleles with your query, identifying the closest relatives
- **Batch querying**: pass a directory of `.faa` files to query hundreds of genomes in one run
- **Sequence export**: optionally include the amino acid sequence in the output
- **Sharded LMDB backend**: 8 parallel shards with batched reads for fast lookups
- **Optional sequence storage**: store a representative amino acid sequence per allele hash in the database
- **Allele counts export**: dump the full allele frequency table as a TSV

## Installation

```bash
conda create -n whatsgnu-atb python=3.12
conda activate whatsgnu-atb
pip install numpy lmdb pandas
```

For publication figure generation, also install:

```bash
pip install matplotlib seaborn networkx adjustText scipy
```

## Quick Start (Query)

If you just want to query genomes against the pre-built AllTheBacteria database:

### 1. Download the database from OSF

```bash
# Download the database and metadata from https://osf.io/6jr4u/
# You need:
#   - WGNU_ATB_DB/ directory (lmdb_counts/ and lmdb_postings/ with 8 shards each)
#   - samples_with_ids.tsv (for species names and genome names)
```

### 2. Query a single genome

Your input must be a protein FASTA (`.faa`) file. See the [AllTheBacteria Bakta documentation](https://allthebacteria.readthedocs.io/en/latest/annotation.html) or the [Bakta GitHub](https://github.com/oschwengers/bakta) if you need to annotate your genome first.

**Basic query** (GNU scores only — fast, no postings needed):

```bash
python Query_WhatsGNU_ATB.py \
    --db_dir WGNU_ATB_DB/ \
    --shards 8 \
    --faa your_genome.bakta.faa \
    --out_dir results/
```

**Full query** (GNU scores + species breakdown + genome similarity):

```bash
python Query_WhatsGNU_ATB.py \
    --db_dir WGNU_ATB_DB/ \
    --shards 8 \
    --faa your_genome.bakta.faa \
    --include_sequence \
    --with_postings \
    --samples_tsv samples_with_ids.tsv \
    --species_names_tsv samples_with_ids.tsv \
    --top_k_species 5 \
    --top_k_genomes 10 \
    --out_dir results/
```

### 3. Query a batch of genomes

Pass a directory instead of a single file:

```bash
python Query_WhatsGNU_ATB.py \
    --db_dir WGNU_ATB_DB/ \
    --shards 8 \
    --faa directory_of_faa_files/ \
    --include_sequence \
    --with_postings \
    --out_dir results_batch/
```

## Query Output Files

### `<sample>.whatsgnu.tsv`

Per-protein results with one row per protein:

| Column | Description |
|---|---|
| `protein_id` | Protein identifier from the FASTA header |
| `allele_hash` | 128-bit BLAKE2b hash of the amino acid sequence (hex) |
| `sequence` | Amino acid sequence from the query genome (if `--include_sequence`) |
| `GNU_count` | Number of genomes containing this exact allele |
| `top5_species_names` | Top 5 species carrying this allele (if `--with_postings`) |
| `top5_species_counts` | Counts per species (if `--with_postings`) |
| `total_db_hits` | Total genomes in posting list |
| `hits_checked` | Number of postings actually decoded |

### `<sample>.similarity.tsv`

Genome similarity ranking (if `--with_postings`):

| Column | Description |
|---|---|
| `rank` | Rank by shared alleles (1 = most similar) |
| `genome_id` | Integer genome ID |
| `sample_name` | Sample accession (if `--samples_tsv` provided) |
| `species_id` | Species integer ID |
| `species_name` | Species name (if `--species_names_tsv` provided) |
| `shared_alleles` | Number of identical proteins shared with query |
| `percent_of_query` | Shared alleles as percentage of query proteome |

## Query Options Reference

| Option | Description | Default |
|---|---|---|
| `--db_dir` | Root database directory (required) | — |
| `--shards` | Number of shards, must be power of 2 (required) | — |
| `--faa` | Input `.faa` file or directory of `.faa` files (required) | — |
| `--out_dir` | Output directory (required) | — |
| `--with_postings` | Enable species breakdown and genome similarity | off |
| `--include_sequence` | Include amino acid sequence in output | off |
| `--top_k_species` | Number of top species to report per protein | 5 |
| `--top_k_genomes` | Number of top similar genomes to report | 10 |
| `--postings_limit` | Max genome IDs to decode per allele (0 = all) | 0 |
| `--species_names_tsv` | TSV mapping SpeciesID → species name | none |
| `--samples_tsv` | TSV mapping SampleID → sample accession | none |

## Interpreting GNU Scores

| GNU Score Range | Interpretation |
|---|---|
| >100,000 | Highly conserved ubiquitous allele |
| 1000–10,000 | Common allele |
| 1–100 | Rare allele, likely strain-specific |
| 0 | Unique to the query genome — not in any AllTheBacteria genome |

## Building a Database

To build a new database from scratch (e.g., for a custom genome set):

### Input Requirements

A sample table TSV with these columns:

| Column | Description |
|---|---|
| `SampleID` | Unique integer ID per genome |
| `Sample` | Sample name (used to find `.faa` file) |
| `SpeciesID` | Integer species ID |

Optional column: `faa_path` (full path to FAA file). If absent, uses `--faa_dir/<Sample><faa_suffix>`.

### Build Command

```bash
python WhatsGNU_ATB_DB.py \
    --sample_table samples_with_ids.tsv \
    --faa_dir /path/to/faa_files/ \
    --out_dir WGNU_ATB_DB/ \
    --tmp_dir /scratch/tmp/ \
    --shards 8 \
    --with_postings \
    --sort_mem_mb 65536 \
    --lmdb_map_gb_counts_per_shard 24 \
    --lmdb_map_gb_postings_per_shard 160 \
    --export_allele_counts allele_counts.tsv \
    --log_file build.log \
    --log_level INFO
```

### Build with Sequences

To also store representative amino acid sequences per allele hash:

```bash
python WhatsGNU_ATB_DB.py \
    --sample_table samples_with_ids.tsv \
    --faa_dir /path/to/faa_files/ \
    --out_dir WGNU_ATB_DB/ \
    --shards 8 \
    --with_postings \
    --with_sequences \
    --lmdb_map_gb_sequences_per_shard 25 \
    --log_level INFO
```

## Build Options Reference

| Option | Description | Default |
|---|---|---|
| `--sample_table` | Sample table TSV (required) | — |
| `--faa_dir` | Directory of `.faa` files | none |
| `--out_dir` | Output directory (required) | — |
| `--tmp_dir` | Temp directory for intermediate files | `<out_dir>/tmp` |
| `--reduce_tmp_dir` | Local scratch for sort/reduce (faster I/O) | none |
| `--shards` | Number of shards, power of 2 | 16 |
| `--with_postings` | Build posting lists (genome IDs per allele) | off |
| `--with_sequences` | Store representative AA sequence per allele | off |
| `--faa_suffix` | Suffix appended to Sample name for FAA lookup | `.bakta.faa` |
| `--sort_mem_mb` | RAM for external sort per shard (MB) | 65536 |
| `--lmdb_map_gb_counts_per_shard` | LMDB map size for counts (GB) | 24 |
| `--lmdb_map_gb_postings_per_shard` | LMDB map size for postings (GB) | 160 |
| `--lmdb_map_gb_sequences_per_shard` | LMDB map size for sequences (GB) | 25 |
| `--export_allele_counts` | Path to write allele frequency TSV | none |
| `--parse_only` | Only parse FAA → record bins, skip reduce | off |
| `--reduce_only` | Only reduce existing record bins → LMDB | off |
| `--resume` | Auto-detect: skip parse if record bins exist | off |
| `--skip_existing_shards` | Skip shards with existing LMDB output | off |
| `--log_file` | Log file path | `<out_dir>/build.log` |
| `--log_level` | Logging level | INFO |

## Database Structure

```
WGNU_ATB_DB/
├── lmdb_counts/
│   ├── shard_00/         # LMDB: hash → (func_id, GNU_count)
│   ├── shard_01/
│   └── ...
├── lmdb_postings/        # (if --with_postings)
│   ├── shard_00/         # LMDB: hash → varint-encoded genome IDs
│   ├── shard_01/
│   └── ...
├── lmdb_sequences/       # (if --with_sequences)
│   ├── shard_00/         # LMDB: hash → amino acid sequence (UTF-8)
│   └── ...
├── indexes/
│   └── genome_species.u32   # Binary array: genome_id → species_id
└── metadata/
    ├── build_info.json       # Build parameters, stats, version
    └── functions.tsv.gz      # Function ID → function description
```

## Technical Details

- **Hashing**: BLAKE2b with 128-bit (16-byte) digest of the amino acid sequence
- **Sharding**: `shard_id = first_byte(hash) & (num_shards - 1)`
- **GNU count**: number of genomes containing an allele at least once (deduplicated within each genome)
- **Postings**: delta + varint encoded sorted unique genome IDs
- **External sort**: numpy structured arrays for memory-efficient sorting; batched multi-pass merge with fanin of 64
- **Query optimizations**: batched LMDB reads (one transaction per shard), numpy-vectorized species lookups, partial argsort for top-K genome ranking

## Resource Requirements

### Building (2.4M genomes, 8 shards)

| Resource | Recommendation |
|---|---|
| RAM | 250–500 GB |
| CPUs | 4–6 cores |
| Disk (tmp) | ~2 TB scratch |
| Wall time | 6–24 hours (I/O dependent) |

### Querying

| Resource | Recommendation |
|---|---|
| RAM | ~2 GB (basic) / ~4 GB (with postings) |
| Wall time | ~5–150 seconds per genome |

## Citation

If you use WhatsGNU-ATB in your research, please cite:

> Moustafa AM and Planet PJ. WhatsGNU: a tool for identifying proteomic novelty. *Genome Biology*, 2020. [doi:10.1186/s13059-020-01965-w](https://doi.org/10.1186/s13059-020-01965-w)
> Hunt M, Lima L, Shen W, Lees J, Iqbal Z. AllTheBacteria - all bacterial genomes assembled, available and searchable. *bioRxiv*, 2024.[https://doi.org/10.1101/2024.03.08.584059](https://doi.org/10.1101/2024.03.08.584059)

> AllTheBacteria — all bacterial genomes assembled, available and searchable. *bioRxiv*, 2024. [doi:10.1101/2024.03.08.584059](https://doi.org/10.1101/2024.03.08.584059)

## License

MIT
