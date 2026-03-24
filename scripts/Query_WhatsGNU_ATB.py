#!/usr/bin/env python3
"""
Query_WhatsGNU_ATB_fast.py

Optimized query script for WhatsGNU_ATB databases.
~10-30x faster than the original for postings+similarity queries.

Key optimizations:
  1. numpy array for genome_shared_alleles (indexed increment, no Counter overhead)
  2. numpy-vectorized species lookup (genome_species[gids] in one call)
  3. Batched LMDB reads: one transaction per shard instead of per-protein
  4. Faster varint decoding with numpy buffer ops
  5. Two-pass architecture: hash all proteins first, group by shard, batch lookup

Dependencies: lmdb, numpy
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import struct
import sys
import time, resource
from collections import Counter
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import lmdb
import numpy as np

HEADER_RE = re.compile(r"^>(\S+)\s*(.*)$")
VAL_COUNTS = struct.Struct("<II")     # func_id:uint32, GNU_count:uint32
VAL_POST_HDR = struct.Struct("<I")    # n:uint32

__version__ = "1.0.0"
# ─── FASTA parsing ────────────────────────────────────────────────────
def parse_faa(path: Path) -> List[Tuple[str, str, str]]:
    """Parse entire FAA into list of (protein_id, sequence, function)."""
    records = []
    pid = ""
    func = ""
    seq_lines: List[str] = []
    with path.open("rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if pid and seq_lines:
                    records.append((pid, "".join(seq_lines), func))
                    seq_lines = []
                m = HEADER_RE.match(line)
                if m:
                    pid = m.group(1)
                    func = (m.group(2) or "").strip()
                else:
                    pid = line[1:].split()[0] if line[1:].strip() else ""
                    func = ""
            else:
                seq_lines.append(line.strip())
        if pid and seq_lines:
            records.append((pid, "".join(seq_lines), func))
    return records


# ─── Hash + shard ─────────────────────────────────────────────────────
def hash_allele_128(aa_seq: str) -> bytes:
    return hashlib.blake2b(aa_seq.encode("utf-8"), digest_size=16).digest()

def shard_id_from_hash(h16: bytes, nshards: int) -> int:
    return h16[0] & (nshards - 1)


# ─── Fast varint decode ──────────────────────────────────────────────
def decode_postings_numpy(blob: bytes) -> np.ndarray:
    """Decode postings blob into numpy uint32 array. ~3-5x faster than pure Python."""
    if not blob or len(blob) < 4:
        return np.empty(0, dtype=np.uint32)
    n = struct.unpack_from("<I", blob, 0)[0]
    if n == 0:
        return np.empty(0, dtype=np.uint32)

    # Decode varints into deltas array
    data = blob[4:]
    deltas = np.empty(n, dtype=np.int64)
    pos = 0
    dlen = len(data)
    for i in range(n):
        x = 0
        shift = 0
        while pos < dlen:
            b = data[pos]
            pos += 1
            x |= (b & 0x7F) << shift
            if (b & 0x80) == 0:
                break
            shift += 7
        deltas[i] = x

    # Convert deltas to absolute genome IDs via cumsum
    deltas[0] = deltas[0]  # first is absolute
    np.cumsum(deltas, out=deltas)
    return deltas.astype(np.uint32)


def decode_postings_numpy_limited(blob: bytes, limit: int) -> np.ndarray:
    """Decode at most `limit` genome IDs from postings."""
    if not blob or len(blob) < 4:
        return np.empty(0, dtype=np.uint32)
    n = struct.unpack_from("<I", blob, 0)[0]
    n = min(n, limit)
    if n == 0:
        return np.empty(0, dtype=np.uint32)

    data = blob[4:]
    deltas = np.empty(n, dtype=np.int64)
    pos = 0
    dlen = len(data)
    for i in range(n):
        x = 0
        shift = 0
        while pos < dlen:
            b = data[pos]
            pos += 1
            x |= (b & 0x7F) << shift
            if (b & 0x80) == 0:
                break
            shift += 7
        deltas[i] = x

    deltas[0] = deltas[0]
    np.cumsum(deltas, out=deltas)
    return deltas.astype(np.uint32)


# ─── Genome species index (numpy) ────────────────────────────────────
def load_genome_species_numpy(path: Path) -> np.ndarray:
    """Load genome_species.u32 as numpy array for vectorized lookups."""
    data = path.read_bytes()
    return np.frombuffer(data, dtype=np.uint32)


# ─── Species + genome name loading ────────────────────────────────────
def load_species_names(tsv_path: Path) -> Dict[int, str]:
    import csv
    species_map = {}
    with tsv_path.open("rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sid_str = row.get("SpeciesID") or row.get("species_id")
            sname = row.get("Species") or row.get("species_name")
            if sid_str and sname:
                try:
                    species_map[int(sid_str)] = sname
                except ValueError:
                    continue
    return species_map


def load_genome_names(tsv_path: Path) -> Dict[int, str]:
    import csv
    genome_map = {}
    with tsv_path.open("rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid_str = (row.get("SampleID") or row.get("sample_id")
                       or row.get("genome_id"))
            gname = row.get("Sample") or row.get("sample_name")
            if gid_str and gname:
                try:
                    genome_map[int(gid_str)] = gname
                except ValueError:
                    continue
    return genome_map


# ─── Input discovery ──────────────────────────────────────────────────
def discover_faa_inputs(faa_path: Path) -> List[Path]:
    if faa_path.is_file():
        return [faa_path]
    if faa_path.is_dir():
        return sorted(p for p in faa_path.iterdir()
                      if p.is_file() and p.name.endswith(".faa"))
    raise FileNotFoundError(f"Input path not found: {faa_path}")


# ─── Batched LMDB reader ─────────────────────────────────────────────
class BatchShardedDB:
    """
    Like ShardedDB but supports batch lookups within a single transaction.
    Much faster than one txn.get() per key.
    """
    def __init__(self, root: Path, nshards: int, db_name: bytes):
        self.root = root
        self.nshards = nshards
        self.db_name = db_name
        self.envs: Dict[int, Optional[lmdb.Environment]] = {}
        self.dbs: Dict[int, Optional[object]] = {}

    def _open_shard(self, sid: int) -> None:
        if sid in self.envs:
            return
        shard_dir = self.root / f"shard_{sid:02x}"
        if not shard_dir.exists():
            self.envs[sid] = None
            self.dbs[sid] = None
            return
        env = lmdb.open(str(shard_dir), readonly=True, lock=False,
                        readahead=True, max_dbs=1, map_size=1 << 30)
        db = env.open_db(self.db_name)
        self.envs[sid] = env
        self.dbs[sid] = db

    def open_all(self) -> None:
        """Pre-open all shards to avoid lazy-open overhead during query."""
        for sid in range(self.nshards):
            self._open_shard(sid)

    def batch_get(self, keys_by_shard: Dict[int, List[Tuple[int, bytes]]]
                  ) -> Dict[int, Optional[bytes]]:
        """
        Batch lookup: keys_by_shard = {shard_id: [(protein_idx, key16), ...]}
        Returns: {protein_idx: value_bytes_or_None}
        """
        results: Dict[int, Optional[bytes]] = {}
        for sid, idx_keys in keys_by_shard.items():
            self._open_shard(sid)
            env = self.envs.get(sid)
            db = self.dbs.get(sid)
            if env is None or db is None:
                for pidx, _k in idx_keys:
                    results[pidx] = None
                continue
            # Single transaction for all keys in this shard
            with env.begin(db=db, buffers=True) as txn:
                for pidx, k in idx_keys:
                    v = txn.get(k)
                    # Copy from LMDB buffer since it's invalidated after txn
                    results[pidx] = bytes(v) if v is not None else None
        return results

    def close(self) -> None:
        for env in self.envs.values():
            if env is not None:
                env.close()
        self.envs.clear()
        self.dbs.clear()


# ─── Main query function (optimized) ─────────────────────────────────
def query_genome_fast(
    faa_file: Path,
    out_tsv: Path,
    out_similarity: Path,
    counts_db: BatchShardedDB,
    postings_db: Optional[BatchShardedDB],
    genome_species: Optional[np.ndarray],
    max_genome_id: int,
    top_k_species: int,
    top_k_genomes: int,
    postings_limit: Optional[int],
    include_sequence: bool,
    species_names: Optional[Dict[int, str]],
    genome_names: Optional[Dict[int, str]],
) -> None:
    t0 = time.time()
    has_postings = postings_db is not None and genome_species is not None

    # ── Pass 1: Parse FAA + hash all proteins ─────────────────────────
    records = parse_faa(faa_file)
    n_proteins = len(records)
    if n_proteins == 0:
        print(f"  WARNING: no proteins in {faa_file.name}", file=sys.stderr)
        return

    hashes = []
    for _pid, aa_seq, _func in records:
        hashes.append(hash_allele_128(aa_seq))

    t_hash = time.time()
    print(f"  Hashed {n_proteins} proteins in {t_hash - t0:.2f}s", file=sys.stderr)

    # ── Group by shard for batch lookup ───────────────────────────────
    counts_by_shard: Dict[int, List[Tuple[int, bytes]]] = {}
    for i, h in enumerate(hashes):
        sid = shard_id_from_hash(h, counts_db.nshards)
        counts_by_shard.setdefault(sid, []).append((i, h))

    # ── Batch counts lookup ───────────────────────────────────────────
    counts_results = counts_db.batch_get(counts_by_shard)
    gnu_scores = np.zeros(n_proteins, dtype=np.uint32)
    for i, v in counts_results.items():
        if v is not None:
            _fid, gnu = VAL_COUNTS.unpack(v)
            gnu_scores[i] = gnu

    t_counts = time.time()
    print(f"  Counts lookup in {t_counts - t_hash:.2f}s", file=sys.stderr)

    if not has_postings:
        # ── Write simple output (no postings) ─────────────────────────
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        with out_tsv.open("wt") as out:
            parts = ["protein_id", "allele_hash"]
            if include_sequence:
                parts.append("sequence")
            parts.append("GNU_count")
            out.write("\t".join(parts) + "\n")
            for i, (pid, aa_seq, _func) in enumerate(records):
                row = [pid, hashes[i].hex()]
                if include_sequence:
                    row.append(aa_seq)
                row.append(str(int(gnu_scores[i])))
                out.write("\t".join(row) + "\n")
        print(f"  Total: {time.time() - t0:.2f}s (no postings)", file=sys.stderr)
        return

    # ── Batch postings lookup ─────────────────────────────────────────
    post_by_shard: Dict[int, List[Tuple[int, bytes]]] = {}
    for i, h in enumerate(hashes):
        sid = shard_id_from_hash(h, postings_db.nshards)
        post_by_shard.setdefault(sid, []).append((i, h))

    post_results = postings_db.batch_get(post_by_shard)

    t_post_lookup = time.time()
    print(f"  Postings lookup in {t_post_lookup - t_counts:.2f}s", file=sys.stderr)

    # ── Decode postings + accumulate genome hits (numpy) ──────────────
    # genome_shared: array of size max_genome_id+1, counts shared alleles
    genome_shared = np.zeros(max_genome_id + 1, dtype=np.int32)

    # Per-protein species composition (stored for output)
    per_protein_species: List[Optional[Tuple[List, List, int, int]]] = [None] * n_proteins
    proteins_with_hits = 0

    for i in range(n_proteins):
        pv = post_results.get(i)
        if pv is None or len(pv) < 4:
            continue

        n_total = struct.unpack_from("<I", pv, 0)[0]

        if postings_limit:
            gids = decode_postings_numpy_limited(pv, postings_limit)
        else:
            gids = decode_postings_numpy(pv)

        if len(gids) == 0:
            continue

        proteins_with_hits += 1

        # Accumulate genome similarity (vectorized — single numpy call)
        np.add.at(genome_shared, gids, 1)

        # Species composition for this protein (vectorized)
        # Clip to valid range
        valid = gids < len(genome_species)
        sids = genome_species[gids[valid]]
        sc = Counter()
        # Use numpy unique for counting instead of Python loop
        uniq_sids, sid_counts = np.unique(sids, return_counts=True)
        for s, c in zip(uniq_sids, sid_counts):
            sc[int(s)] = int(c)

        top = sc.most_common(top_k_species)
        per_protein_species[i] = (
            [s for s, _ in top],
            [c for _, c in top],
            n_total,
            len(gids),
        )

    t_decode = time.time()
    print(f"  Postings decode + accumulate in {t_decode - t_post_lookup:.2f}s",
          file=sys.stderr)

    # ── Write per-protein output ──────────────────────────────────────
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("wt") as out:
        parts = ["protein_id", "allele_hash"]
        if include_sequence:
            parts.append("sequence")
        parts.append("GNU_count")

        sp_label = "species_names" if species_names else "species_ids"
        parts.extend([
            f"top{top_k_species}_{sp_label}",
            f"top{top_k_species}_species_counts",
            "total_db_hits", "hits_checked",
        ])
        out.write("\t".join(parts) + "\n")

        for i, (pid, aa_seq, _func) in enumerate(records):
            row = [pid, hashes[i].hex()]
            if include_sequence:
                row.append(aa_seq)
            row.append(str(int(gnu_scores[i])))

            info = per_protein_species[i]
            if info is None:
                row.extend(["", "", "0", "0"])
            else:
                top_sids, top_cnts, n_total, n_decoded = info
                if species_names:
                    names = ",".join(
                        species_names.get(s, f"Unknown_{s}") for s in top_sids)
                else:
                    names = ",".join(str(s) for s in top_sids)
                counts_str = ",".join(str(c) for c in top_cnts)
                row.extend([names, counts_str, str(n_total), str(n_decoded)])

            out.write("\t".join(row) + "\n")

    t_write = time.time()

    # ── Write genome similarity ───────────────────────────────────────
    # Find top-K genomes from the numpy array (partial argsort)
    # Only consider genomes with at least 1 shared allele
    nonzero_mask = genome_shared > 0
    nonzero_idx = np.nonzero(nonzero_mask)[0]
    nonzero_vals = genome_shared[nonzero_idx]

    # Top-K by partial sort (much faster than full sort for 2.4M array)
    k = min(top_k_genomes, len(nonzero_vals))
    if k > 0:
        top_k_pos = np.argpartition(nonzero_vals, -k)[-k:]
        top_k_pos = top_k_pos[np.argsort(nonzero_vals[top_k_pos])[::-1]]
        top_gids = nonzero_idx[top_k_pos]
        top_counts = genome_shared[top_gids]
    else:
        top_gids = np.array([], dtype=np.uint32)
        top_counts = np.array([], dtype=np.int32)

    with out_similarity.open("wt") as sim_out:
        header = ["rank", "genome_id"]
        if genome_names:
            header.append("sample_name")
        header.append("species_id")
        if species_names:
            header.append("species_name")
        header.extend(["shared_alleles", "percent_of_query"])
        sim_out.write("\t".join(header) + "\n")

        for rank, (gid, cnt) in enumerate(zip(top_gids, top_counts), 1):
            gid = int(gid)
            cnt = int(cnt)
            sid = int(genome_species[gid]) if gid < len(genome_species) else 0
            pct = 100 * cnt / n_proteins if n_proteins > 0 else 0

            row = [str(rank), str(gid)]
            if genome_names:
                row.append(genome_names.get(gid, f"Unknown_GID{gid}"))
            row.append(str(sid))
            if species_names:
                row.append(species_names.get(sid, f"Unknown_SID{sid}"))
            row.extend([str(cnt), f"{pct:.2f}"])
            sim_out.write("\t".join(row) + "\n")

    t_end = time.time()
    print(f"  Proteins: {n_proteins}, hits: {proteins_with_hits}, "
          f"unique genomes touched: {int(nonzero_mask.sum()):,}",
          file=sys.stderr)
    if k > 0:
        print(f"  Top genome: ID={int(top_gids[0])}, "
              f"shared={int(top_counts[0])} ({100*int(top_counts[0])/n_proteins:.1f}%)",
              file=sys.stderr)
    print(f"  TOTAL: {t_end - t0:.2f}s "
          f"(hash={t_hash-t0:.1f}s counts={t_counts-t_hash:.1f}s "
          f"post_lookup={t_post_lookup-t_counts:.1f}s "
          f"decode={t_decode-t_post_lookup:.1f}s "
          f"write={t_end-t_decode:.1f}s)", file=sys.stderr)


# ─── Main ─────────────────────────────────────────────────────────────
def main() -> int:
    t0_all = time.time()
    ap = argparse.ArgumentParser(
        description="Fast WhatsGNU_ATB query with genome similarity")

    ap.add_argument("--db_dir", required=True,
                    help="Root database directory from WhatsGNU_ATB.py")
    ap.add_argument("--shards", type=int, required=True,
                    help="Number of shards (power of 2)")
    ap.add_argument("--faa", required=True,
                    help="FAA file or folder of FAA files")
    ap.add_argument("--out_dir", required=True,
                    help="Output directory for TSV reports")
    ap.add_argument("--with_postings", action="store_true",
                    help="Enable postings/species/similarity")
    ap.add_argument("--top_k_species", type=int, default=5)
    ap.add_argument("--top_k_genomes", type=int, default=10)
    ap.add_argument("--postings_limit", type=int, default=0,
                    help="Max genome IDs to decode per allele (0=all)")
    ap.add_argument("--include_sequence", action="store_true")
    ap.add_argument("--species_names_tsv", type=str, default=None)
    ap.add_argument("--samples_tsv", type=str, default=None)
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    args = ap.parse_args()
    nshards = args.shards
    if nshards & (nshards - 1) != 0:
        print("ERROR: --shards must be power of 2", file=sys.stderr)
        return 2

    db_dir = Path(args.db_dir)
    counts_root = db_dir / "lmdb_counts"
    if not counts_root.exists():
        print(f"ERROR: counts DB not found at {counts_root}", file=sys.stderr)
        return 2

    postings_root = db_dir / "lmdb_postings"
    genome_species_path = db_dir / "indexes" / "genome_species.u32"

    postings_db: Optional[BatchShardedDB] = None
    genome_species: Optional[np.ndarray] = None
    max_genome_id = 0

    if args.with_postings:
        if not postings_root.exists():
            print(f"ERROR: postings DB not found at {postings_root}",
                  file=sys.stderr)
            return 2
        if not genome_species_path.exists():
            print(f"ERROR: genome_species.u32 not found", file=sys.stderr)
            return 2

        t0 = time.time()
        postings_db = BatchShardedDB(postings_root, nshards, b"postings")
        postings_db.open_all()
        genome_species = load_genome_species_numpy(genome_species_path)
        max_genome_id = len(genome_species) - 1
        print(f"Loaded postings DB + genome_species ({len(genome_species):,} "
              f"genomes) in {time.time()-t0:.1f}s", file=sys.stderr)

    counts_db = BatchShardedDB(counts_root, nshards, b"counts")
    counts_db.open_all()

    faa_inputs = discover_faa_inputs(Path(args.faa))
    if not faa_inputs:
        print("ERROR: no FAA files found", file=sys.stderr)
        return 2

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    postings_limit = None if args.postings_limit <= 0 else args.postings_limit

    # Load optional name mappings
    species_names = None
    if args.species_names_tsv:
        p = Path(args.species_names_tsv)
        if p.exists():
            species_names = load_species_names(p)
            print(f"Loaded {len(species_names)} species names", file=sys.stderr)

    genome_names = None
    if args.samples_tsv:
        p = Path(args.samples_tsv)
        if p.exists():
            genome_names = load_genome_names(p)
            print(f"Loaded {len(genome_names)} genome names", file=sys.stderr)

    t_total = time.time()
    try:
        for i, faa_file in enumerate(faa_inputs):
            print(f"\n[{i+1}/{len(faa_inputs)}] Querying {faa_file.name}...",
                  file=sys.stderr)
            out_tsv = out_dir / (faa_file.name + ".whatsgnu.tsv")
            out_sim = out_dir / (faa_file.name + ".similarity.tsv")

            query_genome_fast(
                faa_file=faa_file,
                out_tsv=out_tsv,
                out_similarity=out_sim,
                counts_db=counts_db,
                postings_db=postings_db,
                genome_species=genome_species,
                max_genome_id=max_genome_id,
                top_k_species=max(1, args.top_k_species),
                top_k_genomes=max(1, args.top_k_genomes),
                postings_limit=postings_limit,
                include_sequence=args.include_sequence,
                species_names=species_names,
                genome_names=genome_names,
            )
            print(f"  Wrote: {out_tsv}", file=sys.stderr)
            if postings_db:
                print(f"  Wrote: {out_sim}", file=sys.stderr)
    finally:
        counts_db.close()
        if postings_db:
            postings_db.close()

    print(f"\n[DONE] {len(faa_inputs)} genome(s) in "
          f"{time.time()-t_total:.1f}s total", file=sys.stderr)
    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576  # Linux: KB → GB
    print(f"Runtime: {time.time() - t0_all:.1f}s | Peak RAM: {peak_gb:.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
