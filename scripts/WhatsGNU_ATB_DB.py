#!/usr/bin/env python3
"""
WhatsGNU_ATB.py (counts + postings + sequences, stable genome IDs + progress logging)

Builds sharded LMDB databases keyed by 128-bit hash of protein AA sequence.

INPUTS
------
--sample_table (TSV) with columns:
  SampleID (int), Sample (str), SpeciesID (int)
Optional:
  faa_path (str)  # full path to FAA file; if absent uses --faa_dir/<Sample>.faa

--faa_dir (required if faa_path not in sample_table)

OUTPUTS (out_dir)
-----------------
lmdb_counts/ (sharded)
  key: 16-byte hash
  val: func_id:uint32, GNU_count:uint32

lmdb_postings/ (sharded, if --with_postings)
  key: 16-byte hash
  val: n:uint32 + delta+varint encoded sorted unique genome_ids (genome_id = SampleID)

lmdb_sequences/ (sharded, if --with_sequences)
  key: 16-byte hash
  val: representative AA sequence (UTF-8)

metadata/
  functions.tsv.gz
  build_info.json

indexes/
  genome_species.u32  (index by genome_id (=SampleID): species_id)

SUMMARY LOGGING
---------------
- genomes processed / skipped
- total raw protein sequences parsed (across all FAA files)
- total allele-genome records written (unique alleles per genome)
- total unique alleles/hashes (sum of shard keys)

NOTES
-----
GNU_count = number of genomes containing allele at least once (dedup within each FAA).
Function conflicts across genomes: currently "first-seen" func_id wins.
Representative sequence stored: the first AA sequence observed for that hash during reduce
(deterministic due to sorting; derived from the first time we see the allele in build phase).

DEPENDENCIES
------------
pip install lmdb pandas
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import re
import shutil
import struct
import sys
import time
from heapq import merge
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

try:
    import lmdb  # type: ignore
except ImportError:
    lmdb = None

try:
    import pandas as pd  # type: ignore
except ImportError:
    pd = None

try:
    import numpy as np  # type: ignore
except ImportError:
    np = None


# -----------------------------
# FASTA parsing + hashing
# -----------------------------
HEADER_RE = re.compile(r"^>(\S+)\s*(.*)$")

def hash_allele_128(aa_seq: str) -> bytes:
    return hashlib.blake2b(aa_seq.encode("utf-8"), digest_size=16).digest()

def parse_faa(path: Path) -> Iterator[Tuple[str, str]]:
    """
    Yield (aa_sequence, function_string) from a .faa.
    function_string is header description after first whitespace.
    """
    seq_lines: List[str] = []
    func = ""
    with path.open("rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if seq_lines:
                    yield ("".join(seq_lines), func)
                    seq_lines = []
                m = HEADER_RE.match(line)
                func = (m.group(2) if m else "") or ""
                func = func.strip()
            else:
                seq_lines.append(line.strip())
        if seq_lines:
            yield ("".join(seq_lines), func)


# -----------------------------
# Varint encoding (postings)
# -----------------------------
def varint_encode(x: int) -> bytes:
    if x < 0:
        raise ValueError("varint_encode expects non-negative int")
    out = bytearray()
    while True:
        b = x & 0x7F
        x >>= 7
        if x:
            out.append(b | 0x80)
        else:
            out.append(b)
            break
    return bytes(out)

def encode_postings_delta_varint(sorted_unique_ids: List[int]) -> bytes:
    if not sorted_unique_ids:
        return b""
    out = bytearray()
    prev = 0
    for i, gid in enumerate(sorted_unique_ids):
        delta = gid if i == 0 else (gid - prev)
        out += varint_encode(delta)
        prev = gid
    return bytes(out)


# -----------------------------
# Records & sharding
# -----------------------------
# Temp record: 16-byte hash + genome_id:uint32 + func_id:uint32
REC_STRUCT = struct.Struct("<16sII")

# Counts LMDB value
VAL_COUNTS = struct.Struct("<II")  # func_id:uint32, GNU_count:uint32

# Postings LMDB value: n:uint32 + postings_bytes
VAL_POST_HDR = struct.Struct("<I")

def shard_id_from_hash(h16: bytes, nshards: int) -> int:
    return h16[0] & (nshards - 1)


# -----------------------------
# External sort (per shard)
# -----------------------------
def iter_records_from_bin(bin_path: Path) -> Iterator[Tuple[bytes, int, int]]:
    with bin_path.open("rb") as f:
        while True:
            chunk = f.read(REC_STRUCT.size)
            if not chunk:
                break
            yield REC_STRUCT.unpack(chunk)

def write_records_to_bin(records: Iterable[Tuple[bytes, int, int]], out_path: Path) -> None:
    with out_path.open("wb") as f:
        for h, gid, fid in records:
            f.write(REC_STRUCT.pack(h, gid, fid))

def external_sort_bin(in_path: Path, tmp_dir: Path, sort_mem_mb: int) -> Path:
    """External sort records by (hash16, genome_id, func_id).

    IMPORTANT PERFORMANCE NOTE
    --------------------------
    The original tuple-based implementation creates huge Python object overhead,
    which explodes the number of runs and makes merge painfully slow at ATB scale.
    This implementation uses numpy structured arrays (fixed-width 24-byte records)
    to keep RAM usage close to sort_mem_mb and reduce run count.

    The output is a single sorted binary file with the same record packing as REC_STRUCT.
    """
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # If numpy isn't available, fall back to the slower pure-Python sort.
    if np is None:
        recs_per_chunk = max(1, (sort_mem_mb * 1024 * 1024) // REC_STRUCT.size)
        runs: List[Path] = []
        buf: List[Tuple[bytes, int, int]] = []

        def flush_run(local_buf: List[Tuple[bytes, int, int]]) -> None:
            local_buf.sort(key=lambda x: (x[0], x[1], x[2]))
            run_path = tmp_dir / f"run_{len(runs):06d}.bin"
            write_records_to_bin(local_buf, run_path)
            runs.append(run_path)

        for rec in iter_records_from_bin(in_path):
            buf.append(rec)
            if len(buf) >= recs_per_chunk:
                flush_run(buf)
                buf = []
        if buf:
            flush_run(buf)

        if len(runs) == 1:
            return runs[0]

        iters = [iter_records_from_bin(p) for p in runs]
        merged_iter = merge(*iters, key=lambda x: (x[0], x[1], x[2]))

        sorted_path = tmp_dir / "sorted_merged.bin"
        write_records_to_bin(merged_iter, sorted_path)

        for p in runs:
            try:
                p.unlink()
            except OSError:
                pass
        return sorted_path

    # Numpy fast path
    dtype = np.dtype([("h", "S16"), ("gid", "<u4"), ("fid", "<u4")], align=False)
    rec_size = dtype.itemsize  # 24
    assert rec_size == REC_STRUCT.size

    # Reserve ~90% of sort_mem_mb for the numpy array itself.
    mem_bytes = int(sort_mem_mb) * 1024 * 1024
    target_bytes = max(rec_size, int(mem_bytes * 0.90))
    recs_per_chunk = max(1, target_bytes // rec_size)

    runs: List[Path] = []

    with in_path.open("rb") as f:
        run_idx = 0
        while True:
            raw = f.read(recs_per_chunk * rec_size)
            if not raw:
                break
            # Truncate to full records
            raw = raw[: (len(raw) // rec_size) * rec_size]
            arr = np.frombuffer(raw, dtype=dtype)
            # lexsort uses last key as primary -> provide (fid, gid, h) for primary h
            order = np.lexsort((arr["fid"], arr["gid"], arr["h"]))
            arr_sorted = arr[order]
            run_path = tmp_dir / f"run_{run_idx:06d}.bin"
            with run_path.open("wb") as out:
                out.write(arr_sorted.tobytes(order="C"))
            runs.append(run_path)
            run_idx += 1

    if not runs:
        # empty input
        outp = tmp_dir / "sorted_merged.bin"
        outp.write_bytes(b"")
        return outp

    if len(runs) == 1:
        return runs[0]

    def merge_some(run_paths: List[Path], out_path: Path) -> None:
        iters = [iter_records_from_bin(p) for p in run_paths]
        merged_iter = merge(*iters, key=lambda x: (x[0], x[1], x[2]))
        write_records_to_bin(merged_iter, out_path)

    # Batched multi-pass merge to avoid thousands of open files.
    MERGE_FANIN = 64
    round_idx = 0
    cur_runs = runs
    while len(cur_runs) > 1:
        next_runs: List[Path] = []
        for i in range(0, len(cur_runs), MERGE_FANIN):
            chunk = cur_runs[i : i + MERGE_FANIN]
            if len(chunk) == 1:
                next_runs.append(chunk[0])
                continue
            out_run = tmp_dir / f"merge_r{round_idx:02d}_{len(next_runs):06d}.bin"
            merge_some(chunk, out_run)
            # cleanup merged inputs
            for p in chunk:
                try:
                    p.unlink()
                except OSError:
                    pass
            next_runs.append(out_run)
        cur_runs = next_runs
        round_idx += 1

    final_sorted = tmp_dir / "sorted_merged.bin"
    # Move final run to canonical name
    if cur_runs[0] != final_sorted:
        try:
            cur_runs[0].replace(final_sorted)
        except OSError:
            shutil.copy2(cur_runs[0], final_sorted)
            try:
                cur_runs[0].unlink()
            except OSError:
                pass

    return final_sorted


# -----------------------------
# LMDB helpers
# -----------------------------
def lmdb_open_env(path: Path, map_size_bytes: int, readonly: bool) -> "lmdb.Environment":
    assert lmdb is not None
    path.mkdir(parents=True, exist_ok=True)
    return lmdb.open(
        str(path),
        map_size=map_size_bytes,
        subdir=True,
        create=not readonly,
        readonly=readonly,
        lock=not readonly,
        readahead=readonly,
        writemap=False,
        metasync=False,
        sync=False,
        max_dbs=1,
    )

def reduce_sorted_records(
    sorted_bin: Path,
    out_counts_shard: Path,
    out_postings_shard: Optional[Path],
    out_sequences_shard: Optional[Path],
    map_size_bytes_counts: int,
    map_size_bytes_postings: int,
    map_size_bytes_sequences: int,
    with_postings: bool,
    with_sequences: bool,
) -> Tuple[int, int]:
    """
    Consume sorted records (hash, genome_id, func_id) and write:
      counts:    hash -> (func_id, GNU_count)
      postings:  hash -> n:uint32 + delta+varint(genome_ids)   (optional)
      sequences: hash -> representative AA sequence UTF-8       (optional; stored at reduce time)

    IMPORTANT:
    - Postings list is unique+sorted due to sorted input + per-allele unique genome_id counting.
    - Sequences are stored from a per-shard hash->sequence map produced during parse.
    """
    assert lmdb is not None

    env_c = lmdb_open_env(out_counts_shard, map_size_bytes_counts, readonly=False)
    db_c = env_c.open_db(b"counts")

    env_p = None
    db_p = None
    if with_postings and out_postings_shard is not None:
        env_p = lmdb_open_env(out_postings_shard, map_size_bytes_postings, readonly=False)
        db_p = env_p.open_db(b"postings")

    env_s = None
    db_s = None
    if with_sequences and out_sequences_shard is not None:
        env_s = lmdb_open_env(out_sequences_shard, map_size_bytes_sequences, readonly=False)
        db_s = env_s.open_db(b"sequences")

    n_keys = 0
    n_records = 0

    current_h: Optional[bytes] = None
    current_fid: int = 0
    last_gid: Optional[int] = None
    gnu_count: int = 0
    postings: List[int] = []

    txn_c = env_c.begin(write=True, db=db_c)
    txn_p = env_p.begin(write=True, db=db_p) if env_p is not None else None
    txn_s = env_s.begin(write=True, db=db_s) if env_s is not None else None

    # We will not have direct access to sequences in this reducer unless we provide a map.
    # The main() function will provide a per-shard sequences "side file" for lookup.
    # That file stores: hash16 + seq_len:uint32 + seq_bytes
    # We'll read it into a dict iterator and keep in sync with allele order (sorted by hash).
    seq_iter: Optional[Iterator[Tuple[bytes, str]]] = None
    seq_next: Optional[Tuple[bytes, str]] = None

    # Load sequence side file if requested
    if with_sequences:
        seq_side = sorted_bin.parent / "seq_side.bin"
        if not seq_side.exists():
            raise RuntimeError(f"with_sequences requested but missing sequence side file: {seq_side}")

        def iter_seq_side(p: Path) -> Iterator[Tuple[bytes, str]]:
            with p.open("rb") as f:
                while True:
                    h = f.read(16)
                    if not h:
                        break
                    ln_bytes = f.read(4)
                    if len(ln_bytes) != 4:
                        raise IOError("Corrupt seq_side.bin (length)")
                    (ln,) = struct.unpack("<I", ln_bytes)
                    seqb = f.read(ln)
                    if len(seqb) != ln:
                        raise IOError("Corrupt seq_side.bin (payload)")
                    yield (h, seqb.decode("utf-8"))

        seq_iter = iter_seq_side(seq_side)
        try:
            seq_next = next(seq_iter)
        except StopIteration:
            seq_next = None

    def get_seq_for_hash(h: bytes) -> Optional[str]:
        # seq_side.bin is sorted by hash (we will generate it that way), so we can advance
        nonlocal seq_next, seq_iter
        if not with_sequences or seq_iter is None:
            return None
        while seq_next is not None and seq_next[0] < h:
            try:
                seq_next = next(seq_iter)
            except StopIteration:
                seq_next = None
                break
        if seq_next is not None and seq_next[0] == h:
            return seq_next[1]
        return None

    def flush_current() -> None:
        nonlocal n_keys, postings, gnu_count, current_h, current_fid, txn_c, txn_p, txn_s
        if current_h is None:
            return
        txn_c.put(current_h, VAL_COUNTS.pack(current_fid, gnu_count))
        if with_postings and txn_p is not None:
            pb = encode_postings_delta_varint(postings)
            txn_p.put(current_h, VAL_POST_HDR.pack(len(postings)) + pb)
        if with_sequences and txn_s is not None:
            seq = get_seq_for_hash(current_h)
            if seq is not None:
                # Put only if absent to be safe; duplicates shouldn't exist.
                txn_s.put(current_h, seq.encode("utf-8"), overwrite=False)
        n_keys += 1

    for h, gid, fid in iter_records_from_bin(sorted_bin):
        n_records += 1
        if current_h is None:
            current_h = h
            current_fid = fid
            last_gid = gid
            gnu_count = 1
            postings = [gid] if with_postings else []
            continue

        if h == current_h:
            if gid != last_gid:
                gnu_count += 1
                last_gid = gid
                if with_postings:
                    postings.append(gid)
        else:
            flush_current()
            current_h = h
            current_fid = fid
            last_gid = gid
            gnu_count = 1
            postings = [gid] if with_postings else []

        if n_keys > 0 and (n_keys % 5_000_000 == 0):
            txn_c.commit()
            txn_c = env_c.begin(write=True, db=db_c)
            if txn_p is not None:
                txn_p.commit()
                txn_p = env_p.begin(write=True, db=db_p)
            if txn_s is not None:
                txn_s.commit()
                txn_s = env_s.begin(write=True, db=db_s)

    flush_current()

    txn_c.commit()
    if txn_p is not None:
        txn_p.commit()
    if txn_s is not None:
        txn_s.commit()

    env_c.sync()
    env_c.close()
    if env_p is not None:
        env_p.sync()
        env_p.close()
    if env_s is not None:
        env_s.sync()
        env_s.close()

    return n_keys, n_records


# -----------------------------
# Export allele counts TSV
# -----------------------------
def export_allele_counts_tsv(counts_root: Path, nshards: int, out_tsv: Path) -> None:
    assert lmdb is not None
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("wt", encoding="utf-8") as out:
        out.write("allele_hash\tGNU_score\tfunc_id\n")
        for sid in range(nshards):
            shard_dir = counts_root / f"shard_{sid:02x}"
            if not shard_dir.exists():
                continue

            env = lmdb_open_env(shard_dir, map_size_bytes=1 << 30, readonly=True)
            db = env.open_db(b"counts")              # <-- THIS is the key fix
            with env.begin(db=db) as txn:            # <-- iterate named DB
                cur = txn.cursor()
                for k, v in cur:
                    fid, gnu = VAL_COUNTS.unpack(v)  # v is now guaranteed 8 bytes
                    out.write(f"{k.hex()}\t{gnu}\t{fid}\n")
            env.close()

# -----------------------------
# Logging + progress
# -----------------------------
def setup_logger(log_path: Path, level: str = "INFO") -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("WhatsGNU_ATB")
    logger.setLevel(getattr(logging, level))
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")

    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    logger.propagate = False
    return logger

def should_report(genomes_done: int) -> bool:
    if genomes_done == 10_000:
        return True
    if genomes_done == 100_000:
        return True
    if genomes_done > 100_000 and (genomes_done % 100_000 == 0):
        return True
    return False


# -----------------------------
# Misc helpers
# -----------------------------
def write_gzip_tsv(path: Path, rows: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as gz:
        for r in rows:
            gz.write(r)
            if not r.endswith("\n"):
                gz.write("\n")


def have_all_rec_files(tmp_dir: Path, shards: int) -> bool:
    """Return True if all rec_XX.bin files exist and are non-empty."""
    for sid in range(shards):
        p = tmp_dir / f"rec_{sid:02x}.bin"
        if (not p.exists()) or p.stat().st_size == 0:
            return False
    return True


def main() -> int:
    ap = argparse.ArgumentParser()

    ap.add_argument("--sample_table", required=True, help="TSV with SampleID, Sample, SpeciesID (+ optional faa_path).")
    ap.add_argument("--faa_dir", default=None, help="Directory holding <Sample>.faa (used if faa_path not in table).")
    ap.add_argument("--out_dir", required=True, help="Output directory.")
    ap.add_argument("--tmp_dir", default=None, help="Temp dir (default: out_dir/tmp).")

    ap.add_argument("--reduce_tmp_dir", default=None,
                    help="Optional local scratch directory for reduce/sort work (e.g. /tmp/WGNU_$SLURM_JOB_ID). "
                         "If provided, each shard's rec file (and sequence side file if enabled) is copied here and "
                         "sorting/merging happens locally for speed. Outputs still go to out_dir.")

    ap.add_argument("--parse_only", action="store_true",
                    help="Only parse FAA files and write rec_*.bin (+ functions.tsv.gz, genome_species.u32). Do not reduce.")
    ap.add_argument("--reduce_only", action="store_true",
                    help="Only reduce existing rec_*.bin into LMDB. Skip FAA parsing.")
    ap.add_argument("--resume", action="store_true",
                    help="If rec_*.bin already exist in tmp_dir, skip parsing and start reduce; otherwise parse then reduce.")
    ap.add_argument("--skip_existing_shards", action="store_true",
                    help="When reducing, skip shards whose LMDB output already exists (useful for resuming after timeout).")

    ap.add_argument("--shards", type=int, default=16, help="Number of shards (power of 2). Default=16.")
    ap.add_argument("--sort_mem_mb", type=int, default=65536,
                    help="RAM for external sort run creation per shard (MB). Default=65536.")

    ap.add_argument("--lmdb_map_gb_counts_per_shard", type=int, default=24,
                    help="LMDB map size per shard for counts DB (GB). Default=24.")
    ap.add_argument("--lmdb_map_gb_postings_per_shard", type=int, default=160,
                    help="LMDB map size per shard for postings DB (GB). Default=160.")
    ap.add_argument("--lmdb_map_gb_sequences_per_shard", type=int, default=25,
                    help="LMDB map size per shard for sequences DB (GB). Default=25.")

    ap.add_argument("--with_postings", action="store_true", help="Also write postings LMDB.")
    ap.add_argument("--with_sequences", action="store_true", help="Also write sequences LMDB (hash->AA seq).")
    ap.add_argument("--faa_suffix", default=".bakta.faa",
                help="FAA filename suffix appended to Sample (default: .faa). Example: .bakta.faa")

    ap.add_argument("--export_allele_counts", default=None,
                    help="Optional path to write allele_hash + GNU_score (+ func_id) TSV at the end.")

    ap.add_argument("--log_file", default=None, help="Progress log file (default: <out_dir>/build.log).")
    ap.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = ap.parse_args()

    if lmdb is None:
        print("ERROR: missing lmdb. Install: pip install lmdb", file=sys.stderr)
        return 2
    if pd is None:
        print("ERROR: missing pandas. Install: pip install pandas", file=sys.stderr)
        return 2
    if np is None:
        print("ERROR: missing numpy. Install: pip install numpy", file=sys.stderr)
        return 2

    nshards = args.shards
    if nshards & (nshards - 1) != 0:
        print("ERROR: --shards must be a power of 2 (e.g., 32, 64).", file=sys.stderr)
        return 2

    out_dir = Path(args.out_dir)
    tmp_dir = Path(args.tmp_dir) if args.tmp_dir else (out_dir / "tmp")
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Phase control
    phase_flags = int(bool(args.parse_only)) + int(bool(args.reduce_only)) + int(bool(args.resume))
    if phase_flags > 1:
        print("ERROR: Use at most one of --parse_only, --reduce_only, --resume.", file=sys.stderr)
        return 2

    rec_ready = have_all_rec_files(tmp_dir, args.shards)
    if args.reduce_only and not rec_ready:
        print(f"ERROR: --reduce_only set but rec_*.bin not found (or empty) under tmp_dir={tmp_dir}", file=sys.stderr)
        return 2

    if args.reduce_only:
        do_parse, do_reduce = False, True
    elif args.parse_only:
        do_parse, do_reduce = True, False
    elif args.resume:
        do_parse, do_reduce = (not rec_ready), True
    else:
        do_parse, do_reduce = True, True

    reduce_tmp_root = Path(args.reduce_tmp_dir) if args.reduce_tmp_dir else None
    if reduce_tmp_root:
        reduce_tmp_root.mkdir(parents=True, exist_ok=True)

    log_path = Path(args.log_file) if args.log_file else (out_dir / "build.log")
    logger = setup_logger(log_path, args.log_level)

    start_time = time.time()
    last_time = start_time
    last_genomes = 0
    last_records = 0

    logger.info(
        f"Starting build | shards={nshards} | postings={args.with_postings} | sequences={args.with_sequences} "
        f"| sort_mem_mb={args.sort_mem_mb} | parse={do_parse} | reduce={do_reduce}"
    )
    logger.info(f"Output dir: {out_dir}")
    logger.info(f"Temp dir: {tmp_dir}")
    if reduce_tmp_root:
        logger.info(f"Reduce tmp dir: {reduce_tmp_root} (local scratch)")
    logger.info(f"Sample table: {args.sample_table}")
    if args.faa_dir:
        logger.info(f"FAA dir: {args.faa_dir}")
    logger.info(
        f"LMDB map sizes per shard (GB): counts={args.lmdb_map_gb_counts_per_shard}, "
        f"postings={args.lmdb_map_gb_postings_per_shard}, sequences={args.lmdb_map_gb_sequences_per_shard}"
    )

    sample_table = Path(args.sample_table)
    df = pd.read_csv(sample_table, sep="\t", dtype=str)

    required_cols = {"SampleID", "Sample", "SpeciesID"}
    missing = required_cols - set(df.columns)
    if missing:
        logger.error(f"sample_table missing columns: {sorted(missing)}")
        return 2

    df["SampleID"] = df["SampleID"].astype(int)
    df["SpeciesID"] = df["SpeciesID"].astype(int)
    df["Sample"] = df["Sample"].astype(str).str.strip()

    if "faa_path" in df.columns:
        df["faa_path"] = df["faa_path"].astype(str)
    else:
        if not args.faa_dir:
            logger.error("Provide --faa_dir if sample_table lacks faa_path column.")
            return 2
        faa_dir = Path(args.faa_dir)
        df["faa_path"] = df["Sample"].apply(lambda s: str(faa_dir / f"{s}{args.faa_suffix}"))

    df = df.sort_values("SampleID").reset_index(drop=True)

    # Build genome_species.u32 indexed by genome_id (=SampleID)
    idx_dir = out_dir / "indexes"
    idx_dir.mkdir(parents=True, exist_ok=True)
    max_gid = int(df["SampleID"].max())
    genome_species_path = idx_dir / "genome_species.u32"

    if do_parse or (not genome_species_path.exists()) or genome_species_path.stat().st_size == 0:
        logger.info(f"Writing genome_species.u32 for max genome_id={max_gid:,}")
        with genome_species_path.open("wb") as f:
            f.write(struct.pack("<I", 0))  # index 0 unused
            for _ in range(max_gid):
                f.write(struct.pack("<I", 0))
        with genome_species_path.open("r+b") as f:
            for _, row in df.iterrows():
                gid = int(row["SampleID"])
                sid = int(row["SpeciesID"])
                f.seek(gid * 4)
                f.write(struct.pack("<I", sid))
    else:
        logger.info(f"Using existing genome_species.u32: {genome_species_path}")

    # Record bins (always defined; may have been produced earlier if --resume/--reduce_only)
    shard_bins = [tmp_dir / f"rec_{i:02x}.bin" for i in range(nshards)]
    seq_side_bins = [tmp_dir / f"seq_{i:02x}.bin" for i in range(nshards)]

    # Parse phase (FAA -> rec_*.bin + functions.tsv.gz)
    func_list: List[str] = []
    total_records_written = 0              # allele-genome records (after per-genome dedup)
    total_protein_seqs_seen = 0            # raw proteins parsed (including duplicates within genome)
    n_genomes_processed = 0
    n_genomes_skipped = 0

    meta_dir = out_dir / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)
    functions_path = meta_dir / "functions.tsv.gz"

    if do_parse:
        # Start fresh rec files in tmp_dir for this run
        for p in shard_bins:
            try:
                p.unlink()
            except OSError:
                pass
        if args.with_sequences:
            for p in seq_side_bins:
                try:
                    p.unlink()
                except OSError:
                    pass

        shard_fhs = [open(p, "wb", buffering=16 * 1024 * 1024) for p in shard_bins]
        seq_side_fhs = [open(p, "wb", buffering=16 * 1024 * 1024) for p in seq_side_bins] if args.with_sequences else []

        func_to_id: Dict[str, int] = {}

        def get_func_id(func: str) -> int:
            f = (func or "").strip()
            if f not in func_to_id:
                func_to_id[f] = len(func_list)
                func_list.append(f)
            return func_to_id[f]

        logger.info(f"Processing {len(df):,} genomes (rows in sample_table)")

        for _, row in df.iterrows():
            gid = int(row["SampleID"])
            sample = str(row["Sample"])
            faa_path = Path(str(row["faa_path"]))

            if not faa_path.exists():
                n_genomes_skipped += 1
                if n_genomes_skipped <= 20:
                    logger.warning(f"Missing FAA for {sample} at {faa_path} (skipping)")
                continue
            if faa_path.name.endswith(".gz"):
                n_genomes_skipped += 1
                if n_genomes_skipped <= 20:
                    logger.warning(f"Compressed FAA not supported: {faa_path} (decompress first; skipping)")
                continue

            # per-genome dedup: hash -> func_id
            seen: Dict[bytes, int] = {}
            seen_seq: Dict[bytes, str] = {} if args.with_sequences else {}

            for aa_seq, func in parse_faa(faa_path):
                if not aa_seq:
                    continue
                total_protein_seqs_seen += 1
                h = hash_allele_128(aa_seq)
                if h in seen:
                    continue
                fid = get_func_id(func)
                seen[h] = fid
                if args.with_sequences:
                    seen_seq[h] = aa_seq

            for h, fid in seen.items():
                sid_shard = shard_id_from_hash(h, nshards)
                shard_fhs[sid_shard].write(REC_STRUCT.pack(h, gid, fid))
                total_records_written += 1

            if args.with_sequences:
                for h, seq in seen_seq.items():
                    sid_shard = shard_id_from_hash(h, nshards)
                    sb = seq.encode("utf-8")
                    seq_side_fhs[sid_shard].write(h)
                    seq_side_fhs[sid_shard].write(struct.pack("<I", len(sb)))
                    seq_side_fhs[sid_shard].write(sb)

            n_genomes_processed += 1

            if should_report(n_genomes_processed):
                now = time.time()
                elapsed = now - start_time
                interval = now - last_time if now > last_time else 1e-9

                dG = n_genomes_processed - last_genomes
                dR = total_records_written - last_records

                g_per_s = dG / interval
                r_per_s = dR / interval

                logger.info(
                    f"Progress genomes={n_genomes_processed:,} "
                    f"skipped={n_genomes_skipped:,} "
                    f"proteins_raw={total_protein_seqs_seen:,} "
                    f"allele_genome_records={total_records_written:,} "
                    f"elapsed={elapsed/3600:.2f}h "
                    f"rate={g_per_s:,.2f} genomes/s "
                    f"({g_per_s*3600:,.0f} genomes/h) "
                    f"rec_rate={r_per_s:,.0f} rec/s"
                )

                last_time = now
                last_genomes = n_genomes_processed
                last_records = total_records_written

        for fh in shard_fhs:
            fh.close()
        if args.with_sequences:
            for fh in seq_side_fhs:
                fh.close()

        logger.info(
            f"Finished parsing genomes | processed={n_genomes_processed:,} skipped={n_genomes_skipped:,} "
            f"proteins_raw={total_protein_seqs_seen:,} allele_genome_records={total_records_written:,}"
        )

        # Write functions metadata
        write_gzip_tsv(
            functions_path,
            ["func_id\tfunction"] + [f"{i}\t{func_list[i]}" for i in range(len(func_list))]
        )
        logger.info(f"Wrote functions.tsv.gz with {len(func_list):,} unique functions")
        n_functions_val = len(func_list)
    else:
        logger.info("Skipping FAA parsing (reduce-only/resume).")
        if functions_path.exists():
            # Count functions without loading whole file into memory
            n_functions = 0
            with gzip.open(functions_path, "rt", encoding="utf-8", errors="replace") as gz:
                for i, _line in enumerate(gz):
                    pass
            # i is last index; subtract header (line 0)
            try:
                n_functions = max(0, i)
            except UnboundLocalError:
                n_functions = 0
            logger.info(f"Using existing functions.tsv.gz with {n_functions:,} unique functions")
            # Keep func_list empty to avoid loading; build_info will use n_functions
            n_functions_val = n_functions
        else:
            logger.warning("functions.tsv.gz not found; build_info will report n_functions=0")
            n_functions_val = 0

    if not do_reduce:
        # Parse-only mode: write a minimal build_info.json and exit.
        build_info = {
            "version": "v5",
            "mode": "parse_only",
            "hash": "blake2b_128(digest_size=16)",
            "shards": nshards,
            "shard_rule": "shard_id = first_byte(hash16) & (shards-1)",
            "gnu_count_definition": "number_of_genomes_containing_allele (dedup within each .faa)",
            "sort_mem_mb": args.sort_mem_mb,
            "n_functions": n_functions_val,
            "n_genomes_processed": n_genomes_processed,
            "n_genomes_skipped": n_genomes_skipped,
            "proteins_raw": total_protein_seqs_seen,
            "allele_genome_records": total_records_written,
            "unique_alleles": 0,
            "per_shard": [],
            "note": "Reduce not executed (--parse_only).",
        }
        (meta_dir / "build_info.json").write_text(json.dumps(build_info, indent=2), encoding="utf-8")
        logger.info("Parse-only complete; wrote build_info.json and exiting.")
        return 0

    # Build LMDB shards
    counts_root = out_dir / "lmdb_counts"
    postings_root = out_dir / "lmdb_postings"
    sequences_root = out_dir / "lmdb_sequences"

    counts_root.mkdir(parents=True, exist_ok=True)
    if args.with_postings:
        postings_root.mkdir(parents=True, exist_ok=True)
    if args.with_sequences:
        sequences_root.mkdir(parents=True, exist_ok=True)

    map_counts = int(args.lmdb_map_gb_counts_per_shard) * (1024 ** 3)
    map_posts = int(args.lmdb_map_gb_postings_per_shard) * (1024 ** 3)
    map_seqs = int(args.lmdb_map_gb_sequences_per_shard) * (1024 ** 3)

    per_shard_stats = []
    total_unique_alleles = 0

    for sid in range(nshards):
        bin_path = shard_bins[sid]
        if not bin_path.exists() or bin_path.stat().st_size == 0:
            logger.info(f"[reduce] shard {sid:02x}: empty, skipping")
            continue

        out_counts_shard = counts_root / f"shard_{sid:02x}"
        out_posts_shard = (postings_root / f"shard_{sid:02x}") if args.with_postings else None
        out_seqs_shard = (sequences_root / f"shard_{sid:02x}") if args.with_sequences else None

        if args.skip_existing_shards:
            # If counts shard looks present, assume shard completed.
            if (out_counts_shard / "data.mdb").exists() and (out_counts_shard / "lock.mdb").exists():
                logger.info(f"[reduce] shard {sid:02x}: output exists, skipping (--skip_existing_shards)")
                continue

        work_root = reduce_tmp_root if reduce_tmp_root else tmp_dir
        shard_tmp = work_root / f"shard_{sid:02x}"
        if shard_tmp.exists():
            shutil.rmtree(shard_tmp)
        shard_tmp.mkdir(parents=True, exist_ok=True)

        # If using a separate local work_root, copy the rec file (and seq side file) locally for faster I/O.
        rec_for_sort = bin_path
        if work_root != tmp_dir:
            local_rec = shard_tmp / bin_path.name
            if (not local_rec.exists()) or (local_rec.stat().st_size != bin_path.stat().st_size):
                shutil.copy2(bin_path, local_rec)
            rec_for_sort = local_rec

        # If sequences enabled, we need to sort+dedup seq_side.bin by hash within this shard
        # and write shard_tmp/seq_side.bin (sorted unique) for the reducer.
        if args.with_sequences:
            seq_in = seq_side_bins[sid]
            if work_root != tmp_dir:
                local_seq = shard_tmp / seq_in.name
                if seq_in.exists() and ((not local_seq.exists()) or (local_seq.stat().st_size != seq_in.stat().st_size)):
                    shutil.copy2(seq_in, local_seq)
                seq_in = local_seq
            if not seq_in.exists() or seq_in.stat().st_size == 0:
                logger.warning(f"[reduce] shard {sid:02x}: sequences enabled but seq_side empty")
            else:
                # Convert variable-length records to in-memory chunks then external sort by hash.
                # We do a simple external sort into runs of (hash, seq) using Python; mem bound by sort_mem_mb.
                logger.info(f"[reduce] shard {sid:02x}: sorting+deduping sequences side file")
                t_seq = time.time()

                # Read seq records into runs
                recs_per_chunk = max(1, (args.sort_mem_mb * 1024 * 1024) // 512)  # heuristic
                runs: List[Path] = []
                buf: List[Tuple[bytes, bytes]] = []

                def flush_seq_run(local_buf: List[Tuple[bytes, bytes]]) -> None:
                    local_buf.sort(key=lambda x: x[0])
                    # dedup within run (keep first)
                    deduped: List[Tuple[bytes, bytes]] = []
                    last_h: Optional[bytes] = None
                    for h, sb in local_buf:
                        if last_h is None or h != last_h:
                            deduped.append((h, sb))
                            last_h = h
                    run_path = shard_tmp / f"seq_run_{len(runs):06d}.bin"
                    with run_path.open("wb") as f:
                        for h, sb in deduped:
                            f.write(h)
                            f.write(struct.pack("<I", len(sb)))
                            f.write(sb)
                    runs.append(run_path)

                with seq_in.open("rb") as f:
                    while True:
                        h = f.read(16)
                        if not h:
                            break
                        ln_bytes = f.read(4)
                        if len(ln_bytes) != 4:
                            raise RuntimeError(f"Corrupt sequence side file: {seq_in}")
                        (ln,) = struct.unpack("<I", ln_bytes)
                        sb = f.read(ln)
                        if len(sb) != ln:
                            raise RuntimeError(f"Corrupt sequence side file (payload): {seq_in}")
                        buf.append((h, sb))
                        if len(buf) >= recs_per_chunk:
                            flush_seq_run(buf)
                            buf = []
                    if buf:
                        flush_seq_run(buf)
                        buf = []

                # Merge runs -> shard_tmp/seq_side.bin with global dedup
                def iter_seq_run(p: Path) -> Iterator[Tuple[bytes, bytes]]:
                    with p.open("rb") as f:
                        while True:
                            h = f.read(16)
                            if not h:
                                break
                            ln_bytes = f.read(4)
                            (ln,) = struct.unpack("<I", ln_bytes)
                            sb = f.read(ln)
                            yield (h, sb)

                if len(runs) == 0:
                    # no sequences for this shard
                    pass
                elif len(runs) == 1:
                    # already sorted+deduped
                    (shard_tmp / "seq_side.bin").write_bytes(runs[0].read_bytes())
                else:
                    iters = [iter_seq_run(p) for p in runs]
                    merged = merge(*iters, key=lambda x: x[0])
                    out_seq = shard_tmp / "seq_side.bin"
                    with out_seq.open("wb") as out:
                        last_h: Optional[bytes] = None
                        for h, sb in merged:
                            if last_h is not None and h == last_h:
                                continue
                            out.write(h)
                            out.write(struct.pack("<I", len(sb)))
                            out.write(sb)
                            last_h = h

                # cleanup seq runs
                for p in runs:
                    try:
                        p.unlink()
                    except OSError:
                        pass

                logger.info(f"[reduce] shard {sid:02x}: sequences side ready in {(time.time()-t_seq)/60:.1f} min")

        logger.info(f"[reduce] shard {sid:02x}: external sort allele records (mem={args.sort_mem_mb} MB)")
        t0 = time.time()
        sorted_bin = external_sort_bin(rec_for_sort, shard_tmp, sort_mem_mb=args.sort_mem_mb)
        logger.info(f"[reduce] shard {sid:02x}: sort done in {(time.time()-t0)/60:.1f} min")

        logger.info(f"[reduce] shard {sid:02x}: writing LMDB (counts{' + postings' if args.with_postings else ''}{' + sequences' if args.with_sequences else ''})")
        t1 = time.time()
        n_keys, n_recs = reduce_sorted_records(
            sorted_bin=sorted_bin,
            out_counts_shard=out_counts_shard,
            out_postings_shard=out_posts_shard,
            out_sequences_shard=out_seqs_shard,
            map_size_bytes_counts=map_counts,
            map_size_bytes_postings=map_posts,
            map_size_bytes_sequences=map_seqs,
            with_postings=args.with_postings,
            with_sequences=args.with_sequences,
        )
        logger.info(f"[reduce] shard {sid:02x}: LMDB write done in {(time.time()-t1)/60:.1f} min | keys={n_keys:,} records={n_recs:,}")

        per_shard_stats.append({
            "shard": f"{sid:02x}",
            "records": n_recs,
            "keys": n_keys,
            "sorted_bytes": sorted_bin.stat().st_size if sorted_bin.exists() else None,
        })

        total_unique_alleles += n_keys
        shutil.rmtree(shard_tmp, ignore_errors=True)

    # Final summary numbers
    logger.info(
        f"SUMMARY | genomes_processed={n_genomes_processed:,} genomes_skipped={n_genomes_skipped:,} "
        f"proteins_raw={total_protein_seqs_seen:,} allele_genome_records={total_records_written:,} "
        f"unique_alleles={total_unique_alleles:,}"
    )

    build_info = {
        "version": "v5",
        "mode": "counts_postings_sequences" if (args.with_postings and args.with_sequences)
                else ("counts_and_postings" if args.with_postings else ("counts_and_sequences" if args.with_sequences else "counts_only")),
        "hash": "blake2b_128(digest_size=16)",
        "shards": nshards,
        "shard_rule": "shard_id = first_byte(hash16) & (shards-1)",
        "gnu_count_definition": "number_of_genomes_containing_allele (dedup within each .faa)",
        "postings_definition": "sorted unique genome_ids (SampleID) encoded as delta+varint" if args.with_postings else None,
        "sequences_definition": "representative AA sequence stored for each allele hash (UTF-8)" if args.with_sequences else None,
        "sort_mem_mb": args.sort_mem_mb,
        "lmdb_value_counts": "func_id:uint32, GNU_count:uint32 (little-endian)",
        "lmdb_value_postings": "n:uint32 + delta+varint genome_ids" if args.with_postings else None,
        "lmdb_value_sequences": "utf-8 AA sequence" if args.with_sequences else None,
        "n_functions": n_functions_val,
        "n_genomes_processed": n_genomes_processed,
        "n_genomes_skipped": n_genomes_skipped,
        "proteins_raw": total_protein_seqs_seen,
        "allele_genome_records": total_records_written,
        "unique_alleles": total_unique_alleles,
        "per_shard": per_shard_stats,
        "reduce_tmp_dir": str(reduce_tmp_root) if reduce_tmp_root else None,
        "parse_executed": do_parse,
        "reduce_executed": do_reduce,
    }
    (meta_dir / "build_info.json").write_text(json.dumps(build_info, indent=2), encoding="utf-8")
    logger.info("Wrote build_info.json")

    if args.export_allele_counts:
        export_path = Path(args.export_allele_counts)
        logger.info(f"[export] writing allele counts TSV to {export_path}")
        export_allele_counts_tsv(counts_root, nshards, export_path)
        logger.info("[export] done")

    total_elapsed_h = (time.time() - start_time) / 3600
    logger.info(f"Build complete in {total_elapsed_h:.2f} hours")
    logger.info(f"Output: {out_dir}")
    logger.info(f"Log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
