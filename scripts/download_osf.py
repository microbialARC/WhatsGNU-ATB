#!/usr/bin/env python3
"""
download_osf.py — Download specific folders from the WhatsGNU-ATB OSF project.

Usage:
    # Download the database (required for querying)
    python download_osf.py --folder WGNU_ATB_DB --out-dir ./whatsgnu_db

    # Download sample tables
    python download_osf.py --folder Sample_tables --out-dir ./whatsgnu_db

    # Download multiple folders
    python download_osf.py --folder WGNU_ATB_DB Sample_tables --out-dir ./whatsgnu_db

    # Download everything
    python download_osf.py --all --out-dir ./whatsgnu_db

    # List available folders
    python download_osf.py --list

Available folders:
    WGNU_ATB_DB                  Pre-built LMDB database (required for querying)
    Sample_tables                Sample IDs, species mapping, genome lists
    ATB_hash_seq                 Hash-to-sequence lookup (20 compressed parts)
    ATB_summary_figures_tables   Publication figures, tables, and cache
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from typing import List, Optional


OSF_PROJECT = "6jr4u"
API_BASE = "https://api.osf.io/v2"

FOLDER_DESCRIPTIONS = {
    "WGNU_ATB_DB": "Pre-built LMDB database (8 count + 8 posting shards, indexes, metadata, Sample-to-ID mapping). Required for querying.",
    "Sample_tables": "Species stats, genome lists.",
    "ATB_hash_seq": "Hash-to-amino-acid-sequence lookup table (20 xz-compressed parts).",
    "ATB_summary_figures_tables": "Publication figures, per-species tables, allele analysis results, and counts cache.",
}


def api_get(url: str, token: Optional[str] = None) -> dict:
    """GET request to OSF API, returns parsed JSON."""
    req = urllib.request.Request(url)
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read().decode("utf-8"))


def download_file(url: str, dest: Path, token: Optional[str] = None) -> None:
    """Download a file with progress reporting."""
    req = urllib.request.Request(url)
    if token:
        req.add_header("Authorization", f"Bearer {token}")

    dest.parent.mkdir(parents=True, exist_ok=True)

    with urllib.request.urlopen(req) as resp:
        total = resp.headers.get("Content-Length")
        total = int(total) if total else None
        downloaded = 0
        block_size = 1024 * 1024  # 1 MB

        with open(dest, "wb") as f:
            while True:
                chunk = resp.read(block_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)

                if total:
                    pct = downloaded / total * 100
                    mb_done = downloaded / 1e6
                    mb_total = total / 1e6
                    print(f"\r    {mb_done:.1f} / {mb_total:.1f} MB ({pct:.0f}%)",
                          end="", flush=True)
                else:
                    print(f"\r    {downloaded / 1e6:.1f} MB",
                          end="", flush=True)
        print()


def list_folder_contents(folder_href: str, token: Optional[str] = None) -> List[dict]:
    """Recursively list all files in an OSF folder."""
    files = []
    url = folder_href

    while url:
        data = api_get(url, token)
        for item in data.get("data", []):
            kind = item["attributes"]["kind"]
            name = item["attributes"]["name"]

            if kind == "file":
                files.append({
                    "name": name,
                    "path": item["attributes"].get("materialized_path", name),
                    "size": item["attributes"].get("size", 0),
                    "download": item["links"]["download"],
                })
            elif kind == "folder":
                sub_href = item["relationships"]["files"]["links"]["related"]["href"]
                files.extend(list_folder_contents(sub_href, token))

        url = data.get("links", {}).get("next")

    return files


def find_folder(folder_name: str, token: Optional[str] = None) -> Optional[str]:
    """Find a top-level folder in the OSF project and return its files API URL."""
    url = f"{API_BASE}/nodes/{OSF_PROJECT}/files/osfstorage/"
    data = api_get(url, token)

    for item in data.get("data", []):
        if (item["attributes"]["kind"] == "folder" and
                item["attributes"]["name"] == folder_name):
            return item["relationships"]["files"]["links"]["related"]["href"]
    return None


def list_top_folders(token: Optional[str] = None) -> List[str]:
    """List all top-level folders in the OSF project."""
    url = f"{API_BASE}/nodes/{OSF_PROJECT}/files/osfstorage/"
    data = api_get(url, token)
    folders = []
    for item in data.get("data", []):
        if item["attributes"]["kind"] == "folder":
            folders.append(item["attributes"]["name"])
    return folders


def download_folder(folder_name: str, out_dir: Path,
                    token: Optional[str] = None) -> int:
    """Download all files in an OSF folder to a local directory."""
    print(f"\nLooking for folder: {folder_name} ...")
    href = find_folder(folder_name, token)
    if href is None:
        print(f"ERROR: Folder '{folder_name}' not found in OSF project {OSF_PROJECT}")
        return 1

    print(f"Listing contents ...")
    files = list_folder_contents(href, token)
    if not files:
        print(f"  No files found in {folder_name}/")
        return 0

    total_size = sum(f["size"] or 0 for f in files)
    print(f"  Found {len(files)} files ({total_size / 1e9:.2f} GB)")
    print()

    downloaded = 0
    for i, f in enumerate(files, 1):
        rel_path = f["path"].lstrip("/")
        dest = out_dir / rel_path
        size_mb = (f["size"] or 0) / 1e6

        if dest.exists() and f["size"] and dest.stat().st_size == f["size"]:
            print(f"  [{i}/{len(files)}] SKIP (exists) {rel_path}")
            continue

        print(f"  [{i}/{len(files)}] Downloading {rel_path} ({size_mb:.1f} MB)")

        retries = 3
        for attempt in range(retries):
            try:
                download_file(f["download"], dest, token)
                downloaded += 1
                break
            except (urllib.error.URLError, OSError) as e:
                if attempt < retries - 1:
                    print(f"    Retry {attempt + 2}/{retries} ...")
                    time.sleep(2 ** attempt)
                else:
                    print(f"    FAILED: {e}")

    print(f"\nDone. Downloaded {downloaded} files to {out_dir}/")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Download WhatsGNU-ATB data from OSF (project: 6jr4u)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download the database (required for querying)
  python download_osf.py --folder WGNU_ATB_DB --out-dir ./whatsgnu_db

  # Download database + sample tables
  python download_osf.py --folder WGNU_ATB_DB Sample_tables --out-dir ./whatsgnu_db

  # Download everything
  python download_osf.py --all --out-dir ./whatsgnu_db

  # List available folders
  python download_osf.py --list
""")

    ap.add_argument("--folder", nargs="+", metavar="NAME",
                    help="Folder(s) to download: WGNU_ATB_DB, Sample_tables, "
                         "ATB_hash_seq, ATB_summary_figures_tables")
    ap.add_argument("--all", action="store_true",
                    help="Download all folders")
    ap.add_argument("--list", action="store_true",
                    help="List available folders and exit")
    ap.add_argument("--out-dir", type=str, default="./whatsgnu_atb_data",
                    help="Local directory to save files (default: ./whatsgnu_atb_data)")
    ap.add_argument("--token", type=str, default=None,
                    help="OSF personal access token (or set OSF_TOKEN env var). "
                         "Not required for public projects.")

    args = ap.parse_args()
    token = args.token or os.environ.get("OSF_TOKEN")

    if args.list:
        print(f"OSF project: https://osf.io/{OSF_PROJECT}/")
        print(f"\nAvailable folders:\n")
        folders = list_top_folders(token)
        for name in sorted(folders):
            desc = FOLDER_DESCRIPTIONS.get(name, "")
            print(f"  {name}")
            if desc:
                print(f"      {desc}")
            print()
        return 0

    if not args.folder and not args.all:
        ap.print_help()
        print("\nERROR: Specify --folder NAME(s) or --all")
        return 1

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.all:
        folders = list_top_folders(token)
    else:
        folders = args.folder

    print("=" * 50)
    print("  WhatsGNU-ATB OSF Downloader")
    print("=" * 50)
    print(f"  Project:  https://osf.io/{OSF_PROJECT}/")
    print(f"  Output:   {out_dir}")
    print(f"  Folders:  {', '.join(folders)}")

    rc = 0
    for folder_name in folders:
        result = download_folder(folder_name, out_dir, token)
        if result != 0:
            rc = result

    if rc == 0 and "WGNU_ATB_DB" in folders:
        print("\n" + "=" * 50)
        print("  Database ready! Query with:")
        print("=" * 50)

        db_path = out_dir / "WGNU_ATB_DB"
        st_path = out_dir / "Sample_tables" / "samples_with_ids.tsv"
        st_flag = ""
        if st_path.exists():
            st_flag = (f" \\\n      --samples_tsv {st_path}"
                       f" \\\n      --species_names_tsv {st_path}")

        print(f"""
  python Query_WhatsGNU_ATB.py \\
      --db_dir {db_path} \\
      --shards 8 \\
      --faa your_genome.bakta.faa \\
      --with_postings{st_flag} \\
      --out_dir results/
""")

    return rc


if __name__ == "__main__":
    raise SystemExit(main())
