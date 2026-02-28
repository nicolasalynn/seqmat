#!/usr/bin/env python3
"""
Convert conservation.pkl to conservation.db (SQLite, same pattern as genes.db).
Usage:
  python scripts/convert_conservation_to_db.py [PATH]
  If PATH is a dir, looks for conservation.pkl in it and writes conservation.db there.
  If PATH is a file, treats it as conservation.pkl and writes same dir with .db.
"""
import sys
from pathlib import Path

# Run from repo root so seqmat is importable
repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repo_root))

from seqmat.utils import convert_conservation_pkl_to_db


def main():
    if len(sys.argv) < 2:
        print("Usage: convert_conservation_to_db.py <path to dir or conservation.pkl>", file=sys.stderr)
        sys.exit(1)
    path = Path(sys.argv[1]).resolve()
    if path.is_file():
        if path.suffix.lower() != ".pkl":
            print("Expected a .pkl file", file=sys.stderr)
            sys.exit(1)
        out = convert_conservation_pkl_to_db(path)
    else:
        pkl_path = path / "conservation.pkl"
        if not pkl_path.exists():
            print(f"Not found: {pkl_path}", file=sys.stderr)
            sys.exit(1)
        out = convert_conservation_pkl_to_db(pkl_path, path / "conservation.db")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
