#!/usr/bin/env bash
# Upload locally built genes.db and FASTA to an S3 bucket.
# Uses SEQMAT_DATA_DIR as data root (or first arg); bucket from second arg or SEQMAT_PREBUILT_BUCKET.
# Requires: aws cli, gzip (for FASTA if not already .gz)

set -euo pipefail

DATA_ROOT="${1:-${SEQMAT_DATA_DIR:-}}"
BUCKET="${2:-${SEQMAT_PREBUILT_BUCKET:-}}"

if [[ -z "$DATA_ROOT" ]]; then
  echo "Usage: $0 [DATA_ROOT] [BUCKET]" >&2
  echo "  DATA_ROOT: path to seqmat data root (e.g. /path/to/seqmat_base), or set SEQMAT_DATA_DIR" >&2
  echo "  BUCKET: S3 bucket name (e.g. seqmat-prebuilt-public), or set SEQMAT_PREBUILT_BUCKET" >&2
  exit 1
fi

if [[ -z "$BUCKET" ]]; then
  echo "BUCKET required (second argument or SEQMAT_PREBUILT_BUCKET)" >&2
  exit 1
fi

if ! command -v aws &>/dev/null; then
  echo "aws cli not found. Install it and configure credentials." >&2
  exit 1
fi

for org in hg38 mm39; do
  dir="$DATA_ROOT/$org"
  if [[ ! -d "$dir" ]]; then
    echo "Skip $org (no dir $dir)"
    continue
  fi

  db="$dir/genes.db"
  if [[ -f "$db" ]]; then
    echo "Uploading $org/genes.db ..."
    aws s3 cp "$db" "s3://$BUCKET/$org/genes.db"
  fi

  fa="$dir/$org.fa"
  fa_gz="$dir/$org.fa.gz"
  if [[ -f "$fa" ]]; then
    echo "Uploading $org/$org.fa.gz (from $org.fa) ..."
    gzip -c "$fa" | aws s3 cp - "s3://$BUCKET/$org/$org.fa.gz"
  elif [[ -f "$fa_gz" ]]; then
    echo "Uploading $org/$org.fa.gz ..."
    aws s3 cp "$fa_gz" "s3://$BUCKET/$org/$org.fa.gz"
  fi

  cons_db="$dir/conservation.db"
  cons_pkl="$dir/conservation.pkl"
  if [[ -f "$cons_db" ]]; then
    echo "Uploading $org/conservation.db ..."
    aws s3 cp "$cons_db" "s3://$BUCKET/$org/conservation.db"
  fi
  if [[ -f "$cons_pkl" ]]; then
    echo "Uploading $org/conservation.pkl ..."
    aws s3 cp "$cons_pkl" "s3://$BUCKET/$org/conservation.pkl"
  fi

  echo "Done $org."
done

echo "Upload complete."
