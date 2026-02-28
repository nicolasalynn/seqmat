# Scripts

Helper scripts for maintainers: uploading prebuilt data to S3 and converting conservation format.

## upload_prebuilt_to_s3.sh

Upload `genes.db`, FASTA, and conservation data to the SeqMat prebuilt S3 bucket (e.g. `seqmat-prebuilt-public`).

**Usage:**
```bash
./scripts/upload_prebuilt_to_s3.sh [DATA_ROOT] [BUCKET]
```
- `DATA_ROOT`: path to seqmat data root (where `hg38/`, `mm39/` live), or set `SEQMAT_DATA_DIR`
- `BUCKET`: S3 bucket name in the account you’re authenticated to

Uploads per organism: `genes.db` (if present), `{org}.fa.gz` (from `{org}.fa` or `{org}.fa.gz`), and when present both `conservation.db` and `conservation.pkl`.

**Example:**
```bash
export SEQMAT_PREBUILT_BUCKET=seqmat-prebuilt-public
./scripts/upload_prebuilt_to_s3.sh
# or
./scripts/upload_prebuilt_to_s3.sh /path/to/seqmat_base seqmat-prebuilt-public
```

See repo README (Data Setup) for S3 bucket setup: public read, restricted upload.

---

## copy_mm39_conservation_from_old_bucket.sh

One-time copy of mm39 conservation from the old bucket to the new one.

```bash
./scripts/copy_mm39_conservation_from_old_bucket.sh [NEW_BUCKET]
# Default NEW_BUCKET is seqmat-prebuilt-public. Requires read on old bucket, write on new.
```

If the object is at a different key on the old bucket, edit `KEY` in the script.

---

## convert_conservation_to_db.py

Convert `conservation.pkl` to `conservation.db` (SQLite, same pattern as genes.db).

**Usage:**
```bash
# From repo root; path = directory containing conservation.pkl, or path to the .pkl file
python scripts/convert_conservation_to_db.py /path/to/hg38
# or
python scripts/convert_conservation_to_db.py /path/to/hg38/conservation.pkl
```
Writes `conservation.db` in the same directory. Prebuilt download tries `.db` first, then `.pkl`.
