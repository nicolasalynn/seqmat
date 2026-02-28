#!/usr/bin/env bash
# Copy mm39 conservation.pkl from old SeqMat bucket to new bucket.
# Old: genome-data-public-access (eu-north-1), New: seqmat-prebuilt-public (us-east-1).
# Run once; requires read access to old bucket and write access to new bucket.

set -euo pipefail

OLD_BUCKET="genome-data-public-access"
NEW_BUCKET="${1:-seqmat-prebuilt-public}"
KEY="mm39/conservation.pkl"

echo "Copying s3://$OLD_BUCKET/$KEY -> s3://$NEW_BUCKET/$KEY"
aws s3 cp "s3://$OLD_BUCKET/$KEY" "s3://$NEW_BUCKET/$KEY" \
  --source-region eu-north-1 \
  --region us-east-1

echo "Done."
