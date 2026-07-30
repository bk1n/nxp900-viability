#!/usr/bin/env bash

set -euo pipefail
source .env

cp -r "$DATA_PATH" .
echo "Data fetched from $DATA_PATH on $(date)"