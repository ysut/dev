#!/bin/bash

SCRIPT_DIR=$(cd $(dirname $0); pwd)

AM_HG19_URL="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz"
AM_HG38_URL="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
HG19_FILE="AlphaMissense_hg19.tsv.gz"
HG38_FILE="AlphaMissense_hg38.tsv.gz"

set -eu

echo "Downloading files......"
mkdir -P AlphaMissense; cd $_
curl -sSL ${AM_HG19_URL} -o ${HG19_FILE}
curl -sSL ${AM_HG38_URL} -o ${HG38_FILE}

cd $SCRIPT_DIR