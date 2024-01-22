#!/bin/bash

set -eu

readonly BASE_URI="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
readonly VERSION="release_43"
readonly FILE_NAME="GRCh37_mapping/gencode.v43lift37.annotation.gff3.gz"
readonly GENCODE_URI="${BASE_URI}/${VERSION}/${FILE_NAME}"

wget "${GENCODE_URI}"

gunzip -c gencode.v43lift37.annotation.gff3.gz \
  >gencode.v43lift37.annotation.gff3

mv gencode.v43lift37.annotation.gff3.gz \
  gencode.v43lift37.annotation.gff3.original.gz

bedtools sort -i gencode.v43lift37.annotation.gff3 \
  >gencode.v43lift37.annotation.sort.gff3

rm -rf gencode.v43lift37.annotation.gff3

bgzip -@ 4 gencode.v43lift37.annotation.sort.gff3

tabix gencode.v43lift37.annotation.sort.gff3.gz