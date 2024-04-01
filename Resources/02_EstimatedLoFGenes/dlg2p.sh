#!/bin/bash

readonly BASE_URI="https://www.ebi.ac.uk/gene2phenotype/downloads"
readonly URIS_FILE="./uris.txt"

# Remove old uris.txt file if it exists
function remove_old_uris_file() {
  if [ -e "$URIS_FILE" ]; then
    echo "Old URI file exists."
    rm -rf "$URIS_FILE"
    echo "Removed old uris.txt"
  fi
}

# Download function using curl
function download_with_curl() {
  echo "Downloading with curl: $1"
  curl -O -s "$1"
}

# Download function using wget
function download_with_wget() {
  echo "Downloading with wget: $1"
  wget -nv "$1"
}

# Determine which download function to use
function download_file() {
  if command -v curl > /dev/null; then
    download_with_curl "$1"
  elif command -v wget > /dev/null; then
    download_with_wget "$1"
  else
    echo "Neither curl nor wget is installed."
    exit 1
  fi
}

# Download all CSV files
function download_csvs() {
  remove_old_uris_file
  touch "$URIS_FILE"
  echo "${BASE_URI}/DDG2P.csv.gz" >> "$URIS_FILE"
  echo "${BASE_URI}/EyeG2P.csv.gz" >> "$URIS_FILE"
  echo "${BASE_URI}/SkinG2P.csv.gz" >> "$URIS_FILE"
  echo "${BASE_URI}/CancerG2P.csv.gz" >> "$URIS_FILE"
  echo "${BASE_URI}/CardiacG2P.csv.gz" >> "$URIS_FILE"
  echo "${BASE_URI}/SkeletalG2P.csv.gz" >> "$URIS_FILE"

  while IFS= read -r uri; do
    download_file "$uri"
  done < "$URIS_FILE"
  
  rm "$URIS_FILE"
}

function main() {
  download_csvs
}

main
