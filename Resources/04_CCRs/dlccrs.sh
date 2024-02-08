#!/bin/bash

cd $(dirname "$0")

readonly BASE_URI="https://s3.us-east-2.amazonaws.com/ccrs/ccrs"
readonly autosomes="ccrs.autosomes.v2.20180420.bed.gz"
readonly xchrom="ccrs.xchrom.v2.20180420.bed.gz"


# Remove old uris.txt file if it exists
function remove_old_files() {
  if [ -e "$autosomes" ]; then
    echo "Old autosomes file exists."
    rm -rf "$autosomesE"
    echo "Removed old files"
  fi
  if [ -e "$xchrom" ]; then
    echo "Old xchrom file exists."
    rm -rf "$xchrom"
    echo "Removed old files"
  fi
}

# Download function using curl
function download_with_curl() {
  echo "Downloading with curl: $1"
  curl -O -sL "$1"
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

# Download all BED files
function download_ccrs() {
  remove_old_files
  download_file "${BASE_URI}/${autosomes}"
  download_file "${BASE_URI}/${xchrom}"
}

function main() {
  download_ccrs
}

main
