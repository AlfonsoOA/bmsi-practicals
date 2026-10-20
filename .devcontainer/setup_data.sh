#!/bin/bash
# Runs once, when a student's codespace is created: puts the practical data and
# scripts into ~/BMSi, the folders the scripts use by default. Files already
# there are never overwritten, so running it again is harmless.
#   - data/ForStudents_PC3.zip, data/ForStudents_PC4.zip (in the repository)
#   - ForStudents_PC1.zip, ForStudents_PC2.zip (too large for the repository:
#     downloaded from the release "data-2026-27" of this repository)
set -e
REPO=$(cd "$(dirname "$0")/.." && pwd)
RELEASE=https://github.com/AlfonsoOA/bmsi-practicals/releases/download/data-2026-27
mkdir -p ~/BMSi

for z in "$REPO"/data/ForStudents_*.zip; do
  [ -e "$z" ] || continue
  echo "Unpacking $(basename "$z")"
  unzip -n -q "$z" -d ~/BMSi
done

# Large packages: a zip with a ForStudents/ folder, unpacked inside the practical's folder
get_release() {   # zip name, practical folder
  local dest=~/BMSi/$2
  if [ -d "$dest/ForStudents" ]; then echo "$2/ForStudents already there"; return; fi
  mkdir -p "$dest"
  echo "Downloading $1 (this can take a few minutes)..."
  if curl -fsSL --retry 3 -o "$dest/$1" "$RELEASE/$1"; then
    unzip -n -q "$dest/$1" -d "$dest" && rm -f "$dest/$1"
  else
    rm -f "$dest/$1"; echo "WARNING: $1 could not be downloaded; run this script again later: bash $0"
  fi
}
get_release ForStudents_PC1.zip PC1_assembly
get_release ForStudents_PC2.zip PC2_RNAseq

mkdir -p ~/BMSi/scripts
cp -n "$REPO"/scripts/* ~/BMSi/scripts/ 2>/dev/null || true

# RStudio starts in ~/BMSi (Files pane and working directory)
mkdir -p ~/.config/rstudio
[ -f ~/.config/rstudio/rstudio-prefs.json ] || \
  echo '{ "initial_working_directory": "~/BMSi" }' > ~/.config/rstudio/rstudio-prefs.json
echo "Data and scripts are in ~/BMSi:"; ls ~/BMSi
