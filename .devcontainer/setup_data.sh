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
# The data go to the large temporary disk of the codespace (/tmp, >100 GB free) through
# a link: the image almost fills the 32 GB main disk, and downloading and unpacking the
# data there can run out of space ("No space left on device").
if [ ! -e ~/BMSi ] && [ "$(df --output=avail -k /tmp 2>/dev/null | tail -1)" -gt 20000000 ] 2>/dev/null; then
  mkdir -p /tmp/BMSi && ln -s /tmp/BMSi ~/BMSi
fi
mkdir -p ~/BMSi

# First of all (RStudio can open before the downloads finish): RStudio starts in ~/BMSi (Files pane and working directory)
mkdir -p ~/.config/rstudio
[ -f ~/.config/rstudio/rstudio-prefs.json ] || \
  echo '{ "initial_working_directory": "~/BMSi" }' > ~/.config/rstudio/rstudio-prefs.json

# Terminal: pasting multi-line blocks must work as in a local terminal. The RStudio
# terminal in the browser marks pasted text with escape codes ("bracketed paste")
# that break the first line of if/for blocks, so it is switched off.
grep -q "enable-bracketed-paste" ~/.inputrc 2>/dev/null || echo "set enable-bracketed-paste off" >> ~/.inputrc

# conda ready in every new terminal, so that "conda activate bmsi_pc1" (or bmsi_pc2) works
# directly, as the protocols say
grep -q "conda shell.bash hook" ~/.bashrc 2>/dev/null || \
  echo 'command -v conda > /dev/null 2>&1 && eval "$(conda shell.bash hook)"' >> ~/.bashrc

# The preference above is not always honoured (RStudio can still open in
# /workspaces/bmsi-practicals), so interactive R sessions also move to ~/BMSi
# on start. Rscript (the automatic tests) is not interactive and is unaffected.
grep -q "BMSi" ~/.Rprofile 2>/dev/null || \
  echo 'if (interactive() && dir.exists("~/BMSi")) setwd("~/BMSi")' >> ~/.Rprofile

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
# Disk: the codespace has 32 GB and the image takes most of them (~4 GB free after the
# data). The trimmed reads provided in ForStudents of PC1 (1.4 GB) are not used: every
# student creates their own with fastp (section 4.2).
rm -f ~/BMSi/PC1_assembly/ForStudents/3_*_paired.fastq

mkdir -p ~/BMSi/scripts
cp -n "$REPO"/scripts/* ~/BMSi/scripts/ 2>/dev/null || true

# In a codespace, RStudio starts now (PC1, PC3 and PC4 are ready) and the large package of
# PC2 (2.1 GB) is downloaded afterwards in the background, so nobody waits for it.
# Elsewhere (automatic tests) everything is done in order.
if [ "$CODESPACES" = "true" ]; then
  bash "$REPO/.devcontainer/start_rstudio.sh" || true
  setsid nohup bash -c "$(declare -f get_release); RELEASE=$RELEASE; get_release ForStudents_PC2.zip PC2_RNAseq; echo PC2 data ready" \
    > ~/BMSi/PC2_download.log 2>&1 < /dev/null &
  echo "PC2 data are being downloaded in the background (see ~/BMSi/PC2_download.log)."
else
  get_release ForStudents_PC2.zip PC2_RNAseq
fi
echo "Data and scripts are in ~/BMSi:"; ls ~/BMSi
