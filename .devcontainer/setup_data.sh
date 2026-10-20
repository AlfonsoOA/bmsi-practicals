#!/bin/bash
# Runs once, when a student's codespace is created: unpacks the practical data
# (data/ForStudents_*.zip) and copies the scripts into ~/BMSi, the folder the
# scripts use by default. Files already there are never overwritten, so running
# it again is harmless.
set -e
REPO=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p ~/BMSi
for z in "$REPO"/data/ForStudents_*.zip; do
  [ -e "$z" ] || continue
  echo "Unpacking $(basename "$z")"
  unzip -n -q "$z" -d ~/BMSi
done
mkdir -p ~/BMSi/scripts
cp -n "$REPO"/scripts/* ~/BMSi/scripts/ 2>/dev/null || true
# RStudio starts in ~/BMSi (Files pane and working directory)
mkdir -p ~/.config/rstudio
[ -f ~/.config/rstudio/rstudio-prefs.json ] || \
  echo '{ "initial_working_directory": "~/BMSi" }' > ~/.config/rstudio/rstudio-prefs.json
echo "Data and scripts are in ~/BMSi:"; ls ~/BMSi
