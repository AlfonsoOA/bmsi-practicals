#!/bin/bash
# Test run of the practicals inside the image, as a student would run them with
# the ForStudents data (heavy steps precomputed). Prints a markdown table with
# exit status, time and peak memory of each script. Logs go to test_logs/.
REPO=$(pwd)
LOGS=$REPO/test_logs; mkdir -p $LOGS
bash .devcontainer/setup_data.sh > $LOGS/setup.log 2>&1

# PC2 (R part): its input is produced by the terminal part; for the test, the
# instructor's quantification files are used
if [ ! -d ~/BMSi/PC2_RNAseq_R ] && [ -f tests/PC2_RNAseq_R.tar.gz ]; then
  mkdir -p ~/BMSi/PC2_RNAseq_R && tar xzf tests/PC2_RNAseq_R.tar.gz -C ~/BMSi/PC2_RNAseq_R
fi

# PC3: second pass (reinjections excluded), the configuration the caches belong to
sed 's|^EXCLUDE_RUNS <- c()|EXCLUDE_RUNS <- c("SOL2023001F_1_S2-A2_1_12740", "SOL2023001F_5_S2-A7_1_12744", "SOL2023001F_7_S2-A10_1_12737")|' \
  ~/BMSi/scripts/PC3_Proteomics.R > /tmp/PC3_Proteomics_pass2.R

run() {   # name, script
  /usr/bin/time -v -o $LOGS/$1.time Rscript "$2" > $LOGS/$1.log 2>&1
  local st=$?
  local mem=$(awk -F': ' '/Maximum resident/ {printf "%.2f", $2/1048576}' $LOGS/$1.time)
  local t=$(awk -F': ' '/Elapsed/ {print $2}' $LOGS/$1.time)
  echo "| $1 | $([ $st -eq 0 ] && echo OK || echo "FAILED ($st)") | $t | $mem GB |"
}
echo "| Practical | Result | Time (h:mm:ss) | Peak memory |"
echo "| --- | --- | --- | --- |"
run PC2_RNAseq      ~/BMSi/scripts/PC2_RNAseq.R
run PC3_Proteomics  /tmp/PC3_Proteomics_pass2.R
run PC4_Enrichment  ~/BMSi/scripts/PC4_Enrichment.R
run PC4b_MOFA       ~/BMSi/scripts/PC4b_MOFA.R

# Terminal tools of PC1 and PC2
source /opt/conda/etc/profile.d/conda.sh
for e in bmsi_pc1 bmsi_pc2; do
  conda env list | grep -q "^$e " && r=OK || r=MISSING
  echo "| conda env $e | $r | | |"
done
for t in megahit fastp salmon quast.py prefetch blastn; do
  e=bmsi_pc1; [ $t = salmon ] && e=bmsi_pc2
  conda run -n $e which $t > /dev/null 2>&1 && r=OK || r=MISSING
  echo "| $t ($e) | $r | | |"
done
command -v rstudio-server > /dev/null && echo "| RStudio Server | OK | | |" || echo "| RStudio Server | MISSING | | |"
