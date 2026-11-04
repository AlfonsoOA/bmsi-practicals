#!/bin/bash
# Packs the results of one practical into a single zip file in ~/BMSi (reports,
# tables and figures; no raw data and no internal caches), so that everything
# needed to answer the questions can be kept on your own computer:
#
#   bash ~/BMSi/scripts/pack_results.sh PC1      (or PC2, PC3, PC4)
#
# Then, in RStudio: Files pane → folder BMSi → tick PC1_results.zip → More → Export.
# The IGV files of Practical 1 (BAM files) are not included: download them as the
# protocol explains (section 4.11).
B=~/BMSi
PC=$(echo "$1" | tr '[:lower:]' '[:upper:]')
case "$PC" in
  PC1) DIR=$B/PC1_assembly
       ITEMS="multiqc_report quast_report_megahit qc_raw qc_clean fastp.html fastp.json strategy_comparison_mqc.tsv
              contigs_vs_ref.flagstat ref_aligned_final.stats megahit_output_final/contigs_clean.fasta *.fasta" ;;
  PC2) DIR=$B
       ITEMS="PC2_RNAseq/multiqc_report PC2_RNAseq/busco_*/short_summary*.txt PC2_RNAseq_R/results*" ;;
  PC3) DIR=$B/PC3_Proteomics; ITEMS="results" ;;
  PC4) DIR=$B/PC4_Enrichment; ITEMS="results" ;;
  *)   echo "Write the practical: bash $0 PC1   (PC1, PC2, PC3 or PC4)"; exit 1 ;;
esac
cd "$DIR" 2>/dev/null || { echo "Folder $DIR not found: run the practical first."; exit 1; }
FOUND=""
for i in $ITEMS; do [ -e "$i" ] && FOUND="$FOUND $i"; done
if [ -z "$FOUND" ]; then echo "No results of $PC found in $DIR: run the practical first."; exit 1; fi
OUTZIP=$B/${PC}_results.zip
rm -f "$OUTZIP"
# Internal caches and trained models are left out: they are not needed to read the results
zip -q -r "$OUTZIP" $FOUND -x '*.rds' '*.stan' '*.key' '*Bayesian_*_*.csv' '*.hdf5' '*.bam' '*.bai' '*.fastq*'
echo "Results of $PC packed in: $OUTZIP ($(du -h "$OUTZIP" | cut -f1))"
echo "Download it: RStudio → Files → BMSi → tick ${PC}_results.zip → More → Export."
