#!/bin/bash
# Packs the results of one practical into a single zip file in ~/BMSi (reports,
# tables and figures; no raw data and no internal caches), so that everything
# needed to answer the questions can be kept on your own computer:
#
#   bash ~/BMSi/scripts/pack_results.sh PC1      (or PC2, PC3, PC4)
#
# Then, in RStudio: Files pane → folder BMSi → tick PC1_results.zip → More → Export.
# Practical 1 also includes the files for IGV (reference, contigs, BAM files and their
# indexes, coverage and contig tracks): about 300 MB.
B=~/BMSi
PC=$(echo "$1" | tr '[:lower:]' '[:upper:]')
case "$PC" in
  PC1) DIR=$B/PC1_assembly
       ITEMS="multiqc_report quast_report_megahit qc_raw qc_clean fastp.html fastp.json strategy_comparison_mqc.tsv
              contigs_vs_ref.flagstat ref_aligned_final.stats megahit_output_final/contigs_clean.fasta
              megahit_output_final/contigs_clean.fasta.fai *.fasta *.fna *.fai
              coverage_1kb.bedgraph contigs_on_ref.bed
              ref_aligned_final.bam ref_aligned_final.bam.bai contigs_vs_ref.bam contigs_vs_ref.bam.bai
              ref_vs_contigs.bam ref_vs_contigs.bam.bai reads_vs_contigs.bam reads_vs_contigs.bam.bai"
       KEEP_BAM=yes ;;   # the IGV files (sections 4.11 and 4.12) go in too
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
EXCL="'*.rds' '*.stan' '*.key' '*Bayesian_*_*.csv' '*.hdf5' '*.fastq*'"
[ "$KEEP_BAM" = yes ] || EXCL="$EXCL '*.bam' '*.bai'"
eval zip -q -r "$OUTZIP" $FOUND -x $EXCL
echo "Results of $PC packed in: $OUTZIP ($(du -h "$OUTZIP" | cut -f1))"
echo "Download it: RStudio → Files → BMSi → tick ${PC}_results.zip → More → Export."
