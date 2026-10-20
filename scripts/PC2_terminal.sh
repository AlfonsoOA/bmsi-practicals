#!/usr/bin/env bash
# =============================================================================
# PC2_terminal.sh - Systems Molecular Biology, Practical 2 (RNA-seq), 2026-2027
# Complete terminal workflow: installation, quality control, quantification
# against the Ensembl reference and against a de novo Trinity assembly, BUSCO,
# Trinity -> Ensembl gene mapping, MultiQC and export of the files for R.
#
# Usage (Windows/WSL, macOS or Linux):
#   bash PC2_terminal.sh 2>&1 | tee PC2_terminal.log
#
# Long steps (download, Trinity, BUSCO) are skipped when their outputs exist.
# To analyse other data, change only the lines between the dashed comments
# (WORKDIR, ENS_URL, UNIPROT_PROTEOME, BUSCO_LINEAGE, MEM) and samples.tsv.
# Same commands, in the same order, as Appendix A of the guide.
# =============================================================================

# --- Install the tools (once) ---
MAIN_PKGS="python=3.11 fastqc=0.12.1 fastp=1.3.7 multiqc=1.35 salmon=2.8.0 minimap2=2.31 seqkit=2.14.0 sra-tools=3.4.1"
TRINITY_PKGS="trinity=2.15.2"
BUSCO_PKGS="busco=6.1.0"
ENVS_DIR=~/BMSi/envs   # exact environment files from the course materials, if available

# 1. Install conda (Miniforge) only if there is no conda at all
if ! command -v conda > /dev/null 2>&1 && [ ! -d "$HOME/miniforge3" ]; then
  curl -L -o Miniforge3.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  bash Miniforge3.sh -b -p "$HOME/miniforge3"
  rm Miniforge3.sh
  "$HOME/miniforge3/bin/conda" init bash
fi

# 2. Make conda usable in this terminal
if command -v conda > /dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
else
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
fi

# 3. Macs with Apple chips use the Intel versions of the tools (through Rosetta)
if [ "$(uname)-$(uname -m)" = "Darwin-arm64" ]; then
  export CONDA_SUBDIR=osx-64
fi

# 4. Create each environment only if it does not exist yet
make_env() {
  if conda env list | grep -q "^$1 "; then
    echo "Environment $1 already exists: nothing to install."
  elif [ -f "$ENVS_DIR/$1_linux-64.txt" ] && [ "$(uname)-$(uname -m)" = "Linux-x86_64" ]; then
    conda create -y -n "$1" --file "$ENVS_DIR/$1_linux-64.txt"   # exact copy of the instructor environment
  else
    conda create -y -n "$1" --override-channels -c conda-forge -c bioconda $2
  fi
}
make_env bmsi_pc2 "$MAIN_PKGS"
make_env bmsi_pc2_trinity "$TRINITY_PKGS"
make_env bmsi_pc2_busco "$BUSCO_PKGS"
conda activate bmsi_pc2

# --- Session: every new terminal ---
if ! command -v salmon > /dev/null 2>&1; then
  if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
  else
    eval "$(conda shell.bash hook)"
  fi
  conda activate bmsi_pc2
fi

# ---- Folder and references: the only lines to change for other data ----
WORKDIR=~/BMSi/PC2_RNAseq
ENS_URL=https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/chlamydomonas_reinhardtii/cdna/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.cdna.all.fa.gz
UNIPROT_PROTEOME=UP000006906
BUSCO_LINEAGE=chlorophyceae_odb12.2
MEM=24G
# (the samples are listed in the file samples.tsv, section 4.2)
# -------------------------------------------------------------------------
THREADS=$(getconf _NPROCESSORS_ONLN)
mkdir -p $WORKDIR && cd $WORKDIR

# Folder that R will read: in Windows (WSL), inside your Windows user folder
if grep -qi microsoft /proc/version 2> /dev/null; then
  WINUSER=$(cmd.exe /c "echo %USERNAME%" 2> /dev/null | tr -d '\r')
  if [ -z "$WINUSER" ]; then   # Windows commands unavailable (e.g. inside tmux): use the user folder
    WINUSER=$(ls /mnt/c/Users | grep -vxE 'Public|Default|Default User|All Users|defaultuser0|WsiAccount|desktop.ini' | head -1)
  fi
  R_DIR=/mnt/c/Users/$WINUSER/BMSi/PC2_RNAseq_R
else
  R_DIR=~/BMSi/PC2_RNAseq_R
fi
echo "Folder: $WORKDIR | Threads: $THREADS | Trinity memory: $MEM | R folder: $R_DIR"

# --- Samples ---
if [ ! -f samples.tsv ]; then
  printf 'sample\tcondition\trun\n' > samples.tsv
  printf 'C1\tControl\tSRR28217878\n' >> samples.tsv
  printf 'C2\tControl\tSRR28217875\n' >> samples.tsv
  printf 'C3\tControl\tSRR28217869\n' >> samples.tsv
  printf 'NaCl1\tNaCl\tSRR28217842\n' >> samples.tsv
  printf 'NaCl2\tNaCl\tSRR28217877\n' >> samples.tsv
  printf 'NaCl3\tNaCl\tSRR28217881\n' >> samples.tsv
fi
cat samples.tsv

# --- Provided files (students) ---
if [ -d ForStudents ]; then
  mv -n ForStudents/* .
fi

# --- [LONG] Raw reads ---
mkdir -p raw
for S in $(tail -n +2 samples.tsv | cut -f1); do
  RUN=$(awk -F'\t' -v s="$S" '$1==s {print $3}' samples.tsv)
  if [ -s raw/$S.fastq.gz ]; then
    echo "$S: raw reads found, skipping download."
  else
    if time prefetch -p $RUN && time fasterq-dump -p -e $THREADS $RUN -o raw/$S.fastq; then
      gzip raw/$S.fastq
      rm -rf $RUN
    else
      echo "ERROR: $S ($RUN) could not be downloaded. Check its accession in samples.tsv."
    fi
  fi
done
ls -lh raw

# --- Quality control and cleaning ---
zcat raw/$(tail -n +2 samples.tsv | head -1 | cut -f1).fastq.gz | head -4
mkdir -p qc_raw
fastqc -t $THREADS -o qc_raw raw/*.fastq.gz
mkdir -p clean
for S in $(tail -n +2 samples.tsv | cut -f1); do
  fastp -i raw/$S.fastq.gz -o clean/$S.fastq.gz \
        --cut_front --cut_tail --cut_mean_quality 20 \
        --qualified_quality_phred 20 --unqualified_percent_limit 40 \
        --length_required 36 --thread $THREADS \
        --html clean/$S.fastp.html --json clean/$S.fastp.json
done
mkdir -p qc_clean
fastqc -t $THREADS -o qc_clean clean/*.fastq.gz
seqkit stats -j $THREADS raw/*.fastq.gz clean/*.fastq.gz

# --- Reference transcriptome ---
mkdir -p reference
if [ ! -s reference/ensembl_cdna.fa ]; then
  wget -O reference/ensembl_cdna.fa.gz $ENS_URL
  gunzip reference/ensembl_cdna.fa.gz
fi
# 1. Transcript -> gene, and the Phytozome name of each gene (gene_symbol)
grep ">" reference/ensembl_cdna.fa | awk -v OFS='\t' '{
  tx = substr($1, 2); gene = ""; sym = ""
  for (i = 2; i <= NF; i++) {
    if ($i ~ /^gene:/) gene = substr($i, 6)
    if ($i ~ /^gene_symbol:/) sym = substr($i, 13)
  }
  print tx, gene, sym }' > reference/ensembl_annotation.tsv
cut -f1,2 reference/ensembl_annotation.tsv > reference/ensembl_tx2gene.tsv

# 2. Protein names and GO terms of the reference proteome from UniProt (downloaded once)
if [ -z "$UNIPROT_PROTEOME" ]; then
  echo "ERROR: UNIPROT_PROTEOME is not defined. Run the session block (4.1) first."
elif [ ! -s reference/uniprot_proteome.tsv ]; then
  curl -L --retry 3 --speed-limit 1000 --speed-time 120 -o reference/uniprot_proteome.tsv.gz \
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&fields=accession,gene_orf,protein_name,go_id&query=proteome:$UNIPROT_PROTEOME"
  if gzip -t reference/uniprot_proteome.tsv.gz 2> /dev/null; then
    gunzip reference/uniprot_proteome.tsv.gz
  else
    echo "ERROR: the UniProt download failed or was too slow. Is UNIPROT_PROTEOME set (block 4.1)?"
    echo "If it is, the server is busy: run this block again later, or use the file provided by the instructor."
    head -c 300 reference/uniprot_proteome.tsv.gz | tr -cd '[:print:]\n'; echo
    rm -f reference/uniprot_proteome.tsv.gz
  fi
fi

# 3. One line per gene: Ensembl gene, Phytozome name, protein name, GO terms
if [ -s reference/uniprot_proteome.tsv ]; then
awk -F'\t' -v OFS='\t' '
  FNR == 1 { f++; if (f == 2) next }
  f == 1 { if (!($2 in sym)) sym[$2] = $3; next }
  { n = split($2, orf, /[ ;]+/)
    for (k = 1; k <= n; k++) { g = orf[k]
      if (!(g in sym)) continue
      if (!(g in name) || name[g] ~ /^Uncharacterized/) name[g] = $3
      if ($4 != "") { if (g in go) go[g] = go[g] "; " $4; else go[g] = $4 } } }
  END { print "gene", "symbol", "description", "GO"
        for (g in sym) print g, sym[g], name[g], go[g] }' \
  reference/ensembl_annotation.tsv reference/uniprot_proteome.tsv > reference/gene_annotation.tsv

awk -F'\t' 'NR > 1 { n++; if ($3 != "") u++; if ($3 != "" && $3 !~ /^Uncharacterized/) d++; if ($4 != "") g++ }
  END { print "Genes:", n, "| in UniProt:", u, "| with an informative name:", d, "| with GO terms:", g }' reference/gene_annotation.tsv
else
  echo "reference/uniprot_proteome.tsv not found: download it (step 2) before building the annotation."
fi

# --- Quantification against the reference ---
if [ ! -s salmon_index_ensembl/versionInfo.json ]; then
  time salmon index -t reference/ensembl_cdna.fa -i salmon_index_ensembl -p $THREADS
fi
mkdir -p quant_ensembl
for S in $(tail -n +2 samples.tsv | cut -f1); do
  time salmon quant -i salmon_index_ensembl -l A -r clean/$S.fastq.gz \
       -p $THREADS -o quant_ensembl/$S
done
for S in $(tail -n +2 samples.tsv | cut -f1); do
  echo "$S: library type $(grep -o '"expected_format": *"[A-Z]*"' quant_ensembl/$S/lib_format_counts.json | cut -d'"' -f4)"
done

# --- [LONG] De novo transcriptome ---
TRINITY_FA=trinity_assembly.Trinity.fasta
FIRST=$(tail -n +2 samples.tsv | head -1 | cut -f1)
LIBTYPE=$(grep -o '"expected_format": *"[A-Z]*"' quant_ensembl/$FIRST/lib_format_counts.json | cut -d'"' -f4)
case $LIBTYPE in
  SR) STRAND="--SS_lib_type R" ;;
  SF) STRAND="--SS_lib_type F" ;;
  *)  STRAND="" ;;
esac
echo "Library type detected by Salmon: $LIBTYPE -> Trinity strand option: ${STRAND:-none (unstranded)}"
if [ -s $TRINITY_FA ]; then
  echo "Trinity assembly found: skipping Trinity."
else
  READS=$(ls clean/*.fastq.gz | paste -sd, -)
  time conda run --no-capture-output -n bmsi_pc2_trinity \
    Trinity --seqType fq --single $READS --CPU $THREADS --max_memory $MEM \
            $STRAND --output trinity_assembly --full_cleanup --no_salmon
fi
seqkit stats -a $TRINITY_FA

# --- [LONG] Completeness ---
for A in ensembl trinity; do
  if [ $A = ensembl ]; then FA=reference/ensembl_cdna.fa; else FA=$TRINITY_FA; fi
  if ls busco_$A/short_summary*.txt > /dev/null 2>&1; then
    echo "BUSCO result for $A found: skipping."
  else
    time conda run --no-capture-output -n bmsi_pc2_busco \
      busco -i $FA -o busco_$A -m transcriptome -l $BUSCO_LINEAGE -c $THREADS -f
  fi
  echo "$A: $(grep -o 'C:.*' busco_$A/short_summary*.txt)"
done
E=$(ls busco_ensembl/run_*/full_table.tsv busco_ensembl/full_table.tsv 2> /dev/null | head -1)
T=$(ls busco_trinity/run_*/full_table.tsv busco_trinity/full_table.tsv 2> /dev/null | head -1)
: > busco_lost_in_trinity.tsv
printf 'Ensembl -> Trinity\tBUSCO genes\n'
awk -F'\t' -v OFS='\t' '
  FNR == 1 { f++ }
  /^# Busco id/ { if (f == 1) for (i = 1; i <= NF; i++) if ($i == "Description") dcol = i; next }
  /^#/ { next }
  f == 1 { if (!($1 in ens)) ens[$1] = $2; if (dcol && $dcol != "" && !($1 in desc)) desc[$1] = $dcol; next }
  { if (!($1 in tri)) tri[$1] = $2 }
  END { for (b in ens) { n[ens[b] " -> " tri[b]]++
          if (ens[b] ~ /Complete|Duplicated/ && tri[b] == "Missing") print b, desc[b] > "busco_lost_in_trinity.tsv" }
        for (k in n) print k, n[k] }' $E $T | sort -t$'\t' -k2,2nr
echo "Complete in the reference but missing in Trinity: $(wc -l < busco_lost_in_trinity.tsv) genes, e.g.:"
head -10 busco_lost_in_trinity.tsv

# --- Quantification against the de novo transcriptome ---
if [ ! -s salmon_index_trinity/versionInfo.json ]; then
  time salmon index -t $TRINITY_FA -i salmon_index_trinity -p $THREADS
fi
mkdir -p quant_trinity
for S in $(tail -n +2 samples.tsv | cut -f1); do
  time salmon quant -i salmon_index_trinity -l A -r clean/$S.fastq.gz \
       -p $THREADS -o quant_trinity/$S
done

# --- Trinity -> Ensembl genes ---
minimap2 -x asm5 -t $THREADS --secondary=no reference/ensembl_cdna.fa $TRINITY_FA > trinity_vs_ensembl.paf

# Best good match of each Trinity transcript, translated into an Ensembl gene
awk -F'\t' -v OFS='\t' '
  FNR == NR { gene[$1] = $2; next }
  { ident = $10 / $11; cover = ($4 - $3) / $2
    if (ident >= 0.95 && cover >= 0.5 && $10 > best[$1]) { best[$1] = $10; hit[$1] = gene[$6] } }
  END { for (t in hit) print t, hit[t] }' reference/ensembl_tx2gene.tsv trinity_vs_ensembl.paf > trinity_to_ensembl.tsv

# Transcript -> gene table for R: the Ensembl gene of the transcript, or else of another
# isoform of the same Trinity gene, or else the Trinity gene itself
awk -F'\t' -v OFS='\t' '
  FNR == 1 { f++ }
  f == 1 { ens[$1] = $2; next }
  f == 2 { if (($2 in ens) && !($1 in g2e)) g2e[$1] = ens[$2]; next }
  { print $2, (($2 in ens) ? ens[$2] : (($1 in g2e) ? g2e[$1] : $1)) }' \
  trinity_to_ensembl.tsv $TRINITY_FA.gene_trans_map $TRINITY_FA.gene_trans_map > trinity_tx2gene.tsv

N_TX=$(wc -l < trinity_tx2gene.tsv)
N_MATCHED=$(wc -l < trinity_to_ensembl.tsv)
N_ENS_FOUND=$(cut -f2 trinity_to_ensembl.tsv | sort -u | wc -l)
N_ENS_TOTAL=$(cut -f2 reference/ensembl_tx2gene.tsv | sort -u | wc -l)
N_TRI_ONLY=$(awk -F'\t' '$2 ~ /^TRINITY/ {print $2}' trinity_tx2gene.tsv | sort -u | wc -l)
echo "Trinity transcripts: $N_TX, matched to an Ensembl gene: $N_MATCHED"
echo "Ensembl genes recovered by Trinity: $N_ENS_FOUND of $N_ENS_TOTAL"
echo "Trinity genes without an Ensembl match: $N_TRI_ONLY"

# --- Unified report ---
multiqc qc_raw qc_clean clean quant_ensembl quant_trinity busco_*/short_summary*.txt -d -dd -1 -o multiqc_report -f

# --- Copy the results for R ---
mkdir -p $R_DIR
cp samples.tsv reference/ensembl_tx2gene.tsv reference/gene_annotation.tsv trinity_tx2gene.tsv $R_DIR/
for Q in quant_ensembl quant_trinity; do
  for S in $(tail -n +2 samples.tsv | cut -f1); do
    mkdir -p $R_DIR/$Q/$S
    cp $Q/$S/quant.sf $R_DIR/$Q/$S/
  done
done
ls $R_DIR $R_DIR/quant_ensembl
