ENV_NAME=bmsi_pc1
PKGS="python=3.11 fastqc=0.12.1 fastp=1.3.7 megahit=1.2.9 bwa=0.7.19 samtools=1.24 minimap2=2.31 quast=5.3.0 seqkit=2.14.0 multiqc=1.35 sra-tools=3.4.1 blast=2.16.0"

# 1. Which tools are missing?
MISSING=""
for cmd in fastqc fastp megahit bwa samtools minimap2 quast.py seqkit multiqc prefetch blastn; do
  command -v $cmd > /dev/null 2>&1 || MISSING="$MISSING $cmd"
done

if [ -z "$MISSING" ]; then
  echo "All tools are available: nothing to install."
else
  echo "Missing tools:$MISSING"

  # 2. Install conda (Miniforge) only if there is no conda at all
  if ! command -v conda > /dev/null 2>&1 && [ ! -d "$HOME/miniforge3" ]; then
    curl -L -o Miniforge3.sh \
      "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3.sh -b -p "$HOME/miniforge3"
    rm Miniforge3.sh
    "$HOME/miniforge3/bin/conda" init bash
  fi

  # 3. Make conda usable in this terminal
  if command -v conda > /dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
  else
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
  fi

  # 4. Macs with Apple chips use the Intel versions of the tools (through Rosetta)
  if [ "$(uname)-$(uname -m)" = "Darwin-arm64" ]; then
    export CONDA_SUBDIR=osx-64
  fi

  # 5. Create the environment, or add what is missing if it already exists
  if conda env list | grep -q "^$ENV_NAME "; then
    conda install -y -n $ENV_NAME --override-channels -c conda-forge -c bioconda $PKGS
  else
    conda create -y -n $ENV_NAME --override-channels -c conda-forge -c bioconda $PKGS
  fi
  conda activate $ENV_NAME
fi

fastqc --version; fastp --version; megahit --version; bwa 2>&1 | grep Version
samtools --version | head -1; minimap2 --version; quast.py --version; seqkit version; multiqc --version
blastn -version | head -1

mkdir -p ~/BMSi/PC1_assembly
cd ~/BMSi/PC1_assembly
pwd

ls -lh ForStudents

if ! command -v megahit > /dev/null 2>&1; then
  if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
  else
    eval "$(conda shell.bash hook)"
  fi
  conda activate bmsi_pc1
fi
# ---- Dataset and folder: the only lines to change for other data ----
WORKDIR=~/BMSi/PC1_assembly
SRR=SRR36009015
REF_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
# ----------------------------------------------------------------------
REF=$(basename $REF_URL .gz)
THREADS=$(getconf _NPROCESSORS_ONLN)
mkdir -p $WORKDIR && cd $WORKDIR
echo "Folder: $WORKDIR | Reads: $SRR | Reference: $REF | Threads: $THREADS"

if [ -d ForStudents ]; then
  mkdir -p megahit_output_final
  P=ForStudents
  [ -f $P/1_${SRR}_1.fastq ] && mv $P/1_${SRR}_1.fastq ${SRR}_1.fastq
  [ -f $P/1_${SRR}_2.fastq ] && mv $P/1_${SRR}_2.fastq ${SRR}_2.fastq
  [ -f $P/7_$REF ] && mv $P/7_$REF $REF
  [ -f $P/6_contigs_clean.fasta ] && mv $P/6_contigs_clean.fasta megahit_output_final/contigs_clean.fasta
  [ -f $P/8_ref_aligned_final.bam ] && mv $P/8_ref_aligned_final.bam ref_aligned_final.bam
  [ -f $P/8_ref_aligned_final.bam.bai ] && mv $P/8_ref_aligned_final.bam.bai ref_aligned_final.bam.bai
fi
ls -lh

if [ -s ${SRR}_1.fastq ] && [ -s ${SRR}_2.fastq ]; then
  echo "Raw reads found: skipping download."
else
  time prefetch -p $SRR
  time fasterq-dump -p -e $THREADS --split-files $SRR
  rm -rf $SRR
fi
ls -lh ${SRR}_*.fastq

if [ ! -s $REF ]; then
  wget -O ${REF}.gz $REF_URL
  gunzip ${REF}.gz
fi
samtools faidx $REF
cat ${REF}.fai

mkdir -p qc_raw
fastqc -t 2 -o qc_raw ${SRR}_1.fastq ${SRR}_2.fastq

fastp -i ${SRR}_1.fastq -I ${SRR}_2.fastq \
      -o ${SRR}_1_paired.fastq -O ${SRR}_2_paired.fastq \
      --detect_adapter_for_pe \
      --cut_front --cut_tail --cut_mean_quality 20 \
      --qualified_quality_phred 20 --unqualified_percent_limit 40 \
      --length_required 36 \
      --thread $THREADS \
      --html fastp.html --json fastp.json

mkdir -p qc_clean
fastqc -t 2 -o qc_clean ${SRR}_1_paired.fastq ${SRR}_2_paired.fastq

seqkit stats ${SRR}_1.fastq ${SRR}_2.fastq ${SRR}_1_paired.fastq ${SRR}_2_paired.fastq

if [ -s megahit_output_final/contigs_clean.fasta ]; then
  echo "Contigs found: skipping MEGAHIT."
else
  rm -rf megahit_output_final
  time megahit -1 ${SRR}_1_paired.fastq -2 ${SRR}_2_paired.fastq \
          -o megahit_output_final -t $THREADS --min-contig-len 1000
  seqkit seq -i -w 60 megahit_output_final/final.contigs.fa \
          > megahit_output_final/contigs_clean.fasta
fi

seqkit stats -a megahit_output_final/contigs_clean.fasta

if [ -s ref_aligned_final.bam ] && [ -s ref_aligned_final.bam.bai ]; then
  echo "BAM file found: skipping mapping."
else
  bwa index $REF
  time bwa mem -t $THREADS $REF ${SRR}_1_paired.fastq ${SRR}_2_paired.fastq \
    | samtools sort -@ $THREADS -o ref_aligned_final.bam
  samtools index ref_aligned_final.bam
fi

samtools flagstat ref_aligned_final.bam
samtools stats ref_aligned_final.bam > ref_aligned_final.stats
samtools coverage ref_aligned_final.bam

quast.py megahit_output_final/contigs_clean.fasta -r $REF \
         -o quast_report_megahit -t $THREADS

minimap2 -ax asm20 -t $THREADS $REF megahit_output_final/contigs_clean.fasta \
  | samtools sort -o contigs_vs_ref.bam
samtools index contigs_vs_ref.bam

samtools flagstat contigs_vs_ref.bam > contigs_vs_ref.flagstat
cat contigs_vs_ref.flagstat

Q=quast_report_megahit/report.tsv
q() { awk -F'\t' -v k="$1" '$1==k {print $2}' $Q; }

# Reference-guided: numbers from the reads mapped to the reference
MAPPED=$(samtools flagstat -O tsv ref_aligned_final.bam | awk -F'\t' '$3=="primary mapped %" {sub("%","",$1); print $1}')
read BREADTH DEPTH < <(samtools coverage -H ref_aligned_final.bam | awk '{print $6, $7}')
read LONGEST BLOCKS < <(samtools depth -a ref_aligned_final.bam | awk '
  $3>0 { run++; if (run==1) blocks++; if (run>max) max=run; next }
  { run=0 }
  END { print max, blocks }')
ERR=$(awk -F'\t' '$2=="error rate:" {printf "%.0f", $3*100000}' ref_aligned_final.stats)
UNMAPPED=$(awk -v m="$MAPPED" 'BEGIN {printf "%.1f", 100-m}')

# Write the table (tab-separated; the "#" lines let MultiQC show it as a table)
{
  echo "# id: 'strategy_comparison'"
  echo "# section_name: 'De novo vs. reference-guided'"
  echo "# description: 'The same clean reads analysed with the two strategies.'"
  echo "# plot_type: 'table'"
  printf 'Metric\tDe novo (MEGAHIT)\tReference-guided (BWA)\n'
  printf 'Needs a reference genome\tno\tyes\n'
  printf 'Number of pieces\t%s contigs\t%s covered blocks\n' "$(q '# contigs')" "$BLOCKS"
  printf 'Longest continuous sequence (bp)\t%s\t%s\n' "$(q 'Largest contig')" "$LONGEST"
  printf 'Reference genome recovered (%%)\t%s\t%s\n' "$(q 'Genome fraction (%)')" "$BREADTH"
  printf 'Differences from the reference (per 100 kbp)\t%s\t%s (incl. sequencing errors)\n' "$(q '# mismatches per 100 kbp')" "$ERR"
  printf 'Sequence absent from the reference (bp)\t%s\t0 (unmapped reads discarded: %s%%)\n' "$(q 'Unaligned length')" "$UNMAPPED"
  printf 'Reads used (%%)\t100 (all clean reads)\t%s\n' "$MAPPED"
  printf 'Mean depth (x)\t-\t%s\n' "$DEPTH"
} > strategy_comparison_mqc.tsv

grep -v '^#' strategy_comparison_mqc.tsv | awk -F'\t' '{printf "%-46s %-24s %s\n", $1, $2, $3}'

multiqc . --ignore ForStudents -o multiqc_report -f

# 1. Read coverage averaged in 1 kb windows (reference-guided)
samtools depth -a ref_aligned_final.bam | awk -v OFS='\t' -v W=1000 '
  { b = int(($2-1)/W); k = $1 SUBSEP b
    if (!(k in n)) order[++m] = k
    sum[k] += $3; n[k]++ }
  END { print "track type=bedGraph name=\"Read coverage (1 kb bins)\""
        for (i = 1; i <= m; i++) { split(order[i], a, SUBSEP); s = a[2]*W
          print a[1], s, s + n[order[i]], sum[order[i]] / n[order[i]] } }' > coverage_1kb.bedgraph

# 2. Position of each contig on the reference (de novo)
samtools view -F 260 contigs_vs_ref.bam | awk -v OFS='\t' '{
  len = 0; c = $6
  while (match(c, /[0-9]+[MDN=X]/)) { len += substr(c, RSTART, RLENGTH-1); c = substr(c, RSTART+RLENGTH) }
  strand = (int($2/16) % 2) ? "-" : "+"
  print $3, $4-1, $4-1+len, $1, 0, strand }' | sort -k2,2n > contigs_on_ref.bed

CONTIGS=megahit_output_final/contigs_clean.fasta
samtools faidx $CONTIGS

# The reference aligned onto the contigs: which parts of the contigs also exist in the reference
minimap2 -ax asm20 -t $THREADS $CONTIGS $REF | samtools sort -o ref_vs_contigs.bam
samtools index ref_vs_contigs.bam

# The reads mapped onto the contigs: the evidence that supports each contig
bwa index $CONTIGS
time bwa mem -t $THREADS $CONTIGS ${SRR}_1_paired.fastq ${SRR}_2_paired.fastq \
  | samtools sort -@ $THREADS -o reads_vs_contigs.bam
samtools index reads_vs_contigs.bam

# Contig(s) that do not align to the reference at all
samtools view -f 4 contigs_vs_ref.bam | cut -f1

NEW=$(samtools view -f 4 contigs_vs_ref.bam | cut -f1 | head -1)
samtools faidx $CONTIGS $NEW > ${NEW}.fasta
cat ${NEW}.fasta

time blastn -query ${NEW}.fasta -db nt -remote -evalue 1e-20 -max_target_seqs 10 \
  -outfmt "6 sacc pident length evalue stitle" | awk -F'\t' '!seen[$1]++' | head

