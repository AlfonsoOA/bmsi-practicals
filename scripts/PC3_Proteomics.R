# =============================================================================
# PC3_Proteomics.R
# Systems Molecular Biology - Practical 3 (Proteomics) - 2026-2027
#
# Label-free DIA proteomics (Spectronaut report): two complete pipelines
#   A. step by step: filters, median normalisation, data-driven check of an
#      intensity-dependent bias (cyclic loess), QC drift check, hybrid imputation
#      (random forest + QRILC), exploratory analysis, Welch t-test
#   B. MSstats: same filters, from fragments to proteins, exploratory analysis,
#      linear model
# then biological relevance of each list of significant proteins (fold change
# vs. Bayesian probability) and functional enrichment (ORA and GSEA).
#
# RESULTS (folder "results")
#   1_Data               samples, identifications per run
#   2_PipelineA_ttest    preprocessing, exploratory analysis, t-test
#   3_PipelineB_MSstats  MSstats settings, exploratory analysis, MSstats test
#   4_Relevance          t-test vs MSstats, FC and Bayesian relevance (files *_ttest_* / *_MSstats_*)
#   5_Enrichment         ORA and GSEA (files named after the list they come from)
#
# HOW TO USE
#   1. Check the CONFIGURATION block. With the files of ForStudents_PC3.zip,
#      only EXCLUDE_RUNS needs to be changed (second pass, guide 4.2).
#   2. Run the script section by section in RStudio: place the cursor in a
#      section and press Ctrl+Alt+T (Code > Run Region > Run Code Section).
#   3. All results are written to the folder "results" inside DATA_DIR.
# =============================================================================


# CONFIGURATION: the only part you may need to edit ---------------------------

# Folder with the input files. Empty = automatic:
#   Windows: C:/Users/<you>/BMSi/PC3_Proteomics   macOS/Linux: ~/BMSi/PC3_Proteomics
DATA_DIR <- ""

# Input files (inside DATA_DIR)
REPORT       <- "1_Spectronaut_report.tsv.gz"     # Spectronaut report, fragment level
INTERPRO     <- "2_InterProScanAnnotations.tsv"   # InterProScan output (tsv)
GHOSTKOALA   <- "3_GhostKOALAAnnotations.txt"     # protein -> KO number
FASTA        <- "4_organism_FASTA.fasta"          # proteome used for the search
DESCRIPTIONS <- "5_organisms_proteindescriptions.csv"  # written from the FASTA if missing
KEGG_MAP     <- "6_KEGG_KO_pathways.tsv"          # KO -> KEGG pathway (+ name)

# Conditions: log2FC = log2(TEST / REFERENCE). QC_LABEL = condition of the QC runs.
TEST      <- "Mutante"
REFERENCE <- "WT"
QC_LABEL  <- "QC"

# Runs to exclude (file names, column R.FileName). First pass: none. After the
# exploratory analysis (2.8 and 3.4; guide 4.2), write here the runs to leave out, e.g.:
#   EXCLUDE_RUNS <- c("SOL2023001F_xx_S2-xx_1_xxxxx", ...)
EXCLUDE_RUNS <- c()

# Injection order and sample number, read from the file name with these patterns
INJECTION_PATTERN <- "_([0-9]+)$"          # last number: SOL2023001F_7_S2-A9_1_12741 -> 12741
SAMPLE_PATTERN    <- "^[^_]+_([^_]+)_"     # second field: SOL2023001F_7_S2-A9_1_12741 -> 7

# Normalisation (2.3-2.4): median for all runs; then cyclic loess for ALL runs if an
# intensity-dependent bias is found: "auto" (decided from the data), "median", "loess"
NORMALISATION <- "auto"
DELTA       <- 0.1    # smallest intensity-dependent bias worth correcting (log2; 0.1 = 7%, a tenth of LFC)
CENSOR_TRIM <- 0.15   # fraction of the lowest intensities left out of the bias curves
NPERM       <- 100    # permutations for the noise level of the bias curves

# Protein-specific drift with the QC runs (2.5): "auto" (only if detected and useful), "always", "never"
DRIFT <- "auto"

# Filters and thresholds
QVALUE       <- 0.01   # identification FDR (precursor and protein)
MIN_PEPTIDES <- 2      # minimum number of peptides per protein (pipeline A)
MIN_VALID    <- 0.5    # keep a protein if observed in >= this fraction of the replicates of at least one condition
FDR          <- 0.05   # adjusted p-value for significance
LFC          <- 1      # |log2FC| threshold of relevance (1 = two-fold)
PROB         <- 0.9    # Bayesian probability of a relevant change
ORA_MIN      <- 3      # ORA: minimum number of listed proteins in a term (fewer cannot show an enrichment)

# MSstats settings compared in section 3, and the one carried to the next sections
MSSTATS_SETTINGS <- list(
  default            = list(process = list(normalization = "equalizeMedians", MBimpute = TRUE),  centre = FALSE,
                            label = "equalizeMedians + TMP + imputation of censored values (MSstats default)"),
  no_imputation      = list(process = list(normalization = "equalizeMedians", MBimpute = FALSE), centre = FALSE,
                            label = "equalizeMedians + TMP, no imputation"),
  quantile           = list(process = list(normalization = "quantile", MBimpute = TRUE),         centre = FALSE,
                            label = "quantile normalisation + TMP + imputation"),
  top3               = list(process = list(normalization = "equalizeMedians", MBimpute = TRUE, featureSubset = "top3"), centre = FALSE,
                            label = "equalizeMedians + top 3 features + imputation"),
  default_centred    = list(process = list(normalization = "equalizeMedians", MBimpute = TRUE),  centre = TRUE,
                            label = "default + protein-level median centring"),
  no_imputation_centred = list(process = list(normalization = "equalizeMedians", MBimpute = FALSE), centre = TRUE,
                            label = "no imputation + protein-level median centring")
)
MSSTATS_USE <- "default_centred"

# Functional enrichment: annotations to use (TRUE/FALSE). MetaCyc and Reactome are off: InterPro
# transfers them from protein families, so one generic domain (e.g. a P-loop ATPase) carries
# hundreds of "pathways" that are not pathways of this bacterium. Turn them on to see it.
ENRICH <- c(KEGG = TRUE, GO_BP = TRUE, GO_MF = TRUE, GO_CC = TRUE, InterPro = TRUE,
            Pfam = TRUE, PANTHER = TRUE, CDD = TRUE, SMART = TRUE, NCBIfam = TRUE,
            MetaCyc = FALSE, Reactome = FALSE)

CORES <- 0      # cores for the Bayesian step; 0 = automatic (all but one)
SEED  <- 2026   # makes imputation, permutations and MCMC reproducible


# Section 0. Packages ----------------------------------------------------------
# Installs the packages that are missing (only the first time) and loads them.

if (getRversion() < "4.3.0") {
  stop("This script needs R 4.3.0 or newer. You have R ", getRversion(), ".")
}
options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("data.table", "ggplot2", "ggrepel", "pheatmap", "vegan", "ape", "UpSetR", "VennDiagram",
          "missForest", "pcaMethods", "impute", "imputeLCMD", "rstan", "future.apply", "Matrix",
          "MSstats", "clusterProfiler", "enrichplot", "GO.db")
build_pkgs <- c("BH", "StanHeaders", "RcppEigen", "RcppParallel", "pkgbuild")   # needed to compile the Bayesian model (4.2)
missing_pkgs <- c(pkgs, build_pkgs)[!vapply(c(pkgs, build_pkgs), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)   # CRAN and Bioconductor
}
for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("Package check:\n")
for (p in pkgs) cat(sprintf("  %-16s %s\n", p, packageVersion(p)))
set.seed(SEED)


# Section 1. Data import -------------------------------------------------------
# Report, sample table (condition, sample, injection order), protein descriptions.

if (DATA_DIR == "") {
  home <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
  DATA_DIR <- file.path(home, "BMSi", "PC3_Proteomics")
}
DATA_DIR <- normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE)
needed <- c(REPORT, INTERPRO, GHOSTKOALA, FASTA)
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
if (!all(file.exists(file.path(DATA_DIR, needed)))) {
  stop("Put these files in ", DATA_DIR, " (guide 4.1) and run this section again:\n  ",
       paste(needed[!file.exists(file.path(DATA_DIR, needed))], collapse = "\n  "))
}
OUT <- file.path(DATA_DIR, "results")
out_dir <- function(name) { d <- file.path(OUT, name); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
key_hash <- function(key) {   # short code of a configuration: each one keeps its own result files
  f <- tempfile(); writeBin(charToRaw(enc2utf8(paste(key, collapse = "\n"))), f); h <- substr(unname(tools::md5sum(f)), 1, 8); unlink(f); h
}
save_plot <- function(p, folder, file, w = 7, h = 5.5) ggsave(file.path(folder, file), p, width = w, height = h, dpi = 300)
D1 <- out_dir("1_Data")
cat("Data folder:", DATA_DIR, "\nResults:    ", OUT, "\n")

# 1.1 Report (only the columns used)
cols <- c("R.Condition", "R.FileName", "R.Replicate", "PG.ProteinGroups", "PG.Qvalue", "PG.Quantity",
          "PEP.StrippedSequence", "EG.ModifiedSequence", "EG.Qvalue", "FG.Charge", "F.Charge",
          "F.FrgIon", "F.FrgLossType", "F.ExcludedFromQuantification", "F.PeakArea")
report <- fread(file.path(DATA_DIR, REPORT), select = cols, showProgress = FALSE)
cat("Report:", format(nrow(report), big.mark = ","), "rows,", uniqueN(report$PG.ProteinGroups), "protein groups\n")

# 1.2 Sample table
samples <- unique(report[, .(Run = R.FileName, Condition = R.Condition, Replicate = R.Replicate)])
samples[, Injection := as.integer(sub(paste0(".*", INJECTION_PATTERN), "\\1", Run))]
samples[, Sample := sub(paste0(SAMPLE_PATTERN, ".*"), "\\1", Run)]
samples[, Included := !(Run %in% EXCLUDE_RUNS)]
setorder(samples, Injection)
samples[, Order := seq_len(.N)]
samples[, Label := ifelse(Condition == QC_LABEL, paste0("QC_", Order), paste0(Condition, "_", Sample, "_", Order))]
if (any(!EXCLUDE_RUNS %in% samples$Run)) warning("Not in the report: ", paste(setdiff(EXCLUDE_RUNS, samples$Run), collapse = ", "))
dup <- samples[Condition != QC_LABEL & Included, .N, by = .(Condition, Sample)][N > 1]
if (nrow(dup) > 0) {
  cat("\nWARNING: the same sample was injected more than once (technical replicates):\n")
  print(samples[dup, on = .(Condition, Sample), .(Condition, Sample, Run, Injection)], row.names = FALSE)
  cat("Keep only one injection per sample (EXCLUDE_RUNS), or they count as independent replicates.\n")
}
fwrite(samples, file.path(D1, "Samples.csv"))
print(samples[, .(Label, Run, Condition, Injection, Included)], row.names = FALSE)
run_info <- samples[Included == TRUE]
report <- report[R.FileName %in% run_info$Run]
lab <- setNames(run_info$Label, run_info$Run)
cond <- setNames(run_info$Condition, run_info$Label)
cat("\nRuns analysed:", paste(names(table(run_info$Condition)), table(run_info$Condition), sep = " = ", collapse = ", "), "\n")

# 1.3 Protein descriptions from the FASTA headers (written once)
desc_file <- file.path(DATA_DIR, DESCRIPTIONS)
if (!file.exists(desc_file)) {
  h <- grep("^>", readLines(file.path(DATA_DIR, FASTA)), value = TRUE)
  descriptions <- data.table(ProteinID   = sub("^>[^|]*\\|([^|]+)\\|.*", "\\1", h),
                             Description = trimws(sub("^>\\S+\\s+(.*?)\\s+OS=.*$", "\\1", h)))
  fwrite(descriptions, desc_file)
}
descriptions <- fread(desc_file)
describe <- function(ids) {
  d <- descriptions$Description[match(sub(";.*", "", ids), descriptions$ProteinID)]
  ifelse(is.na(d), "", d)
}

# 1.4 Identifications per run (precursor and protein q-value < QVALUE)
ids <- report[EG.Qvalue < QVALUE & PG.Qvalue < QVALUE,
              .(Proteins = uniqueN(PG.ProteinGroups), Precursors = uniqueN(paste(EG.ModifiedSequence, FG.Charge))),
              by = R.FileName]
ids[, Label := lab[R.FileName]][, Condition := run_info$Condition[match(R.FileName, run_info$Run)]]
ids[, Label := factor(Label, levels = run_info$Label)]
fwrite(ids, file.path(D1, "Identifications_per_run.csv"))
p <- ggplot(ids, aes(Label, Proteins, fill = Condition)) + geom_col() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = "Run (in injection order)", y = "Protein groups identified")
save_plot(p, D1, "Identifications_per_run.png", 8, 5)
print(ids[order(Label), .(Label, Proteins, Precursors)], row.names = FALSE)


# Section 2. Pipeline A: preprocessing, exploratory analysis and t-test ----------
# Protein matrix -> 2-peptide filter -> median normalisation -> intensity-dependent bias
# (cyclic loess for all runs only if needed) -> QC drift check (correction only if needed
# and useful) -> hybrid imputation. The CV is followed after every step.
# Then exploratory analysis of the final matrix and Welch t-test.

D2 <- out_dir("2_PipelineA_ttest")
qc_runs  <- run_info[Condition == QC_LABEL, Label]
smp_runs <- run_info[Condition != QC_LABEL, Label]
inj <- setNames(run_info$Order, run_info$Label)
stages <- list()
protein_matrix <- function(rep) {   # log2 PG.Quantity, proteins x runs, protein q-value < QVALUE
  pg <- unique(rep[, .(R.FileName, PG.ProteinGroups, PG.Quantity, PG.Qvalue)])
  pg[PG.Qvalue >= QVALUE | !is.finite(PG.Quantity) | PG.Quantity <= 0, PG.Quantity := NA]
  pg[, Label := factor(lab[R.FileName], levels = run_info$Label)]
  M <- as.matrix(dcast(pg, PG.ProteinGroups ~ Label, value.var = "PG.Quantity", drop = FALSE), rownames = "PG.ProteinGroups")
  M <- log2(M[, run_info$Label])
  M[rowSums(!is.na(M)) > 0, ]
}
normalise <- function(M) sweep(M, 2, apply(M, 2, median, na.rm = TRUE)) + median(M, na.rm = TRUE)
median_cv <- function(M, runs) {   # median CV (%) per protein, on the intensity scale (not log)
  M <- 2^M[, intersect(runs, colnames(M)), drop = FALSE]
  v <- 100 * apply(M, 1, sd, na.rm = TRUE) / rowMeans(M, na.rm = TRUE)
  median(v[is.finite(v) & rowSums(!is.na(M)) >= 2])
}

# 2.1 Protein matrix as reported by Spectronaut (protein q-value < QVALUE)
X <- protein_matrix(report)
stages[["1 Raw"]] <- X

# 2.2 At least MIN_PEPTIDES peptides per protein IN EACH RUN. A protein identified in a
# run with fewer peptides loses everything it has in that run: its quantity here, and
# its fragments for MSstats (section 3).
npep <- report[EG.Qvalue < QVALUE, .(Peptides = uniqueN(PEP.StrippedSequence)), by = .(R.FileName, PG.ProteinGroups)]
p <- ggplot(npep, aes(pmin(Peptides, 30))) + geom_histogram(binwidth = 1, fill = "#2A6F97", colour = "white") +
  geom_vline(xintercept = MIN_PEPTIDES - 0.5, linetype = "dashed") + theme_bw() +
  labs(x = "Peptides of a protein in a run (30 = 30 or more)", y = "Protein-run pairs", title = "Peptides per protein and run")
save_plot(p, D2, "2.2_Peptides_per_protein_and_run.png", 7, 4.5)
good <- npep[Peptides >= MIN_PEPTIDES, .(R.FileName, PG.ProteinGroups)]
n_before <- sum(!is.na(X))
report <- report[good, on = .(R.FileName, PG.ProteinGroups), nomatch = NULL]
X <- protein_matrix(report)
cat("2-peptide filter: protein values kept", sum(!is.na(X)), "of", n_before, "| proteins", nrow(X), "of", nrow(stages[["1 Raw"]]), "\n")
stages[["2 Filtered (>=2 peptides/run)"]] <- X

# 2.3 Loading or drift? Both move whole runs up or down.
#   Loading: a run received more or less protein; nothing to do with its position in the sequence.
#   Drift: the instrument changes along the sequence; QC runs (identical pooled sample) show it.
# The run medians are split into a trend along the injection order, estimated from the QC runs
# (drift), and the rest (loading). Normalisation removes both when they affect all proteins equally.
run_med <- apply(X, 2, median, na.rm = TRUE)
qc_fit <- if (length(qc_runs) >= 3) lm(run_med[qc_runs] ~ inj[qc_runs]) else NULL
drift_trend <- if (is.null(qc_fit)) 0 * run_med else coef(qc_fit)[1] + coef(qc_fit)[2] * inj[names(run_med)]
shifts <- data.table(Run = names(run_med), Order = inj[names(run_med)], Condition = cond[names(run_med)],
                     Shift = run_med - median(run_med), Drift = drift_trend - mean(drift_trend))
shifts[, Loading := Shift - Drift]
fwrite(shifts, file.path(D2, "2.3_Run_shifts_loading_vs_drift.csv"))
if (length(qc_runs) >= 3) {
  cat(sprintf("Run medians: range %.2f log2 units. QC runs: range %.2f; drift along the sequence (QC trend): %.2f\n",
              diff(range(shifts$Shift)), diff(range(run_med[qc_runs])), diff(range(shifts$Drift))))
} else cat(sprintf("Run medians: range %.2f log2 units. Fewer than 3 QC runs: loading and drift cannot be separated.\n", diff(range(shifts$Shift))))
p <- ggplot(shifts, aes(Order, Shift, colour = Condition)) + geom_point(size = 3) +
  geom_line(aes(y = Drift), colour = "grey40", linetype = "dashed") + theme_bw() +
  labs(x = "Injection order", y = "Run median - overall median (log2)",
       title = "Run shifts: loading (points) and drift estimated from the QC runs (dashed line)")
save_plot(p, D2, "2.3_Loading_vs_drift.png", 8, 4.5)

# Median normalisation: centre the median of every run (removes loading and uniform drift)
X <- normalise(X)
stages[["3 Median normalised"]] <- X

# 2.4 Intensity-dependent bias: median only, or cyclic loess for ALL runs?
# Median normalisation moves whole runs up or down; it cannot remove a bias that depends on
# intensity (e.g. low-abundance proteins measured lower in some runs). Such a bias is a curve
# in the M-A plot of two runs (M = difference, A = mean intensity). How it is judged:
#   - Where: pairs of runs of the SAME condition (and pairs of QC runs). Any curve there is
#     technical; between conditions, a real change of many proteins would also bend it.
#   - The lowest intensities are left out (CENSOR_TRIM): missing values bend the curve there.
#   - Size: amplitude = range (5th to 95th percentile) of the fitted curve, in log2 units.
#   - Noise: the same amplitude after shuffling M (NPERM times): curves made by noise alone.
# One decision for the whole experiment (never run by run: runs normalised with different
# methods would differ in the shape of their distributions, an artefact). Cyclic loess if:
#   1. the curves within conditions are larger than noise, and
#   2. they are not larger between conditions than within (that could be biology: loess
#      would erase it), and
#   3. they are larger than DELTA. As DELTA is a choice, the decision is repeated from
#      DELTA/2 to 2*DELTA; if it changes in that range, a held-out test decides: the correction
#      is learned on half of the proteins and must shrink the curves of the other half
#      (in both halves). A real bias passes this test; a correction of noise does not.
ma_data <- function(x, y) {   # M and A of two runs, lowest intensities left out
  ok <- !is.na(x) & !is.na(y); m <- (x - y)[ok]; a <- ((x + y) / 2)[ok]
  k <- a >= quantile(a, CENSOR_TRIM); list(a = a[k], m = m[k])
}
amplitude <- function(a, m) { f <- lowess(a, m, f = 0.4, iter = 3)$y; unname(diff(quantile(f, c(0.05, 0.95)))) }
bias_pairs <- function(M, nperm = NPERM) {   # amplitude (and noise) of the curves of pairs of runs
  runs <- colnames(M); g <- cond[runs]
  within <- rbindlist(lapply(unique(g), function(x) if (sum(g == x) >= 2) as.data.table(t(combn(runs[g == x], 2)))))
  pr <- rbind(within[, Type := "within"], CJ(V1 = runs[g == REFERENCE], V2 = runs[g == TEST])[, Type := "between"])
  set.seed(SEED)
  res <- t(sapply(seq_len(nrow(pr)), function(i) {
    d <- ma_data(M[, pr$V1[i]], M[, pr$V2[i]])
    c(amplitude(d$a, d$m), if (nperm > 0) quantile(replicate(nperm, amplitude(d$a, sample(d$m))), 0.95) else NA)
  }))
  pr[, `:=`(Amplitude = res[, 1], Noise95 = res[, 2], Group = ifelse(Type == "within", cond[V1], "between conditions"))]
  pr
}
cyclic_loess <- function(M, train = rownames(M), iter = 3) {
  # "Fast" cyclic loess: every run against the average of all runs (pseudo-reference), with the
  # curve fitted only above the CENSOR_TRIM quantile and kept constant outside the fitted range.
  for (it in seq_len(iter)) {
    ref <- rowMeans(M, na.rm = TRUE)
    use <- rowSums(!is.na(M)) >= ncol(M) / 2 & rownames(M) %in% train
    for (j in colnames(M)) {
      d <- ma_data(M[use, j], ref[use])
      fit <- lowess(d$a, d$m, f = 0.4, iter = 3)
      M[, j] <- M[, j] - approx(fit$x, fit$y, xout = (M[, j] + ref) / 2, rule = 2, ties = mean)$y
    }
  }
  M
}
heldout <- function(M) {   # within-condition amplitude of the held-out half, before and after
  set.seed(SEED); fold <- sample(rep(1:2, length.out = nrow(M)))
  sapply(1:2, function(f) {
    test <- fold == f
    Mc <- cyclic_loess(M, train = rownames(M)[!test])
    c(before = median(bias_pairs(M[test, ], 0)[Type == "within", Amplitude]),
      after  = median(bias_pairs(Mc[test, ], 0)[Type == "within", Amplitude]))
  })
}
bias <- bias_pairs(X)
fwrite(bias, file.path(D2, "2.4_Intensity_bias_pairs.csv"))
amp_w <- median(bias[Type == "within", Amplitude]); noise_w <- median(bias[Type == "within", Noise95])
amp_b <- median(bias[Type == "between", Amplitude])
detected <- amp_w > noise_w
biology  <- amp_b > quantile(bias[Type == "within", Amplitude], 0.95)
grid <- data.table(DELTA = round(DELTA * c(0.5, 0.75, 1, 1.5, 2), 3))
grid[, Decision := ifelse(!detected | biology | amp_w < DELTA, "median", "cyclic loess")]
fwrite(grid, file.path(D2, "2.4_DELTA_sensitivity.csv"))
ho <- NULL
reason <- if (!detected) "the curves within conditions are not larger than noise" else
  if (biology) "WARNING: the curves are larger between conditions than within them (possible global biological change); loess would erase it" else
  if (uniqueN(grid$Decision) == 1) sprintf("decision stable for DELTA from %.3f to %.3f", min(grid$DELTA), max(grid$DELTA)) else {
    ho <- heldout(X)
    sprintf("decision depends on DELTA; held-out test: amplitude %.3f -> %.3f and %.3f -> %.3f",
            ho["before", 1], ho["after", 1], ho["before", 2], ho["after", 2])
  }
auto_loess <- detected && !biology && (if (is.null(ho)) amp_w >= DELTA else all(ho["after", ] < ho["before", ]))
use_loess <- switch(NORMALISATION, loess = TRUE, median = FALSE, auto = auto_loess)
norm_note <- c(sprintf("Intensity-dependent bias (median amplitude of the M-A curves, log2): within conditions %.3f (noise %.3f), between conditions %.3f; DELTA = %.3f",
                       amp_w, noise_w, amp_b, DELTA),
               paste("Reason:", reason),
               sprintf("Decision (NORMALISATION = \"%s\"): %s for all runs", NORMALISATION, if (use_loess) "median + cyclic loess" else "median only"))
cat(norm_note, sep = "\n")
writeLines(c(norm_note, "", capture.output(print(grid, row.names = FALSE))), file.path(D2, "2.4_Normalisation_decision.txt"))
p <- ggplot(bias, aes(Group, Amplitude)) + geom_point(aes(colour = "M-A curve"), size = 2.5, position = position_jitter(width = 0.1, seed = 1)) +
  geom_point(aes(y = Noise95, colour = "Noise (95th percentile)"), shape = 4, size = 2, position = position_jitter(width = 0.1, seed = 1)) +
  geom_hline(yintercept = DELTA, linetype = "dashed") + scale_colour_manual(values = c("M-A curve" = "#B5541B", "Noise (95th percentile)" = "grey50")) +
  guides(colour = guide_legend(override.aes = list(shape = c(16, 4)))) +
  theme_bw() + labs(x = NULL, y = "Amplitude of the M-A curve (log2)", colour = NULL,
                    title = "Intensity-dependent bias between pairs of runs", subtitle = paste("Dashed line: DELTA.", norm_note[3]))
save_plot(p, D2, "2.4_Intensity_bias_pairs.png", 8, 5)
if (use_loess) X <- cyclic_loess(X)
stages[["4 Intensity bias step"]] <- X
ma_curves <- function(M, stage) {   # every run against the average of all runs
  ref <- rowMeans(M, na.rm = TRUE)
  rbindlist(lapply(colnames(M), function(j) { d <- ma_data(M[, j], ref); f <- lowess(d$a, d$m, f = 0.4, iter = 3)
    data.table(Stage = stage, Run = j, Condition = cond[j], A = f$x, M = f$y) }))
}
cur <- rbind(ma_curves(stages[["3 Median normalised"]], "Median normalised"), ma_curves(X, if (use_loess) "+ cyclic loess" else "(no loess)"))
cur[, Stage := factor(Stage, levels = unique(Stage))]
p <- ggplot(cur, aes(A, M, group = Run, colour = Condition)) + geom_line() + geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Stage) + theme_bw() + labs(x = "A: mean log2 intensity", y = "M: run - average of all runs (log2)",
  title = "M-A curve of every run (lowest intensities left out)", subtitle = "Curves grouped by condition could be biology, not bias")
save_plot(p, D2, "2.4_MA_curves.png", 10, 4.5)


# 2.5 Protein-specific drift. After normalisation, a drift that is not the same for all
# proteins can remain (e.g. low-abundance proteins fading). For each protein seen in all
# QC runs, the slope of its QC intensity along the injection order is estimated; if the
# slopes depend on abundance (Spearman test), a per-protein correction is tried, and it
# is kept only if it lowers the CV of the biological replicates. (The CV of the QC runs
# cannot judge it: the correction is built from them, so it always makes them look better.)
drift_note <- "Fewer than 3 QC runs: drift not assessed."
if (length(qc_runs) >= 3) {
  full <- rowSums(is.na(X[, qc_runs])) == 0
  qo <- inj[qc_runs]
  slopes <- apply(X[full, qc_runs], 1, function(q) coef(lm(q ~ qo))[2])
  abund <- rowMeans(X[full, qc_runs])
  st <- cor.test(slopes, abund, method = "spearman", exact = FALSE)
  p <- ggplot(data.frame(abund, slopes), aes(abund, slopes)) + geom_point(size = 0.6, alpha = 0.4) +
    geom_smooth(method = "loess", formula = y ~ x, colour = "#B5541B") + geom_hline(yintercept = 0, linetype = "dashed") + theme_bw() +
    labs(x = "Mean log2 intensity in the QC runs", y = "Slope along the injection order (log2 per injection)",
         title = "Protein-specific drift after normalisation",
         subtitle = sprintf("Spearman rho = %.2f, p = %.3g (%d proteins seen in all QC runs)", st$estimate, st$p.value, sum(full)))
  save_plot(p, D2, "2.5_Drift_by_abundance.png", 7, 5)
  corrected <- X
  for (i in which(full)) {
    f <- coef(lm(X[i, qc_runs] ~ qo))
    corrected[i, ] <- X[i, ] - (f[1] + f[2] * inj[colnames(X)]) + mean(X[i, qc_runs])
  }
  corrected <- normalise(corrected)
  bio <- c(median_cv(X, smp_runs[cond[smp_runs] == REFERENCE]), median_cv(X, smp_runs[cond[smp_runs] == TEST]))
  bio_c <- c(median_cv(corrected, smp_runs[cond[smp_runs] == REFERENCE]), median_cv(corrected, smp_runs[cond[smp_runs] == TEST]))
  use_it <- switch(DRIFT, always = TRUE, never = FALSE, auto = st$p.value < 0.05 && mean(bio_c) < mean(bio))
  drift_note <- sprintf(paste0("Abundance-dependent drift: Spearman p = %.3g. Median CV of the replicates (%s, %s): %.1f%%, %.1f%% ",
                               "without correction; %.1f%%, %.1f%% with it. Decision (DRIFT = \"%s\"): %s."),
                        st$p.value, REFERENCE, TEST, bio[1], bio[2], bio_c[1], bio_c[2], DRIFT,
                        if (use_it) "correction APPLIED" else "correction NOT applied")
  if (use_it) X <- corrected
}
cat(drift_note, "\n")
writeLines(drift_note, file.path(D2, "2.5_Drift_decision.txt"))
stages[["5 Drift step"]] <- X

# 2.6 Hybrid imputation (biological samples only; QC runs are no longer needed)
Xs <- X[, smp_runs]
grp <- cond[smp_runs]
obs <- sapply(c(REFERENCE, TEST), function(g) rowSums(!is.na(Xs[, grp == g, drop = FALSE])))
nrep <- table(grp)[c(REFERENCE, TEST)]
keep <- apply(sweep(obs, 2, nrep, "/") >= MIN_VALID, 1, any)
Xs <- Xs[keep, ]; obs <- obs[keep, ]
is_mec <- matrix(FALSE, nrow(Xs), ncol(Xs), dimnames = dimnames(Xs))
for (g in c(REFERENCE, TEST)) is_mec[obs[, g] == 0, grp == g] <- TRUE
is_pov <- is.na(Xs) & !is_mec
A_obs <- Xs   # before imputation: used by the Bayesian model (4.2), which models missing values
cat("Proteins kept (observed in >=", MIN_VALID * 100, "% of the replicates of a condition):", nrow(Xs), "\n")
cat("Missing values: POV =", sum(is_pov), " MEC =", sum(is_mec), " (of", length(Xs), "values)\n")
set.seed(SEED)
Xq <- imputeLCMD::impute.QRILC(Xs, tune.sigma = 1)[[1]]          # left-censored: MEC
Xi <- Xs; Xi[is_mec] <- Xq[is_mec]
if (any(is_pov)) {
  set.seed(SEED)
  Xi <- as.matrix(missForest::missForest(as.data.frame(Xi), verbose = FALSE)$ximp)   # MAR: POV
}
dimnames(Xi) <- dimnames(Xs)
stages[["6 Imputed"]] <- Xi
A <- Xi
imputed <- ifelse(is_mec, "MEC", ifelse(is_pov, "POV", "measured"))
fwrite(data.table(Protein = rownames(A), Description = describe(rownames(A)), round(A, 4)), file.path(D2, "2.6_ttest_final_matrix.csv"))
fwrite(data.table(Protein = rownames(A), imputed), file.path(D2, "2.6_imputed_values.csv"))

# 2.7 What every step did: number of proteins, distributions and CV
cv <- function(M, runs) {   # CV (%) per protein on the intensity scale (not log)
  M <- 2^M[, intersect(runs, colnames(M)), drop = FALSE]
  if (ncol(M) < 2) return(numeric(0))
  v <- 100 * apply(M, 1, sd, na.rm = TRUE) / rowMeans(M, na.rm = TRUE)
  v[is.finite(v) & rowSums(!is.na(M)) >= 2]
}
groups <- list(qc_runs, names(grp)[grp == REFERENCE], names(grp)[grp == TEST])
names(groups) <- c("QC", REFERENCE, TEST)
cv_long <- rbindlist(lapply(names(stages), function(s) rbindlist(lapply(names(groups), function(g) {
  v <- cv(stages[[s]], groups[[g]]); if (length(v)) data.table(Stage = s, Group = g, CV = v) }))))
cv_long[, Group := factor(Group, levels = names(groups))]
cv_sum <- cv_long[, .(Proteins = .N, Median_CV = round(median(CV), 1)), by = .(Stage, Group)]
fwrite(cv_sum, file.path(D2, "2.7_CV_by_stage.csv"))
print(dcast(cv_sum, Stage ~ Group, value.var = "Median_CV"), row.names = FALSE)
p <- ggplot(cv_long, aes(Stage, CV, fill = Group)) + geom_boxplot(outlier.size = 0.3) +
  coord_cartesian(ylim = c(0, 150)) + theme_bw() + theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(x = NULL, y = "CV (%) per protein", title = "Coefficient of variation after every step (pipeline A)",
       subtitle = "QC: technical variation only; conditions: technical + biological")
save_plot(p, D2, "2.7_CV_by_stage.png", 9, 5)
dist_long <- rbindlist(lapply(names(stages), function(s) {
  M <- stages[[s]]; data.table(Stage = s, Run = rep(colnames(M), each = nrow(M)), log2 = as.vector(M))
}))[!is.na(log2)]
dist_long[, Run := factor(Run, levels = run_info$Label)][, Condition := cond[as.character(Run)]]
p <- ggplot(dist_long, aes(Run, log2, fill = Condition)) + geom_boxplot(outlier.size = 0.2) + facet_wrap(~ Stage, ncol = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(x = "Run (in injection order)", y = "log2 intensity", title = "Distributions after every step")
save_plot(p, D2, "2.7_Distributions_by_stage.png", 9, 11)
med <- rbindlist(lapply(names(stages)[1:5], function(s) data.table(Stage = s, Run = colnames(stages[[s]]),
       Median = apply(stages[[s]], 2, median, na.rm = TRUE))))
med[, Order := inj[Run]][, Condition := cond[Run]]
p <- ggplot(med, aes(Order, Median, colour = Condition)) + geom_line(aes(group = 1), colour = "grey70") + geom_point(size = 3) +
  facet_wrap(~ Stage, nrow = 1) + theme_bw() + labs(x = "Injection order", y = "Median log2 intensity", title = "Run medians along the sequence")
save_plot(p, D2, "2.7_Medians_by_injection_order.png", 14, 4)
cat("Imputed values:", sum(is_pov), "POV (random forest) and", sum(is_mec), "MEC (QRILC)\n")



# 2.8 Exploratory analysis of the final matrix --------------------------------------
# The same function is used for pipeline B (3.4), so both can be compared figure by figure.
min_p <- function(g) { n <- table(g); s <- choose(sum(n), n[1]) / ifelse(n[1] == n[2], 2, 1); 1 / s }
explore <- function(M, present, folder, prefix, title) {
  # M: proteins x samples, no missing values; present: proteins observed in each condition
  meta <- data.frame(Condition = factor(cond[colnames(M)], levels = c(REFERENCE, TEST)), row.names = colnames(M))
  f <- function(x) file.path(folder, paste0(prefix, "_", x))
  # Proteins per condition (observed, before imputation): Venn diagram
  invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))
  venn <- VennDiagram::venn.diagram(present, filename = NULL, fill = c("#2A6F97", "#B5541B"), alpha = 0.5, disable.logging = TRUE)
  png(f("Venn_proteins.png"), width = 5, height = 5, units = "in", res = 300); grid::grid.draw(venn); invisible(dev.off())
  # Sample correlation
  pheatmap(cor(M), annotation_col = meta, annotation_row = meta, display_numbers = TRUE, number_format = "%.2f",
           main = paste(title, "- Pearson correlation between samples"), filename = f("Sample_correlation.png"), width = 9, height = 8)
  # PCA, PCoA, NMDS with PERMANOVA, PERMDISP and ANOSIM
  d <- dist(t(M))
  set.seed(SEED)
  perm <- vegan::adonis2(d ~ Condition, data = meta, permutations = 999)
  disp <- vegan::permutest(vegan::betadisper(d, meta$Condition), permutations = 999)
  anos <- vegan::anosim(d, meta$Condition, permutations = 999)
  stats_txt <- sprintf("PERMANOVA R2 = %.2f, p = %.3f | PERMDISP p = %.3f | ANOSIM R = %.2f, p = %.3f | min. p = %.3f",
                       perm$R2[1], perm$`Pr(>F)`[1], disp$tab$`Pr(>F)`[1], anos$statistic, anos$signif, min_p(meta$Condition))
  ord_plot <- function(df, x, y, t, xl, yl) {
    ggplot(df, aes(.data[[x]], .data[[y]], colour = Condition, label = Sample)) + geom_point(size = 3) +
      ggrepel::geom_text_repel(size = 3, show.legend = FALSE) + theme_bw() +
      labs(title = paste(title, "-", t), subtitle = stats_txt, x = xl, y = yl) + theme(plot.subtitle = element_text(size = 7))
  }
  pca <- prcomp(t(M), scale. = TRUE); ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df <- data.frame(pca$x[, 1:2], meta, Sample = colnames(M))
  ggsave(f("PCA.png"), ord_plot(df, "PC1", "PC2", "PCA", paste0("PC1 (", ve[1], "%)"), paste0("PC2 (", ve[2], "%)")), width = 7, height = 5.5, dpi = 300)
  pco <- ape::pcoa(d); re <- round(100 * pco$values$Relative_eig[1:2], 1)
  df <- data.frame(Axis1 = pco$vectors[, 1], Axis2 = pco$vectors[, 2], meta, Sample = colnames(M))
  ggsave(f("PCoA.png"), ord_plot(df, "Axis1", "Axis2", "PCoA (Euclidean distances)", paste0("PCoA 1 (", re[1], "%)"), paste0("PCoA 2 (", re[2], "%)")),
         width = 7, height = 5.5, dpi = 300)
  set.seed(SEED)
  nm <- vegan::metaMDS(d, k = 2, trymax = 50, trace = 0)
  df <- data.frame(nm$points, meta, Sample = colnames(M))
  ggsave(f("NMDS.png"), ord_plot(df, "MDS1", "MDS2", paste0("NMDS (stress = ", round(nm$stress, 3), ")"), "NMDS 1", "NMDS 2"), width = 7, height = 5.5, dpi = 300)
  writeLines(c(paste(title, "-", nrow(M), "proteins"), stats_txt, paste("NMDS stress:", round(nm$stress, 3)), "",
               capture.output(perm), "", capture.output(disp), "", capture.output(anos)), f("Multivariate_tests.txt"))
  # Heatmap of all proteins (row z-scores)
  pheatmap(M, scale = "row", annotation_col = meta, show_rownames = TRUE, fontsize_row = 2, clustering_method = "ward.D2",
           main = paste(title, "-", nrow(M), "proteins (row z-scores)"), filename = f("Heatmap_all_proteins.png"), width = 8, height = 14)
  cat(title, ":", stats_txt, "\n")
}
present_A <- lapply(c(REFERENCE, TEST), function(g) rownames(obs)[obs[, g] > 0]); names(present_A) <- c(REFERENCE, TEST)
explore(A, present_A, D2, "2.8", "t-test pipeline")

# 2.9 Welch t-test + Benjamini-Hochberg -------------------------------------------------
volcano <- function(r, title) {
  r <- r[is.finite(log2FC) & !is.na(adj.pvalue)]
  r[, Status := ifelse(Significant, ifelse(abs(log2FC) >= LFC, "Significant, |log2FC| >= LFC", "Significant, small change"), "Not significant")]
  top <- head(r[Significant == TRUE][order(adj.pvalue)], 12)
  ggplot(r, aes(log2FC, -log10(adj.pvalue), colour = Status)) + geom_point(size = 0.9, alpha = 0.7) +
    geom_hline(yintercept = -log10(FDR), linetype = "dashed") + geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top, aes(label = Protein), size = 2.5, colour = "black", max.overlaps = 30) +
    scale_colour_manual(values = c("Significant, |log2FC| >= LFC" = "#B5541B", "Significant, small change" = "#2A6F97", "Not significant" = "grey70")) +
    theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(title = title, x = paste0("log2FC (", TEST, " / ", REFERENCE, ")"), y = "-log10 adjusted p-value")
}
significance <- function(r, folder, prefix, tag, title) {   # tables (all and significant proteins) and volcano plot
  r[, Description := describe(Protein)]
  r[, Significant := !is.na(adj.pvalue) & adj.pvalue < FDR]
  r[, Direction := ifelse(log2FC > 0, "up", "down")]
  r <- r[order(adj.pvalue)]
  fwrite(r, file.path(folder, paste0(prefix, "_", tag, "_all_proteins.csv")))
  fwrite(r[Significant == TRUE], file.path(folder, paste0(prefix, "_", tag, "_significant_proteins.csv")))
  save_plot(volcano(r, title), folder, paste0(prefix, "_", tag, "_volcano.png"), 7, 6.5)
  cat(title, ": proteins tested", sum(!is.na(r$adj.pvalue)), "| significant", sum(r$Significant),
      "(up", sum(r$Significant & r$log2FC > 0), ", down", sum(r$Significant & r$log2FC < 0), ")\n")
  r
}
ref_cols <- names(grp)[grp == REFERENCE]; test_cols <- names(grp)[grp == TEST]
tt <- t(apply(A, 1, function(x) {
  r <- tryCatch(t.test(x[test_cols], x[ref_cols], var.equal = FALSE), error = function(e) NULL)
  if (is.null(r)) c(NA, NA, NA) else c(mean(x[test_cols]) - mean(x[ref_cols]), r$statistic, r$p.value)
}))
resA <- data.table(Protein = rownames(A), log2FC = tt[, 1], Tvalue = tt[, 2], pvalue = tt[, 3])
resA[, adj.pvalue := p.adjust(pvalue, "BH")]
resA[, MEC := rowSums(is_mec[Protein, , drop = FALSE]) > 0]
results <- list()
results$ttest <- significance(resA, D2, "2.9", "ttest", "Pipeline A: Welch t-test")


# Section 3. Pipeline B: MSstats -------------------------------------------------
# Preprocessing (shared with pipeline A: q-value and >=2 peptides per protein and run,
# applied to the report in 2.2), from fragments to proteins (dataProcess) and linear
# model (groupComparison) with several settings, exploratory analysis and significance.
# A healthy comparison has log2FC centred on 0: most proteins do not change.
# LONG STEP (about 3-5 min per setting): results are saved and reused while the
# analysed runs do not change.

D3 <- out_dir("3_PipelineB_MSstats")

# 3.1 Preprocessing: the report filtered in 2.2 (biological runs only) in MSstats format
ms_runs <- sort(run_info[Condition != QC_LABEL, Run])
key_b <- paste(c("v3", ms_runs, MIN_PEPTIDES, QVALUE), collapse = ";")
cat("MSstats input: the report filtered in 2.2 (q-value <", QVALUE, ", >=", MIN_PEPTIDES, "peptides per protein and run),",
    length(ms_runs), "biological runs\n")

# 3.2 dataProcess + groupComparison with every setting (cached)
centre_runs <- function(pl) {   # protein-level median centring of every run
  m <- tapply(pl$LogIntensities, pl$RUN, median, na.rm = TRUE)
  pl$LogIntensities <- pl$LogIntensities - m[as.character(pl$RUN)] + median(m)
  pl
}
cache_b <- file.path(D3, paste0("MSstats_results_", key_hash(key_b), ".rds"))   # one file per set of runs
old <- file.path(D3, "MSstats_results.rds")   # name used by earlier versions of the script
if (!file.exists(cache_b) && file.exists(old) && identical(readRDS(old)$key, key_b)) file.rename(old, cache_b)
ms <- if (file.exists(cache_b)) readRDS(cache_b) else list()
if (!identical(ms$key, key_b)) ms <- list(key = key_b)
todo <- setdiff(names(MSSTATS_SETTINGS), names(ms))
if (length(todo) > 0) {
  t0 <- Sys.time()
  invisible(capture.output(ms_input <- SpectronauttoMSstatsFormat(as.data.frame(report[R.Condition != QC_LABEL]),
    qvalue_cutoff = QVALUE, use_log_file = FALSE, verbose = FALSE)))
  # BioReplicate comes from R.Replicate, which restarts in each condition (a "2" in WT and a
  # "2" in Mutante). MSstats would read that as the same subject measured in both conditions
  # (repeated measures) and fit a paired model with 0-2 degrees of freedom. Each biological
  # sample is a different subject: make the labels unique.
  ms_input$BioReplicate <- paste(ms_input$Condition, ms_input$BioReplicate, sep = "_")
  processed <- list()
  for (v in todo) {
    set <- MSSTATS_SETTINGS[[v]]
    base <- paste(deparse(set$process), collapse = "")
    if (is.null(processed[[base]])) {
      cat("MSstats dataProcess:", v, "...\n")
      invisible(capture.output(processed[[base]] <- do.call(dataProcess,
        c(list(ms_input, summaryMethod = "TMP", use_log_file = FALSE, verbose = FALSE), set$process))))
    }
    pr <- processed[[base]]
    if (isTRUE(set$centre)) pr$ProteinLevelData <- centre_runs(pr$ProteinLevelData)
    contrast <- matrix(c(-1, 1), nrow = 1, dimnames = list(paste(TEST, "-", REFERENCE), c(REFERENCE, TEST)))
    contrast <- contrast[, levels(pr$ProteinLevelData$GROUP), drop = FALSE]
    invisible(capture.output(cmp <- groupComparison(contrast, pr, use_log_file = FALSE, verbose = FALSE)))
    ms[[v]] <- list(protein = as.data.table(pr$ProteinLevelData)[, .(Protein = as.character(Protein), Run = originalRUN,
                                                                     Condition = GROUP, LogIntensities,
                                                                     LogObserved = ifelse(NumMeasuredFeature > 0, LogIntensities, NA_real_))],
                    result = as.data.table(cmp$ComparisonResult))
    saveRDS(ms, cache_b)
  }
  cat("MSstats done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
} else cat("MSstats results found for these runs: reading them.\n")

# 3.3 Compare the settings: is log2FC centred on 0? are up and down balanced?
ms_summary <- rbindlist(lapply(names(MSSTATS_SETTINGS), function(v) {
  r <- ms[[v]]$result; f <- is.finite(r$log2FC); s <- !is.na(r$adj.pvalue) & r$adj.pvalue < FDR
  data.table(Setting = v, Description = MSSTATS_SETTINGS[[v]]$label, Tested = sum(!is.na(r$adj.pvalue)),
             Median_log2FC = round(median(r$log2FC[f]), 2), Significant = sum(s),
             Up = sum(s & r$log2FC > 0), Down = sum(s & r$log2FC < 0))
}))
ms_summary[, Used := Setting == MSSTATS_USE]
fwrite(ms_summary, file.path(D3, "3.3_MSstats_settings_compared.csv"))
print(ms_summary[, .(Setting, Median_log2FC, Significant, Up, Down, Used)], row.names = FALSE)
lf <- rbindlist(lapply(names(MSSTATS_SETTINGS), function(v) ms[[v]]$result[is.finite(log2FC), .(Setting = v, log2FC)]))
p <- ggplot(lf, aes(log2FC, colour = Setting)) + geom_density() + geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-4, 4)) + theme_bw() + labs(title = "log2FC of all proteins with each MSstats setting",
  subtitle = "Most proteins should not change: the curve should peak at 0", x = paste0("log2FC (", TEST, " / ", REFERENCE, ")"))
save_plot(p, D3, "3.3_MSstats_settings_log2FC.png", 8, 5)
pm <- rbindlist(lapply(names(MSSTATS_SETTINGS), function(v) ms[[v]]$protein[, .(Median = median(LogIntensities, na.rm = TRUE)), by = .(Run, Condition)][, Setting := v]))
pm[, Label := factor(lab[Run], levels = run_info$Label)]
p <- ggplot(pm, aes(Label, Median, colour = Condition)) + geom_point(size = 2.5) + facet_wrap(~ Setting, scales = "free_y") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) + labs(x = NULL, y = "Median protein log2 intensity", title = "Protein-level run medians with each setting")
save_plot(p, D3, "3.3_MSstats_settings_run_medians.png", 11, 7)
B_prot <- ms[[MSSTATS_USE]]$protein
cat("MSstats setting used:", MSSTATS_USE, "\n")

# 3.4 Exploratory analysis of the MSstats protein-level data (setting MSSTATS_USE)
Bm <- dcast(B_prot[is.finite(LogIntensities)], Protein ~ Run, value.var = "LogIntensities")
Bm <- as.matrix(Bm, rownames = "Protein")
colnames(Bm) <- lab[colnames(Bm)]
Bm <- Bm[, intersect(smp_runs, colnames(Bm))]
cat("MSstats protein matrix:", nrow(Bm), "proteins;", sum(complete.cases(Bm)), "with a value in every run (used here)\n")
fwrite(data.table(Protein = rownames(Bm), Description = describe(rownames(Bm)), round(Bm, 4)), file.path(D3, "3.4_MSstats_protein_matrix.csv"))
seen <- unique(report[R.Condition != QC_LABEL & is.finite(F.PeakArea) & F.PeakArea > 0 & EG.Qvalue < QVALUE & PG.Qvalue < QVALUE,
                      .(PG.ProteinGroups, R.Condition)])
present_B <- lapply(c(REFERENCE, TEST), function(g) seen[R.Condition == g, PG.ProteinGroups]); names(present_B) <- c(REFERENCE, TEST)
explore(Bm[complete.cases(Bm), ], present_B, D3, "3.4", "MSstats pipeline")

# 3.5 Significance (MSstats linear model, Benjamini-Hochberg)
B <- ms[[MSSTATS_USE]]$result[, .(Protein = as.character(Protein), log2FC, SE, Tvalue, pvalue, adj.pvalue, issue)]
cat("Proteins missing in one condition (infinite log2FC):", sum(B$issue %in% "oneConditionMissing", na.rm = TRUE), "\n")
results$MSstats <- significance(B, D3, "3.5", "MSstats", "Pipeline B: MSstats")


# Section 4. Biological relevance: fold change and Bayesian probability ----------
# First the two lists of significant proteins are compared; then each one is filtered
# with two criteria of relevance. File names say where each list comes from (ttest or
# MSstats) and which criterion was used (FC or Bayes).
# LONG STEP: one Bayesian model per significant protein (Stan). The result is
# saved and reused while the proteins, the samples and the settings do not change.

D4 <- out_dir("4_Relevance")

# 4.1 t-test vs MSstats: significant proteins and log2FC
sigA <- results$ttest[Significant == TRUE, Protein]; sigB <- results$MSstats[Significant == TRUE, Protein]
both <- intersect(results$ttest$Protein, results$MSstats$Protein)
lfc_ab <- merge(results$ttest[, .(Protein, A = log2FC)], results$MSstats[is.finite(log2FC), .(Protein, B = log2FC)], by = "Protein")
sig_summary <- data.table(
  Item = c("Proteins tested", "Significant", "Significant, up", "Significant, down", "Tested in both pipelines",
           "Significant in both", "Significant only with the t-test", "Significant only with MSstats", "Correlation of log2FC"),
  ttest_MSstats = c(paste(nrow(results$ttest), "/", sum(!is.na(results$MSstats$adj.pvalue))), paste(length(sigA), "/", length(sigB)),
            paste(sum(results$ttest[Significant == TRUE]$log2FC > 0), "/", sum(results$MSstats[Significant == TRUE]$log2FC > 0)),
            paste(sum(results$ttest[Significant == TRUE]$log2FC < 0), "/", sum(results$MSstats[Significant == TRUE]$log2FC < 0)),
            length(both), length(intersect(sigA, sigB)), length(setdiff(sigA, sigB)), length(setdiff(sigB, sigA)),
            round(cor(lfc_ab$A, lfc_ab$B), 3)))
fwrite(sig_summary, file.path(D4, "4.1_Significance_ttest_vs_MSstats.csv"))
print(sig_summary, row.names = FALSE)
p <- ggplot(lfc_ab, aes(A, B)) + geom_point(size = 0.6, alpha = 0.5) + geom_abline(colour = "#B5541B") + theme_bw() +
  labs(x = "log2FC, t-test pipeline", y = "log2FC, MSstats pipeline", title = "Same proteins, two pipelines")
save_plot(p, D4, "4.1_log2FC_ttest_vs_MSstats.png", 6, 6)

# 4.2 Bayesian probability of a relevant change for the significant proteins
# One model per protein (Stan), with three pieces of information taken from ALL proteins:
#   - Missing values are not imputed: each run has a dropout curve (probability of not seeing
#     a protein as a function of its intensity, probit), and a missing value contributes the
#     probability of being missing. A protein absent in a whole condition still gets a posterior.
#   - sigma (one per condition, as in Welch's test): prior centred on the trend of log(SD) with
#     intensity and number of peptides, with the spread that remains around that trend.
#   - log2FC: Student-t(3, 0, s), s = robust spread of the log2FC of all proteins (most proteins
#     do not change; the heavy tail lets large changes through).
#   - The mean of each condition: distribution of the intensities of all proteins. It bounds the
#     log2FC of a protein absent in a whole condition (otherwise only "below detection" is known).
#   - If the dropout curve of a run cannot be estimated (no clear link between intensity and
#     missingness), missing values are taken as censored below the 5th percentile of that run.
# Pipeline t-test: normalised matrix BEFORE imputation. Pipeline MSstats: protein values,
# missing when MSstats measured no feature of the protein in that run (NumMeasuredFeature = 0).
# The model is compiled once (needs Rtools on Windows) and only if some protein must be fitted.
stan_code <- "
data {
  int<lower=0> N_obs;  vector[N_obs] y;  array[N_obs] int<lower=1, upper=2> g_obs;
  int<lower=0> N_mis;  array[N_mis] int<lower=1, upper=2> g_mis;  vector[N_mis] rho;  vector<lower=0>[N_mis] zeta;
  real mu_loc;  real<lower=0> mu_scale;  real<lower=0> s_beta;  real logsig_loc;  real<lower=0> logsig_scale;
}
parameters { real mu;  real beta;  vector[2] log_sigma; }
transformed parameters {
  vector[2] m = [mu, mu + beta]';
  vector[2] sigma = exp(log_sigma);
}
model {
  m ~ normal(mu_loc, mu_scale);   // both condition means: distribution of protein intensities
  beta ~ student_t(3, 0, s_beta);
  log_sigma ~ normal(logsig_loc, logsig_scale);
  y ~ normal(m[g_obs], sigma[g_obs]);
  for (i in 1:N_mis)   // P(missing) = Phi((rho - m) / sqrt(sigma^2 + zeta^2))
    target += normal_lcdf(rho[i] | m[g_mis[i]], sqrt(square(sigma[g_mis[i]]) + square(zeta[i])));
}
"
stan_file <- file.path(D4, "protein_model.stan")
if (!file.exists(stan_file) || !identical(readLines(stan_file), strsplit(stan_code, "\n")[[1]])) writeLines(stan_code, stan_file)
pep_n <- report[EG.Qvalue < QVALUE, .(Peptides = uniqueN(PEP.StrippedSequence)), by = PG.ProteinGroups]
pep_n <- setNames(pep_n$Peptides, pep_n$PG.ProteinGroups)
hyper <- function(M) {   # M: proteins x biological runs, NA = missing. Returns the empirical priors.
  g <- grp[colnames(M)]
  pm <- rowMeans(M, na.rm = TRUE)
  dropout <- t(sapply(colnames(M), function(j) {   # P(missing) = Phi(a + b * mean) = Phi((rho - mean) / zeta)
    f <- suppressWarnings(glm(is.na(M[, j]) ~ pm, family = binomial(link = "probit")))
    b <- coef(f)[2]; a <- coef(f)[1]
    if (!is.finite(b) || b > -0.2) c(rho = unname(quantile(M[, j], 0.05, na.rm = TRUE)), zeta = 1, fitted = 0)   # zeta > 5: not reliable
    else c(rho = unname(-a / b), zeta = unname(-1 / b), fitted = 1)
  }))
  sds <- rbindlist(lapply(c(REFERENCE, TEST), function(x) {
    S <- M[, g == x, drop = FALSE]; n <- rowSums(!is.na(S))
    data.table(Protein = rownames(M), n = n, logsd = log(apply(S, 1, sd, na.rm = TRUE)), mean = pm)
  }))[n >= 3 & is.finite(logsd)]
  sds[, Peptides := log(pmax(1, pep_n[Protein]))]
  sds[is.na(Peptides), Peptides := log(MIN_PEPTIDES)]
  sds[, logsd_c := logsd - (digamma((n - 1) / 2) - log((n - 1) / 2)) / 2]   # bias of log(sample SD)
  trend <- lm(logsd_c ~ mean + Peptides, data = sds)
  tau <- sqrt(max(var(residuals(trend)) - mean(trigamma((sds$n - 1) / 2) / 4), 0.1^2))   # minus sampling noise
  both <- rowSums(!is.na(M[, g == REFERENCE, drop = FALSE])) >= 2 & rowSums(!is.na(M[, g == TEST, drop = FALSE])) >= 2
  lfc <- rowMeans(M[both, g == TEST, drop = FALSE], na.rm = TRUE) - rowMeans(M[both, g == REFERENCE, drop = FALSE], na.rm = TRUE)
  list(M = M, dropout = dropout, trend = trend, tau = tau, s_beta = max(mad(lfc), 0.1),
       mu_loc = median(M, na.rm = TRUE), mu_scale = sd(M, na.rm = TRUE))
}
stan_data <- function(h, prot) {   # data of one protein for the Stan model
  y <- h$M[prot, ]; g <- as.integer(factor(grp[names(y)], levels = c(REFERENCE, TEST))); mis <- is.na(y)
  pep <- if (is.na(pep_n[prot])) log(MIN_PEPTIDES) else log(pep_n[prot])
  list(N_obs = sum(!mis), y = as.array(unname(y[!mis])), g_obs = as.array(g[!mis]),
       N_mis = sum(mis), g_mis = as.array(g[mis]), rho = as.array(unname(h$dropout[names(y)[mis], "rho"])),
       zeta = as.array(unname(h$dropout[names(y)[mis], "zeta"])), mu_loc = h$mu_loc, mu_scale = h$mu_scale, s_beta = h$s_beta,
       logsig_loc = unname(predict(h$trend, data.frame(mean = mean(y, na.rm = TRUE), Peptides = pep))), logsig_scale = h$tau)
}
bayes_fit <- function(mod, d, lfc, seed) {
  invisible(capture.output(fit <- suppressWarnings(rstan::sampling(mod, data = d, chains = 4, iter = 2000, seed = seed, refresh = 0,
                                                                    cores = 1, control = list(adapt_delta = 0.95)))))
  b <- rstan::extract(fit, "beta")$beta
  s <- rstan::summary(fit, pars = c("mu", "beta", "log_sigma"))$summary
  c(mean = mean(b), low = unname(quantile(b, 0.025)), high = unname(quantile(b, 0.975)),
    P_up = mean(b > lfc), P_down = mean(b < -lfc), Rhat = max(s[, "Rhat"]), ESS = min(s[, "n_eff"]),
    Divergent = rstan::get_num_divergent(fit))
}
Bm_obs <- dcast(B_prot, Protein ~ Run, value.var = "LogObserved")
Bm_obs <- as.matrix(Bm_obs, rownames = "Protein"); colnames(Bm_obs) <- lab[colnames(Bm_obs)]
hyp <- list(ttest = hyper(A_obs[, c(ref_cols, test_cols)]), MSstats = hyper(Bm_obs[, intersect(c(ref_cols, test_cols), colnames(Bm_obs))]))
prior_tab <- rbindlist(lapply(names(hyp), function(k) data.table(Pipeline = k, s_log2FC = round(hyp[[k]]$s_beta, 3),
  tau_log_sigma = round(hyp[[k]]$tau, 3), t(round(coef(hyp[[k]]$trend), 3)),
  Dropout_rho_median = round(median(hyp[[k]]$dropout[, "rho"]), 2), Dropout_zeta_median = round(median(hyp[[k]]$dropout[, "zeta"]), 2),
  Dropout_curves_fitted = paste0(sum(hyp[[k]]$dropout[, "fitted"]), "/", nrow(hyp[[k]]$dropout)),
  Intensity_mean = round(hyp[[k]]$mu_loc, 2), Intensity_sd = round(hyp[[k]]$mu_scale, 2))))
fwrite(prior_tab, file.path(D4, "4.2_Empirical_priors.csv"))
cat("Empirical priors (from all proteins):\n"); print(prior_tab, row.names = FALSE)
relevance <- list()
for (k in c("ttest", "MSstats")) {
  sig <- results[[k]][Significant == TRUE]
  # Results are stored per protein, for these samples and settings: only the proteins
  # not yet in the file are fitted (none, when the instructor's file is used).
  key <- paste(c("model-v5", key_b, sort(run_info$Run), LFC, SEED, MSSTATS_USE, NORMALISATION, DELTA, CENSOR_TRIM, DRIFT, MIN_VALID), collapse = ";")
  cache <- file.path(D4, paste0("Bayesian_", k, "_", key_hash(key), ".csv"))   # one file per configuration
  old <- file.path(D4, paste0("Bayesian_", k, ".csv"))   # name used by earlier versions of the script
  if (!file.exists(cache) && file.exists(old) && file.exists(paste0(old, ".key")) && identical(readLines(paste0(old, ".key")), key)) {
    file.rename(old, cache); file.rename(paste0(old, ".key"), paste0(cache, ".key"))
  }
  bay <- data.table(Protein = character(), Bayes_log2FC = numeric(), CrI_low = numeric(), CrI_high = numeric(),
                    P_up = numeric(), P_down = numeric(), Rhat = numeric(), ESS = numeric(), Divergent = numeric())
  if (file.exists(cache) && file.exists(paste0(cache, ".key")) && identical(readLines(paste0(cache, ".key")), key)) bay <- fread(cache)
  todo <- setdiff(intersect(sig$Protein, rownames(hyp[[k]]$M)), bay$Protein)
  if (length(todo) == 0) {
    cat(k, ": Bayesian results found for all", nrow(sig), "significant proteins: reading them.\n")
  } else if (!pkgbuild::has_build_tools(debug = FALSE)) {
    # Without a compiler the missing proteins stay without a Bayesian result (the rest of the
    # script still runs). To fit them: install Rtools (Windows), restart RStudio, run again.
    cat(k, ":", length(todo), "significant proteins have no Bayesian result in the file provided and no C++ compiler",
        "was found to fit them (Windows: install Rtools, guide 3.1). They are left without a Bayesian probability.\n")
  } else {
    cat(k, ": compiling the model (only the first time, 1-2 min)...\n")
    invisible(capture.output(mod <- rstan::stan_model(file = stan_file, auto_write = TRUE)))
    n_cores <- if (CORES > 0) CORES else max(1, future::availableCores() - 1)
    cat(k, ": fitting", length(todo), "Bayesian models on", n_cores, "cores...\n")
    t0 <- Sys.time()
    inputs <- lapply(todo, function(prot) stan_data(hyp[[k]], prot))
    chunks <- split(seq_along(inputs), (seq_along(inputs) - 1) %% min(n_cores, length(inputs)))
    future::plan(future::multisession, workers = n_cores)
    fits <- future.apply::future_lapply(chunks, function(idx) {
      rstan::rstan_options(auto_write = TRUE)
      m <- rstan::stan_model(file = stan_file, auto_write = TRUE)   # reads the compiled model
      lapply(idx, function(i) tryCatch(bayes_fit(m, inputs[[i]], LFC, SEED + i), error = function(e) rep(NA_real_, 8)))
    }, future.seed = TRUE, future.packages = "rstan")
    future::plan(future::sequential)
    new_fits <- data.table(todo[unlist(chunks)], do.call(rbind, unlist(fits, recursive = FALSE)))
    setnames(new_fits, names(bay))
    bay <- rbind(bay, new_fits)
    fwrite(bay, cache); writeLines(key, paste0(cache, ".key"))
    cat("  done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  }
  r <- merge(sig, bay, by = "Protein", all.x = TRUE)
  r[, Missing_values := rowSums(is.na(hyp[[k]]$M[match(Protein, rownames(hyp[[k]]$M)), , drop = FALSE]))]
  r[, FC_relevant := abs(log2FC) >= LFC]
  r[, Bayes_relevant := !is.na(P_up) & (P_up >= PROB | P_down >= PROB)]
  relevance[[k]] <- r[order(-pmax(P_up, P_down, na.rm = TRUE))]
  fwrite(relevance[[k]], file.path(D4, paste0("4.2_", k, "_significant_with_relevance.csv")))
  bad <- r[!is.na(Rhat) & (Rhat > 1.01 | ESS < 400 | Divergent > 0)]
  cat(k, ": significant", nrow(r), "| relevant by FC", sum(r$FC_relevant), "| relevant by Bayes", sum(r$Bayes_relevant),
      "| models with R-hat > 1.01, ESS < 400 or divergences:", nrow(bad), "| absent in one condition (infinite log2FC):", sum(!is.finite(r$log2FC)), "\n")
  r[, Criteria := ifelse(FC_relevant & Bayes_relevant, "Both", ifelse(FC_relevant, "Fold change only", ifelse(Bayes_relevant, "Bayesian only", "Neither")))]
  p <- ggplot(r[!is.na(P_up) & is.finite(log2FC)], aes(log2FC, pmax(P_up, P_down), colour = Criteria)) + geom_point(size = 1.2, alpha = 0.8) +
    geom_hline(yintercept = PROB, linetype = "dashed") + geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed") +
    scale_colour_manual(values = c("Both" = "#14213D", "Fold change only" = "#B5541B", "Bayesian only" = "#2A6F97", "Neither" = "grey70")) +
    theme_bw() + labs(title = paste(k, ": two criteria of relevance (significant proteins)"),
                      x = "log2FC", y = paste0("Bayesian probability of |log2FC| > ", LFC, " in that direction"))
  save_plot(p, D4, paste0("4.2_", k, "_relevance.png"), 7.5, 5.5)
}

# 4.3 The four lists of relevant proteins (ttest/MSstats x FC/Bayes)
lists <- list()
for (k in c("ttest", "MSstats")) for (crit in c("FC", "Bayes")) {
  r <- relevance[[k]][get(paste0(crit, "_relevant")) == TRUE]
  lists[[paste0(k, "_", crit)]] <- r
  fwrite(r, file.path(D4, paste0("4.3_Relevant_", k, "_", crit, ".csv")))
}
rel_summary <- rbindlist(lapply(names(lists), function(n) data.table(List = n, Up = sum(lists[[n]]$Direction == "up"),
                                                                    Down = sum(lists[[n]]$Direction == "down"))))
fwrite(rel_summary, file.path(D4, "4.3_Relevance_summary.csv"))
print(rel_summary, row.names = FALSE)

# 4.4 Comparison of the four lists
ids_list <- lapply(lists, function(r) r$Protein)
png(file.path(D4, "4.4_UpSet_relevant_proteins.png"), width = 8, height = 5, units = "in", res = 300)
if (sum(lengths(ids_list) > 0) >= 2) print(UpSetR::upset(UpSetR::fromList(ids_list[lengths(ids_list) > 0]), nsets = 4, order.by = "freq"))
invisible(dev.off())
all_ids <- unique(unlist(ids_list))
membership <- data.table(Protein = all_ids, Description = describe(all_ids))
for (n in names(ids_list)) membership[, (n) := Protein %in% ids_list[[n]]]
membership[, N_lists := rowSums(.SD), .SDcols = names(ids_list)]
membership <- merge(membership, results$ttest[, .(Protein, log2FC_ttest = round(log2FC, 2))], by = "Protein", all.x = TRUE)
membership <- merge(membership, results$MSstats[, .(Protein, log2FC_MSstats = round(log2FC, 2))], by = "Protein", all.x = TRUE)
setorder(membership, -N_lists, -log2FC_MSstats)
fwrite(membership, file.path(D4, "4.4_Relevant_proteins_membership.csv"))
cat("Relevant proteins in all four lists:", sum(membership$N_lists == 4), "\n")


# Section 5. Functional enrichment (ORA and GSEA) ---------------------------------
# ORA on the up- and down-regulated proteins of each of the four lists; GSEA on all the
# proteins of each pipeline. File names say where the list comes from.

D5 <- out_dir("5_Enrichment")
unlink(file.path(D5, c("ORA_plots", "GSEA_plots")), recursive = TRUE)   # no plots left from earlier runs
acc <- function(x) sub("^[^|]*\\|([^|]+)\\|.*$", "\\1", x)       # tr|W6QSR7|W6QSR7_ECTO5 -> W6QSR7
ipr <- fread(file.path(DATA_DIR, INTERPRO), header = FALSE, sep = "\t", quote = "", fill = TRUE,
             col.names = c("Protein", "MD5", "Length", "Analysis", "Accession", "Description", "Start", "End",
                           "Evalue", "Status", "Date", "IPR", "IPR_desc", "GO", "Pathways"))
ipr[, Protein := acc(Protein)]
sets <- list()
for (a in c("Pfam", "PANTHER", "CDD", "SMART", "NCBIfam")) {
  x <- unique(ipr[Analysis == a, .(Term = Accession, Protein, Name = ifelse(Description %in% c("-", ""), Accession, Description))])
  sets[[a]] <- x
}
sets$InterPro <- unique(ipr[IPR != "-" & IPR != "", .(Term = IPR, Protein, Name = IPR_desc)])
go <- unique(ipr[GO != "-" & GO != "", .(Protein, GO)])[, .(Term = unique(unlist(regmatches(GO, gregexpr("GO:[0-9]{7}", GO))))), by = Protein]
go <- unique(go)
go_info <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(go$Term), columns = c("ONTOLOGY", "TERM"), keytype = "GOID"))
go[, Ontology := go_info$ONTOLOGY[match(Term, go_info$GOID)]][, Name := go_info$TERM[match(Term, go_info$GOID)]]
for (o in c("BP", "MF", "CC")) sets[[paste0("GO_", o)]] <- go[Ontology == o, .(Term, Protein, Name)]
pw <- unique(ipr[Pathways != "-" & Pathways != "", .(Protein, Pathways)])[, .(Term = unique(unlist(strsplit(Pathways, "|", fixed = TRUE)))), by = Protein]
sets$MetaCyc  <- unique(pw[grepl("^MetaCyc:", Term), .(Term = sub("^MetaCyc:", "", Term), Protein)])[, Name := Term]
# Reactome: InterProScan lists each pathway once per species (human, mouse, fly...); keep the
# human version only, otherwise the same pathway is tested ~15 times
sets$Reactome <- unique(pw[grepl("^Reactome:R-HSA-", Term), .(Term = sub("^Reactome:", "", Term), Protein)])[, Name := Term]
ko <- fread(file.path(DATA_DIR, GHOSTKOALA), header = FALSE, sep = "\t", fill = TRUE, col.names = c("Protein", "KO"))
ko <- ko[KO != "", .(Protein = acc(Protein), KO)]
if (!file.exists(file.path(DATA_DIR, KEGG_MAP))) {   # build it once from the KEGG REST service
  km <- tryCatch({
    link <- fread("https://rest.kegg.jp/link/pathway/ko", header = FALSE, col.names = c("KO", "Pathway"))
    link <- link[grepl("^path:map", Pathway), .(KO = sub("^ko:", "", KO), Pathway = sub("^path:", "", Pathway))]
    nm <- fread("https://rest.kegg.jp/list/pathway", header = FALSE, sep = "\t", col.names = c("Pathway", "Name"))
    merge(link, nm[, .(Pathway = sub("^path:", "", Pathway), Name)], by = "Pathway")[, .(KO, Pathway, Name)]
  }, error = function(e) NULL)
  if (!is.null(km)) { fwrite(km, file.path(DATA_DIR, KEGG_MAP), sep = "\t"); cat("KEGG pathways downloaded and saved to", KEGG_MAP, "\n") }
}
if (file.exists(file.path(DATA_DIR, KEGG_MAP))) {
  km <- fread(file.path(DATA_DIR, KEGG_MAP))   # columns KO, Pathway, Name
  sets$KEGG <- unique(merge(ko, km, by = "KO", allow.cartesian = TRUE)[, .(Term = Pathway, Protein, Name)])
} else {
  cat("KEGG: file", KEGG_MAP, "not found and KEGG could not be reached; KEGG enrichment skipped.\n")
}
sets <- sets[intersect(names(ENRICH)[ENRICH], names(sets))]
# Terms that contain (almost) the same proteins are one result, not many: InterPro transfers
# MetaCyc and Reactome pathways from protein families, so one generic domain (e.g. a P-loop
# ATPase) can carry hundreds of "pathways". Terms sharing >= 80% of their proteins (Jaccard)
# are merged into the largest one, whose name says how many it stands for.
collapse_sets <- function(s, min_jaccard = 0.8) {
  s <- unique(s[, .(Term, Protein, Name)])
  tm <- s[, .N, by = Term][order(-N)]
  M <- Matrix::sparseMatrix(i = match(s$Term, tm$Term), j = match(s$Protein, unique(s$Protein)), x = 1)
  inter <- as.matrix(Matrix::tcrossprod(M)); n <- tm$N
  J <- inter / (outer(n, n, "+") - inter)
  rep_of <- rep(NA_integer_, length(n))
  for (t in seq_along(n)) if (is.na(rep_of[t])) { grp <- which(is.na(rep_of) & J[t, ] >= min_jaccard); rep_of[grp] <- t }
  k <- tabulate(rep_of, length(n))
  nm <- unique(s[, .(Term, Name)])[match(tm$Term, Term), Name]
  keep <- data.table(Term = tm$Term[k > 0], Name = ifelse(k[k > 0] > 1, paste0(nm[k > 0], " (+", k[k > 0] - 1, " similar terms)"), nm[k > 0]))
  merge(s[Term %in% keep$Term, .(Term, Protein)], keep, by = "Term")
}
sets <- lapply(sets, collapse_sets)
cat("Annotation sets:\n"); print(sapply(sets, function(s) c(Terms = uniqueN(s$Term), Proteins = uniqueN(s$Protein))))

# 5.1 How many relevant proteins have an annotation?
coverage <- rbindlist(lapply(names(lists), function(n) {
  pr <- lists[[n]]$Protein
  data.table(List = n, Proteins = length(pr), t(sapply(sets, function(s) sum(pr %in% s$Protein))),
             Any = sum(pr %in% unique(unlist(lapply(sets, `[[`, "Protein")))))
}))
fwrite(coverage, file.path(D5, "5.1_Annotation_coverage.csv"))
print(coverage, row.names = FALSE)

# 5.2 ORA: up- and down-regulated proteins of every list vs. all the proteins tested in its pipeline
universe <- list(ttest = results$ttest[!is.na(adj.pvalue), Protein], MSstats = results$MSstats[!is.na(adj.pvalue), Protein])
ora_all <- list()
for (n in names(lists)) for (dirn in c("up", "down")) for (a in names(sets)) {
  genes <- lists[[n]][Direction == dirn, Protein]
  if (length(genes) < 3) next
  e <- tryCatch(suppressMessages(clusterProfiler::enricher(genes, TERM2GENE = sets[[a]][, .(Term, Protein)], TERM2NAME = unique(sets[[a]][, .(Term, Name)]),
                                          universe = universe[[sub("_.*", "", n)]], pvalueCutoff = FDR, qvalueCutoff = 0.2,
                                          minGSSize = 3, maxGSSize = 500)), error = function(e) NULL)
  if (is.null(e)) next
  e@result <- e@result[e@result$Count >= ORA_MIN, ]   # a term found with 1-2 proteins is not an enrichment
  if (nrow(as.data.frame(e)) == 0) next
  ora_all[[paste(n, dirn, a)]] <- data.table(List = n, Direction = dirn, Annotation = a, as.data.frame(e))
  if (nrow(as.data.frame(e)) >= 2) save_plot(enrichplot::dotplot(e, showCategory = 15) + ggtitle(paste("ORA", n, dirn, a)),
                                             out_dir("5_Enrichment/ORA_plots"), paste0("ORA_", n, "_", dirn, "_", a, ".png"), 8, 6)
}
ora_all <- rbindlist(ora_all, fill = TRUE)
if (nrow(ora_all)) fwrite(ora_all, file.path(D5, "5.2_ORA_results.csv"))
for (n in names(lists)) if (nrow(ora_all) && any(ora_all$List == n)) fwrite(ora_all[List == n], file.path(D5, paste0("5.2_ORA_", n, ".csv")))

# 5.3 GSEA: all the proteins of each pipeline ranked by a signed z-score, sign(log2FC) x the normal
# quantile of the p-value: the same scale for both pipelines, and robust to proteins whose t is
# huge only because they have very few degrees of freedom
gsea_all <- list()
for (k in c("ttest", "MSstats")) {
  r <- results[[k]][is.finite(log2FC) & !is.na(pvalue)]
  r[, z := sign(log2FC) * qnorm(pmax(pvalue, 1e-300) / 2, lower.tail = FALSE)]
  ranks <- sort(setNames(r$z, r$Protein), decreasing = TRUE)
  fwrite(data.table(Protein = names(ranks), Statistic = ranks), file.path(D5, paste0("5.3_", k, "_ranked.rnk")), sep = "\t", col.names = FALSE)
  for (a in names(sets)) {
    set.seed(SEED)
    g <- tryCatch(suppressMessages(suppressWarnings(clusterProfiler::GSEA(ranks, TERM2GENE = sets[[a]][, .(Term, Protein)], TERM2NAME = unique(sets[[a]][, .(Term, Name)]),
                                        pvalueCutoff = FDR, minGSSize = 5, maxGSSize = 500, seed = TRUE, verbose = FALSE))), error = function(e) NULL)
    if (is.null(g) || nrow(as.data.frame(g)) == 0) next
    gsea_all[[paste(k, a)]] <- data.table(Pipeline = k, Annotation = a, as.data.frame(g))
    if (nrow(as.data.frame(g)) >= 2) save_plot(enrichplot::dotplot(g, showCategory = 15, split = ".sign") + facet_grid(. ~ .sign) + ggtitle(paste("GSEA", k, a)),
                                               out_dir("5_Enrichment/GSEA_plots"), paste0("GSEA_", k, "_", a, ".png"), 10, 6)
  }
}
gsea_all <- rbindlist(gsea_all, fill = TRUE)
if (nrow(gsea_all)) fwrite(gsea_all, file.path(D5, "5.3_GSEA_results.csv"))
for (k in c("ttest", "MSstats")) if (nrow(gsea_all) && any(gsea_all$Pipeline == k)) fwrite(gsea_all[Pipeline == k], file.path(D5, paste0("5.3_GSEA_", k, ".csv")))

# 5.4 How many independent findings? Significant terms explained by (almost) the same proteins
# are grouped (greedy, Jaccard >= 0.5 of the proteins behind each term, best term first)
cluster_terms <- function(d, genes) {
  if (!nrow(d)) return(d)
  d <- d[order(p.adjust)]; g <- strsplit(d[[genes]], "/", fixed = TRUE); cl <- rep(NA_integer_, nrow(d)); k <- 0
  for (i in seq_len(nrow(d))) if (is.na(cl[i])) { k <- k + 1
    j <- which(is.na(cl) & sapply(g, function(x) length(intersect(x, g[[i]])) / length(union(x, g[[i]]))) >= 0.5); cl[j] <- k }
  d[, Cluster := cl]
  d[, .(Best_term = Description[1], Annotation = Annotation[1], p.adjust = signif(p.adjust[1], 2), Terms = .N,
        Proteins = uniqueN(unlist(strsplit(get(genes), "/"))),
        Main_proteins = paste(head(names(sort(table(unlist(strsplit(get(genes), "/"))), decreasing = TRUE)), 5), collapse = ", ")), by = Cluster]
}
if (nrow(ora_all)) {
  oc <- ora_all[, cluster_terms(.SD, "geneID"), by = .(List, Direction)]
  oc[, Main_proteins := sapply(strsplit(Main_proteins, ", "), function(v) paste(paste0(v, " (", substr(describe(v), 1, 30), ")"), collapse = "; "))]
  fwrite(oc, file.path(D5, "5.4_ORA_independent_findings.csv"))
}
if (nrow(gsea_all)) {
  gc <- gsea_all[, cluster_terms(.SD, "core_enrichment"), by = .(Pipeline, Direction = ifelse(NES > 0, "up", "down"))]
  gc[, Main_proteins := sapply(strsplit(Main_proteins, ", "), function(v) paste(paste0(v, " (", substr(describe(v), 1, 30), ")"), collapse = "; "))]
  fwrite(gc, file.path(D5, "5.4_GSEA_independent_findings.csv"))
  cat("GSEA: significant terms and independent findings\n"); print(gc[, .(Terms = sum(Terms), Findings = .N), by = .(Pipeline, Direction)], row.names = FALSE)
}
if (nrow(ora_all)) fwrite(ora_all[, .(Terms = .N), by = .(List, Direction, Annotation)], file.path(D5, "5_ORA_summary.csv"))
if (nrow(gsea_all)) fwrite(gsea_all[, .(Terms = .N, Up = sum(NES > 0), Down = sum(NES < 0)), by = .(Pipeline, Annotation)], file.path(D5, "5_GSEA_summary.csv"))
cat("ORA: significant terms\n"); if (nrow(ora_all)) print(dcast(ora_all, Annotation ~ List + Direction, fun.aggregate = length, value.var = "ID"))
cat("GSEA: significant terms\n"); if (nrow(gsea_all)) print(dcast(gsea_all, Annotation ~ Pipeline, fun.aggregate = length, value.var = "ID"))


# Section 6. Session information ---------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\nDone. All results are in:", OUT, "\n")
