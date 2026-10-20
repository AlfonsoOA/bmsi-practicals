# =============================================================================
# PC2_RNAseq.R
# Systems Molecular Biology - Practical 2 (RNA-seq) - 2026-2027
#
# Differential expression with two transcriptomes (reference and de novo) and
# two statistical methods (DESeq2 and limma-voom), and comparison of the four
# pipelines.
#
# HOW TO USE
#   1. Check the CONFIGURATION block below. With the files produced by the
#      terminal part of the practical (section 5.10 of the guide), nothing
#      needs to be changed.
#   2. Run the script section by section in RStudio: place the cursor in a
#      section and press Ctrl+Alt+T (Code > Run Region > Run Code Section),
#      or run everything with Ctrl+Shift+Enter.
#   3. All results are written to the folder "results" inside DATA_DIR.
#
# Files expected in DATA_DIR (any number of samples):
#   samples.tsv              columns "sample" and "condition" (others ignored)
#   <quant folder>/<sample>/quant.sf   one folder per sample, per transcriptome
#   <transcript->gene table> two columns, no header: transcript, gene
#   gene_annotation.tsv      optional: gene, symbol, description, GO
# =============================================================================


# CONFIGURATION: the only part you may need to edit ---------------------------

# Folder with the input files. Empty = automatic:
#   Windows: C:/Users/<you>/BMSi/PC2_RNAseq_R   macOS/Linux: ~/BMSi/PC2_RNAseq_R
DATA_DIR <- ""

# Table of samples (tab-separated; columns "sample" and "condition")
SAMPLE_SHEET <- "samples.tsv"

# Transcriptomes to compare: a short name for each, its quantification folder
# and its transcript -> gene table. The FIRST one is the reference. Add or
# remove lines to compare more or fewer transcriptomes.
SOURCES <- list(
  Ensembl = c(quant = "quant_ensembl", tx2gene = "ensembl_tx2gene.tsv"),
  Trinity = c(quant = "quant_trinity", tx2gene = "trinity_tx2gene.tsv")
)

# Gene annotation (gene, symbol, description, GO). "" if there is none.
ANNOTATION <- "gene_annotation.tsv"

# Program that produced the quantification files: "salmon", "kallisto", "rsem"...
QUANT_TYPE <- "salmon"

# Conditions to compare: log2FC = log2(TEST / REFERENCE).
# Empty = first and second condition, in the order of the sample table.
REFERENCE_LEVEL <- ""
TEST_LEVEL      <- ""

# Thresholds that define a differentially expressed gene (DEG)
FDR <- 0.05   # adjusted p-value (Benjamini-Hochberg)
LFC <- 1      # absolute log2 fold change (1 = two-fold)

# FALSE: test log2FC = 0, then keep genes with |log2FC| >= LFC (usual practice).
# TRUE : test directly whether |log2FC| > LFC (DESeq2 lfcThreshold, limma treat).
TEST_AGAINST_LFC <- FALSE

# EXTRA: genes to remove before the differential expression analysis (for
# example rRNA genes found in section 1.6). A text file with one gene ID per
# line, in DATA_DIR. "" = remove nothing. Results then go to a separate folder
# (results_without_<file name>), so both analyses can be compared.
EXCLUDE_GENES <- ""

TOP_HEATMAP <- 50     # number of genes in the heatmaps
SEED        <- 2026   # makes the permutation tests reproducible


# Section 0. Packages ----------------------------------------------------------
# Installs the packages that are missing (only the first time) and loads them.

if (getRversion() < "4.3.0") {
  stop("This script needs R 4.3.0 or newer. You have R ", getRversion(),
       ". Install the current R from https://cran.r-project.org")
}
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c("ggplot2", "ggrepel", "pheatmap", "vegan", "UpSetR")
bioc_pkgs <- c("tximport", "DESeq2", "edgeR", "limma")

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

for (p in c(cran_pkgs, bioc_pkgs)) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
cat("Package check:\n")
for (p in c(bioc_pkgs, cran_pkgs)) cat(sprintf("  %-10s %s\n", p, packageVersion(p)))


# Section 1. Data import -------------------------------------------------------
# Reads the sample table and imports the quantification of every transcriptome,
# summarising transcripts into genes with tximport.

# 1.1 Folders
if (DATA_DIR == "") {
  home <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
  DATA_DIR <- file.path(home, "BMSi", "PC2_RNAseq_R")
}
if (!dir.exists(DATA_DIR)) {
  stop("The data folder does not exist: ", DATA_DIR,
       "\nRun section 5.10 of the guide, or write the right folder in DATA_DIR.")
}
OUT <- file.path(DATA_DIR, if (EXCLUDE_GENES == "") "results" else
  paste0("results_without_", tools::file_path_sans_ext(basename(EXCLUDE_GENES))))
dir.create(OUT, showWarnings = FALSE)
cat("Data folder:   ", DATA_DIR, "\nResults folder:", OUT, "\n")

# 1.2 Samples and conditions
samples <- read.delim(file.path(DATA_DIR, SAMPLE_SHEET), stringsAsFactors = FALSE)
if (!all(c("sample", "condition") %in% names(samples))) {
  stop(SAMPLE_SHEET, " must have the columns 'sample' and 'condition'.")
}
conds <- unique(samples$condition)
if (length(conds) < 2) stop("The sample table has only one condition.")
if (REFERENCE_LEVEL == "") REFERENCE_LEVEL <- conds[1]
if (TEST_LEVEL == "")      TEST_LEVEL      <- setdiff(conds, REFERENCE_LEVEL)[1]
samples <- samples[samples$condition %in% c(REFERENCE_LEVEL, TEST_LEVEL), ]
samples$condition <- factor(samples$condition, levels = c(REFERENCE_LEVEL, TEST_LEVEL))
rownames(samples) <- samples$sample
cat("\nComparison: ", TEST_LEVEL, " vs ", REFERENCE_LEVEL, " (reference)\n", sep = "")
print(table(samples$condition))

# 1.3 Gene annotation (optional)
annot <- NULL
if (ANNOTATION != "" && file.exists(file.path(DATA_DIR, ANNOTATION))) {
  annot <- read.delim(file.path(DATA_DIR, ANNOTATION), quote = "", stringsAsFactors = FALSE)
  annot <- annot[!duplicated(annot[[1]]), ]
  names(annot)[1] <- "gene"
  cat("Annotation:", nrow(annot), "genes\n")
} else if (ANNOTATION != "") {
  warning("Annotation file not found: ", file.path(DATA_DIR, ANNOTATION),
          ". The results will have no gene names. Check section 5.10 of the guide.", immediate. = TRUE)
}

# 1.4 Import every transcriptome
import_source <- function(src) {
  t2g <- read.delim(file.path(DATA_DIR, src[["tx2gene"]]), header = FALSE,
                    col.names = c("TXNAME", "GENEID"), quote = "", stringsAsFactors = FALSE)
  files <- file.path(DATA_DIR, src[["quant"]], samples$sample, "quant.sf")
  names(files) <- samples$sample
  if (!all(file.exists(files))) {
    stop("Quantification files not found:\n  ", paste(files[!file.exists(files)], collapse = "\n  "))
  }
  list(
    de = tximport(files, type = QUANT_TYPE, tx2gene = t2g),     # for DESeq2 (length offsets)
    lv = tximport(files, type = QUANT_TYPE, tx2gene = t2g,      # for limma-voom
                  countsFromAbundance = "lengthScaledTPM"),
    genes = unique(t2g$GENEID)
  )
}
txi <- lapply(SOURCES, import_source)
REF_SOURCE <- names(SOURCES)[1]

import_summary <- data.frame(
  transcriptome = names(txi),
  genes_in_table = sapply(txi, function(x) length(x$genes)),
  genes_with_reads = sapply(txi, function(x) sum(rowSums(x$de$counts) > 0))
)
print(import_summary, row.names = FALSE)

# 1.5 Reads in genes that are absent from the reference transcriptome
# For each other transcriptome: percentage of the reads of each sample that fall
# in genes without a counterpart in the reference.
novel <- do.call(rbind, lapply(setdiff(names(txi), REF_SOURCE), function(s) {
  cnt <- txi[[s]]$de$counts
  is_novel <- !(rownames(cnt) %in% txi[[REF_SOURCE]]$genes)
  data.frame(transcriptome = s, sample = colnames(cnt),
             condition = samples[colnames(cnt), "condition"],
             novel_genes = sum(is_novel),
             pct_reads_in_novel_genes = round(100 * colSums(cnt[is_novel, , drop = FALSE]) / colSums(cnt), 2))
}))
if (!is.null(novel)) {
  write.csv(novel, file.path(OUT, "Reads_in_genes_absent_from_reference.csv"), row.names = FALSE)
  print(novel, row.names = FALSE)
  p <- ggplot(novel, aes(sample, pct_reads_in_novel_genes, fill = condition)) +
    geom_col() + facet_wrap(~ transcriptome) + theme_bw() +
    labs(x = NULL, y = "% of reads in genes absent from the reference",
         title = paste("Reads in genes without a counterpart in", REF_SOURCE))
  ggsave(file.path(OUT, "Reads_in_genes_absent_from_reference.png"), p, width = 7, height = 4.5, dpi = 300)
}

# 1.6 Genes that take the most reads
# A handful of genes with a large share of the reads usually means rRNA or other
# contamination (or an annotation error), and it can distort the normalisation.
top_genes <- do.call(rbind, lapply(names(txi), function(s) {
  cnt <- txi[[s]]$de$counts
  pct <- 100 * sweep(cnt, 2, colSums(cnt), "/")
  ord <- head(order(rowMeans(pct), decreasing = TRUE), 20)
  out <- data.frame(transcriptome = s, rank = seq_along(ord), gene = rownames(cnt)[ord],
                    mean_pct_reads = round(rowMeans(pct)[ord], 2),
                    cumulative_pct = round(cumsum(rowMeans(pct)[ord]), 2),
                    round(pct[ord, , drop = FALSE], 2), check.names = FALSE)
  if (!is.null(annot)) {
    out <- cbind(out[, 1:3], annot[match(out$gene, annot$gene), intersect(c("symbol", "description"), names(annot)), drop = FALSE],
                 out[, -(1:3)])
  }
  out
}))
write.csv(top_genes, file.path(OUT, "Top_genes_by_reads.csv"), row.names = FALSE)
cat("\nShare of the reads taken by the most abundant genes:\n")
for (s in names(txi)) {
  tg <- top_genes[top_genes$transcriptome == s, ]
  cat(sprintf("  %-10s top 1 gene: %5.1f%%   top 10 genes: %5.1f%%   (most abundant: %s)\n",
              s, tg$cumulative_pct[1], tg$cumulative_pct[10], tg$gene[1]))
}

# 1.7 EXTRA: remove the genes listed in EXCLUDE_GENES
if (EXCLUDE_GENES != "") {
  excl <- unique(trimws(readLines(file.path(DATA_DIR, EXCLUDE_GENES))))
  excl <- excl[excl != ""]
  for (s in names(txi)) for (k in c("de", "lv")) {
    keep_rows <- !(rownames(txi[[s]][[k]]$counts) %in% excl)
    for (m in c("abundance", "counts", "length")) txi[[s]][[k]][[m]] <- txi[[s]][[k]][[m]][keep_rows, , drop = FALSE]
  }
  cat("\nRemoved before the analysis:", length(excl), "genes listed in", EXCLUDE_GENES, "\n")
  for (s in names(txi)) cat("  ", s, ":", sum(excl %in% txi[[s]]$genes), "of them present\n")
}

# Section 2. Functions shared by both methods -----------------------------------
# (run this section before sections 3 and 4)

annotate <- function(tab) {
  if (is.null(annot)) return(tab[order(tab$padj), ])
  out <- merge(tab, annot[, intersect(c("gene", "symbol", "description"), names(annot))],
               by = "gene", all.x = TRUE, sort = FALSE)
  out[order(out$padj), ]
}

is_deg <- function(tab) !is.na(tab$padj) & tab$padj < FDR & abs(tab$log2FC) >= LFC

min_permutation_p <- function(groups) {
  n <- table(groups)
  if (length(n) != 2) return(NA)
  splits <- choose(sum(n), n[1]) / ifelse(n[1] == n[2], 2, 1)
  1 / splits
}

save_pipeline_outputs <- function(name, tab, expr, label) {
  tab <- annotate(tab)
  deg <- is_deg(tab)
  write.csv(tab, file.path(OUT, paste0(name, "_all_genes.csv")), row.names = FALSE)
  write.csv(tab[deg, ], file.path(OUT, paste0(name, "_DEGs.csv")), row.names = FALSE)

  # Volcano plot
  tab$status <- ifelse(deg, "DEG",
                ifelse(!is.na(tab$padj) & tab$padj < FDR, "Significant, small change", "Not significant"))
  tab$status <- factor(tab$status, levels = c("DEG", "Significant, small change", "Not significant"))
  tab$name <- if ("symbol" %in% names(tab)) ifelse(is.na(tab$symbol) | tab$symbol == "", tab$gene, tab$symbol) else tab$gene
  top <- head(tab[deg, ], 10)
  p <- ggplot(tab[!is.na(tab$padj), ], aes(log2FC, -log10(padj), colour = status)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed") +
    geom_hline(yintercept = -log10(FDR), linetype = "dashed") +
    ggrepel::geom_text_repel(data = top, aes(label = name), colour = "black", size = 2.5, max.overlaps = 20) +
    scale_colour_manual(values = c("DEG" = "#B5541B", "Significant, small change" = "#2A6F97",
                                   "Not significant" = "grey70"), drop = FALSE) +
    theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(title = paste("Volcano plot -", name),
         subtitle = paste0("DEG: adjusted p < ", FDR, " and |log2FC| >= ", LFC),
         x = paste0("log2 fold change (", TEST_LEVEL, " / ", REFERENCE_LEVEL, ")"), y = "-log10 adjusted p-value")
  ggsave(file.path(OUT, paste0(name, "_Volcano.png")), p, width = 7, height = 6.5, dpi = 300)

  # PERMANOVA: are the samples of the two conditions separated?
  # (same distances that the PCA represents, but using all genes)
  set.seed(SEED)
  perm <- suppressMessages(vegan::adonis2(dist(t(expr)) ~ condition,
                                          data = samples[colnames(expr), , drop = FALSE], permutations = 999))
  min_p <- min_permutation_p(samples$condition)
  perm_txt <- paste0("PERMANOVA (all genes): R2 = ", round(perm$R2[1], 2), ", p = ", signif(perm$`Pr(>F)`[1], 2),
                     " (smallest possible: ", signif(min_p, 2), ")")

  # Histogram of the raw p-values: a flat background (true null genes) plus a
  # peak near 0 (changed genes) is what a well-behaved test gives
  ph <- ggplot(tab[!is.na(tab$pvalue), ], aes(pvalue)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "#2A6F97", colour = "white") +
    theme_bw() + labs(title = paste("p-value histogram -", name),
                      subtitle = if (TEST_AGAINST_LFC) "Tested against |log2FC| > LFC: most genes pile up near p = 1" else
                        "Expected: flat background plus a peak at the left",
                      x = "p-value (before adjustment)", y = "Number of genes")
  ggsave(file.path(OUT, paste0(name, "_pvalue_histogram.png")), ph, width = 6, height = 4.5, dpi = 300)

  # Every gene ranked by its test statistic (sign = direction), with no
  # threshold: input for rank-based enrichment analysis (GSEA, Practical 4)
  rnk <- tab[!is.na(tab$stat), c("gene", "stat")]
  rnk <- rnk[order(rnk$stat, decreasing = TRUE), ]
  write.table(rnk, file.path(OUT, paste0(name, "_ranked.rnk")), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # PCA of the 500 most variable genes
  top_var <- head(order(apply(expr, 1, var), decreasing = TRUE), 500)
  pca <- prcomp(t(expr[top_var, ]))
  ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pdat <- data.frame(pca$x[, 1:2], samples[colnames(expr), ])
  p <- ggplot(pdat, aes(PC1, PC2, colour = condition, label = sample)) +
    geom_point(size = 3) + ggrepel::geom_text_repel(show.legend = FALSE) + theme_bw() +
    labs(title = paste("PCA -", name, "-", label), subtitle = perm_txt,
         x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"))
  ggsave(file.path(OUT, paste0(name, "_PCA.png")), p, width = 6, height = 5, dpi = 300)

  # Heatmap of the most significant DEGs
  hm_genes <- head(tab$gene[deg], TOP_HEATMAP)
  if (length(hm_genes) >= 2) {
    pheatmap(expr[hm_genes, ], scale = "row", show_rownames = length(hm_genes) <= 60,
             labels_row = tab$name[match(hm_genes, tab$gene)],
             annotation_col = samples[colnames(expr), "condition", drop = FALSE],
             main = paste0(name, ": top ", length(hm_genes), " DEGs (", label, ", row z-scores)"),
             filename = file.path(OUT, paste0(name, "_Heatmap.png")), width = 7, height = 9)
  }


  report <- c(
    paste("Pipeline:", name),
    paste("Comparison:", TEST_LEVEL, "vs", REFERENCE_LEVEL, "(reference)"),
    paste("Samples:", paste(names(table(samples$condition)), table(samples$condition), sep = " = ", collapse = ", ")),
    paste("Genes analysed (after filtering low counts):", nrow(tab)),
    paste0("DEGs (adj. p < ", FDR, " and |log2FC| >= ", LFC,
           if (TEST_AGAINST_LFC) ", tested against the threshold" else "", "): ", sum(deg)),
    paste("  up in", TEST_LEVEL, ":", sum(deg & tab$log2FC > 0)),
    paste("  down in", TEST_LEVEL, ":", sum(deg & tab$log2FC < 0)),
    paste("Significant but |log2FC| below the threshold:", sum(!is.na(tab$padj) & tab$padj < FDR & !deg)),
    "",
    paste0("PERMANOVA on ", label, " (Euclidean distances, all genes):"),
    paste0("  R2 = ", round(perm$R2[1], 3), ", p = ", signif(perm$`Pr(>F)`[1], 3)),
    paste0("  Smallest p-value possible with this design: ", signif(min_p, 3),
           " (only ", round(1 / min_p), " different ways to split the samples into the two groups)")
  )
  writeLines(report, file.path(OUT, paste0(name, "_Report.txt")))
  cat(paste(report, collapse = "\n"), "\n\n")
  invisible(tab)
}

results <- list()   # the four pipelines are stored here


# Section 3. Differential expression with DESeq2 -------------------------------
# Negative binomial model on the counts; vst-transformed values for the plots.

for (s in names(txi)) {
  name <- paste0(s, "_DESeq2")
  cat("=====", name, "=====\n")
  dds <- DESeqDataSetFromTximport(txi[[s]]$de, colData = samples, design = ~ condition)
  keep <- rowSums(counts(dds) >= 10) >= min(table(samples$condition))
  dds <- DESeq(dds[keep, ], quiet = TRUE)
  res <- results(dds, contrast = c("condition", TEST_LEVEL, REFERENCE_LEVEL), alpha = FDR,
                 lfcThreshold = if (TEST_AGAINST_LFC) LFC else 0)
  tab <- data.frame(gene = rownames(res), mean_expression = res$baseMean,
                    log2FC = res$log2FoldChange, stat = res$stat, pvalue = res$pvalue, padj = res$padj)
  expr <- assay(vst(dds, blind = TRUE))
  results[[name]] <- save_pipeline_outputs(name, tab, expr, "vst values")
}


# Section 4. Differential expression with limma-voom ---------------------------
# TMM normalisation, voom weights and moderated t-tests on log-CPM values.

for (s in names(txi)) {
  name <- paste0(s, "_limma")
  cat("=====", name, "=====\n")
  y <- DGEList(txi[[s]]$lv$counts)
  design <- model.matrix(~ condition, data = samples)
  y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y, design)
  fit <- lmFit(v, design)
  if (TEST_AGAINST_LFC) {
    tt <- topTreat(treat(fit, lfc = LFC), coef = 2, number = Inf, sort.by = "none")
  } else {
    tt <- topTable(eBayes(fit), coef = 2, number = Inf, sort.by = "none")
  }
  tab <- data.frame(gene = rownames(tt), mean_expression = tt$AveExpr,
                    log2FC = tt$logFC, stat = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val)
  results[[name]] <- save_pipeline_outputs(name, tab, v$E, "voom log-CPM")
}


# Section 5. Comparison of the pipelines ---------------------------------------
# Which choice changes the results more: the transcriptome or the method?

pipes <- names(results)
part <- cbind(source = sub("_[^_]+$", "", pipes), method = sub(".*_", "", pipes))
rownames(part) <- pipes
degs <- lapply(results, function(t) t$gene[is_deg(t)])

lfc_pair <- function(a, b) {
  m <- merge(results[[a]][, c("gene", "log2FC")], results[[b]][, c("gene", "log2FC")], by = "gene")
  m[complete.cases(m), ]
}

# 5.1 Every pair of pipelines that differ in only one choice
pairs <- list()
for (i in seq_along(pipes)) for (j in seq_along(pipes)) if (i < j) {
  same <- part[i, ] == part[j, ]
  if (sum(same) == 1) pairs[[length(pairs) + 1]] <- c(pipes[i], pipes[j], ifelse(same["source"], "method", "transcriptome"))
}
comparison <- do.call(rbind, lapply(pairs, function(pr) {
  a <- pr[1]; b <- pr[2]; m <- lfc_pair(a, b)
  shared <- length(intersect(degs[[a]], degs[[b]]))
  data.frame(changed = pr[3], pipeline_A = a, pipeline_B = b,
             DEGs_A = length(degs[[a]]), DEGs_B = length(degs[[b]]), DEGs_shared = shared,
             Jaccard = round(shared / length(union(degs[[a]], degs[[b]])), 3),
             genes_in_both = nrow(m),
             pearson_log2FC = round(cor(m$log2FC.x, m$log2FC.y), 3),
             spearman_log2FC = round(cor(m$log2FC.x, m$log2FC.y, method = "spearman"), 3))
}))
write.csv(comparison, file.path(OUT, "Comparison_table.csv"), row.names = FALSE)
print(comparison, row.names = FALSE)

# 5.2 Steiger's test: for each pipeline, does changing the method or changing
# the transcriptome alter its log2FC more? The two correlations share that
# pipeline, so they are dependent (Meng, Rosenthal & Rubin 1992).
steiger <- function(r12, r13, r23, n) {
  z1 <- atanh(r12); z2 <- atanh(r13)
  rm2 <- (r12^2 + r13^2) / 2
  f <- min(1, (1 - r23) / (2 * (1 - rm2)))
  h <- (1 - f * rm2) / (1 - rm2)
  z <- (z1 - z2) * sqrt((n - 3) / (2 * (1 - r23) * h))
  c(z = z, p = 2 * pnorm(-abs(z)))
}
steiger_lines <- c("Steiger's test for two dependent correlations (log2FC of the genes shared by the three pipelines)", "")
steiger_table <- NULL
for (a in pipes) {
  m_partner <- pipes[part[, "source"] == part[a, "source"] & part[, "method"] != part[a, "method"]]
  t_partner <- pipes[part[, "method"] == part[a, "method"] & part[, "source"] != part[a, "source"]]
  for (m in m_partner) for (t in t_partner) {
    genes <- Reduce(intersect, lapply(c(a, m, t), function(p) results[[p]]$gene))
    get_lfc <- function(p) results[[p]]$log2FC[match(genes, results[[p]]$gene)]
    x1 <- get_lfc(a); x2 <- get_lfc(m); x3 <- get_lfc(t)
    ok <- complete.cases(x1, x2, x3)
    r12 <- cor(x1[ok], x2[ok]); r13 <- cor(x1[ok], x3[ok]); r23 <- cor(x2[ok], x3[ok])
    st <- steiger(r12, r13, r23, sum(ok))
    verdict <- if (st["p"] >= 0.05) "no significant difference between the two changes" else
      if (r12 > r13) "changing the transcriptome alters the log2FC more than changing the method" else
        "changing the method alters the log2FC more than changing the transcriptome"
    steiger_lines <- c(steiger_lines,
      paste0("Starting pipeline: ", a, "  (", sum(ok), " genes)"),
      paste0("  method changed        r(", a, ", ", m, ") = ", round(r12, 3)),
      paste0("  transcriptome changed r(", a, ", ", t, ") = ", round(r13, 3)),
      paste0("  between the two       r(", m, ", ", t, ") = ", round(r23, 3)),
      paste0("  Z = ", round(st["z"], 2), ", p ", if (st["p"] < 1e-16) "< 1e-16" else paste("=", signif(st["p"], 3))),
      paste0("  -> ", verdict), "")
    steiger_table <- rbind(steiger_table, data.frame(
      start = a, method_changed = m, transcriptome_changed = t, genes = sum(ok),
      r_method = round(r12, 3), r_transcriptome = round(r13, 3),
      R2_method = round(r12^2, 3), R2_transcriptome = round(r13^2, 3),
      Z = round(st["z"], 2), p = signif(st["p"], 3)))
  }
}
if (!is.null(steiger_table)) {
  writeLines(steiger_lines, file.path(OUT, "Comparison_Steiger.txt"))
  write.csv(steiger_table, file.path(OUT, "Comparison_Steiger.csv"), row.names = FALSE)
  cat(paste(steiger_lines, collapse = "\n"), "\n")
}

# 5.3 Overlap of the DEG lists (UpSet plot)
if (sum(lengths(degs) > 0) >= 2) {
  png(file.path(OUT, "Comparison_UpSet.png"), width = 9, height = 5.5, units = "in", res = 300)
  print(UpSetR::upset(UpSetR::fromList(degs[lengths(degs) > 0]), nsets = length(degs),
                      order.by = "freq", nintersects = 30))
  invisible(dev.off())
}

# 5.4 log2FC of the same genes in the pairs of pipelines
sc <- do.call(rbind, lapply(pairs, function(pr) {
  m <- lfc_pair(pr[1], pr[2])
  data.frame(pair = paste0(pr[2], " vs ", pr[1], "\n(", pr[3], " changed)"), x = m$log2FC.x, y = m$log2FC.y)
}))
p <- ggplot(sc, aes(x, y)) + geom_point(size = 0.4, alpha = 0.3) +
  geom_abline(colour = "#B5541B") + facet_wrap(~ pair) + theme_bw() +
  labs(x = "log2FC, pipeline A", y = "log2FC, pipeline B", title = "Same genes, two pipelines")
ggsave(file.path(OUT, "Comparison_log2FC.png"), p, width = 8, height = 7, dpi = 300)

# 5.5 DEGs found only in genes that are absent from the reference transcriptome
novel_degs <- do.call(rbind, lapply(pipes[part[, "source"] != REF_SOURCE], function(p) {
  t <- results[[p]][is_deg(results[[p]]), ]
  t <- t[!(t$gene %in% txi[[REF_SOURCE]]$genes), c("gene", "mean_expression", "log2FC", "padj")]
  if (nrow(t) == 0) return(NULL)
  data.frame(pipeline = p, t)
}))
if (!is.null(novel_degs)) {
  write.csv(novel_degs, file.path(OUT, "DEGs_absent_from_reference.csv"), row.names = FALSE)
  cat("\nDEGs in genes absent from the reference:\n")
  print(table(novel_degs$pipeline, ifelse(novel_degs$log2FC > 0, paste("up in", TEST_LEVEL), paste("down in", TEST_LEVEL))))
}


# Section 6. Session information -------------------------------------------------
# Versions of R and of every package, saved for reproducibility.

writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\nDone. All results are in:", OUT, "\n")
print(list.files(OUT))
