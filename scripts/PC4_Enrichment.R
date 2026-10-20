# =============================================================================
# PC4_Enrichment.R
# Systems Molecular Biology - Practical 4 (Functional enrichment) - 2026-2027
#
# The results of Practical 2 (RNA-seq, Chlamydomonas: 2 transcriptomes x 2 methods)
# and Practical 3 (proteomics, E. oleovorans: 2 pipelines x 2 relevance criteria)
# are analysed again with ONE annotation source for both organisms (eggNOG-mapper),
# with the same ORA and GSEA settings, and then compared:
#   - descriptively: shared terms (Jaccard), NES of one pipeline against another;
#   - differential enrichment: functions where two pipelines disagree systematically;
#   - resampling (leave one sample out): does the pipeline change the enrichment of a
#     function more than the choice of samples does?
# Lists are also exported for Cytoscape/ClueGO.
#
# HOW TO USE
#   1. Unzip ForStudents_PC4.zip in your BMSi folder (guide 4.1).
#   2. Run the script section by section in RStudio (cursor in a section, Ctrl+Alt+T).
#   3. All results are written to the folder "results" inside DATA_DIR.
#
# RESULTS (folder "results")
#   1_Data           what was read, genes/proteins per pipeline
#   2_Statistics     one table per pipeline: log2FC, p-value, z; lists of up/down genes
#   3_Annotation     annotation coverage per organism and source
#   4_Enrichment     ORA and GSEA of every pipeline, same settings
#   5_Comparison     Jaccard of significant terms, NES of pipeline pairs
#   6_Differential   GSEA of the difference between two pipelines
#   7_Resampling     leave-one-sample-out: is the difference larger than sampling noise?
#   8_ClueGO         gene lists ready for Cytoscape/ClueGO
# =============================================================================


# CONFIGURATION: the only part you may need to edit ---------------------------

# Folder with the input files. Empty = automatic:
#   Windows: C:/Users/<you>/BMSi/PC4_Enrichment   macOS/Linux: ~/BMSi/PC4_Enrichment
DATA_DIR <- ""

# Practical 2 (inside DATA_DIR/PC2): Salmon quantifications, as in Practical 2
PC2 <- list(dir = "PC2", samples = "samples.tsv", reference = "Control", test = "NaCl",
            sources = list(Ensembl = c(quant = "quant_ensembl", tx2gene = "ensembl_tx2gene.tsv"),
                           Trinity = c(quant = "quant_trinity", tx2gene = "trinity_tx2gene.tsv")))
# Practical 3 (inside DATA_DIR/PC3): result tables of PC3_Proteomics.R (final analysis)
PC3 <- list(dir = "PC3", reference = "WT", test = "Mutante",
            ttest_matrix = "2.6_ttest_final_matrix.csv", ttest_table = "2.9_ttest_all_proteins.csv",
            msstats_matrix = "3.4_MSstats_protein_matrix.csv", msstats_table = "3.5_MSstats_all_proteins.csv",
            lists = c(ttest_FC = "4.3_Relevant_ttest_FC.csv", ttest_Bayes = "4.3_Relevant_ttest_Bayes.csv",
                      MSstats_FC = "4.3_Relevant_MSstats_FC.csv", MSstats_Bayes = "4.3_Relevant_MSstats_Bayes.csv"))

# Annotation (inside DATA_DIR): eggNOG-mapper output for each proteome, and what links the
# Chlamydomonas proteins to its genes (the protein FASTA of Ensembl Plants, headers with gene:)
EGGNOG_CHLAMY <- "chlamy.emapper.annotations"
EGGNOG_EOLEO  <- "eoleovorans.emapper.annotations"
CHLAMY_FASTA  <- "Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.pep.all.fa"
KEGG_MAP      <- "6_KEGG_KO_pathways.tsv"    # KEGG pathway names (from Practical 3)
ANNOTATIONS   <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Pfam", "COG")   # COG only if eggNOG gives category letters

# Thresholds (the same for every pipeline)
FDR      <- 0.05   # adjusted p-value (genes and terms)
LFC      <- 1      # |log2FC| for a gene of Practical 2 to be in the up/down list
MIN_SIZE <- 5      # smallest gene set tested (GSEA); ORA also needs >= ORA_MIN listed genes
MAX_SIZE <- 500
ORA_MIN  <- 3
DELTA_NES <- 1     # smallest NES difference between two pipelines that counts as relevant (section 7)

# Pairs of pipelines compared (differential enrichment and resampling)
PAIRS <- list(
  "PC2: DESeq2 vs limma (Ensembl)"   = c("Ensembl_DESeq2", "Ensembl_limma"),
  "PC2: DESeq2 vs limma (Trinity)"   = c("Trinity_DESeq2", "Trinity_limma"),
  "PC2: Ensembl vs Trinity (DESeq2)" = c("Ensembl_DESeq2", "Trinity_DESeq2"),
  "PC2: Ensembl vs Trinity (limma)"  = c("Ensembl_limma", "Trinity_limma"),
  "PC3: t-test vs MSstats"           = c("ttest", "MSstats"))

CORES <- 0      # 0 = automatic (all but one)
SEED  <- 2026


# Section 0. Packages ----------------------------------------------------------

if (getRversion() < "4.3.0") stop("This script needs R 4.3.0 or newer. You have R ", getRversion(), ".")
options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("data.table", "ggplot2", "ggrepel", "pheatmap", "UpSetR", "tximport", "DESeq2", "edgeR", "limma",
          "fgsea", "GO.db", "AnnotationDbi", "BiocParallel")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
}
for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("Package check:\n"); for (p in pkgs) cat(sprintf("  %-14s %s\n", p, packageVersion(p)))
set.seed(SEED)


# Section 1. Data import -------------------------------------------------------

if (DATA_DIR == "") {
  home <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
  DATA_DIR <- file.path(home, "BMSi", "PC4_Enrichment")
}
DATA_DIR <- normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE)
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
needed <- c(file.path(PC2$dir, PC2$samples), file.path(PC3$dir, c(PC3$ttest_matrix, PC3$ttest_table, PC3$msstats_matrix, PC3$msstats_table, PC3$lists)),
            EGGNOG_CHLAMY, EGGNOG_EOLEO, CHLAMY_FASTA)
if (!all(file.exists(file.path(DATA_DIR, needed))))
  stop("Put these files in ", DATA_DIR, " (guide 4.1) and run this section again:\n  ",
       paste(needed[!file.exists(file.path(DATA_DIR, needed))], collapse = "\n  "))
OUT <- file.path(DATA_DIR, "results")
out_dir <- function(name) { d <- file.path(OUT, name); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
save_plot <- function(p, folder, file, w = 7, h = 5.5) ggsave(file.path(folder, file), p, width = w, height = h, dpi = 300)
D1 <- out_dir("1_Data")
cat("Data folder:", DATA_DIR, "\n")

# 1.1 Practical 2: transcript quantifications summarised to genes (as in Practical 2)
samples2 <- read.delim(file.path(DATA_DIR, PC2$dir, PC2$samples))
samples2$condition <- factor(samples2$condition, levels = c(PC2$reference, PC2$test))
txi <- lapply(PC2$sources, function(s) {
  t2g <- read.delim(file.path(DATA_DIR, PC2$dir, s["tx2gene"]), header = FALSE)
  files <- setNames(file.path(DATA_DIR, PC2$dir, s["quant"], samples2$sample, "quant.sf"), samples2$sample)
  list(de = suppressMessages(tximport(files, type = "salmon", tx2gene = t2g)),
       lv = suppressMessages(tximport(files, type = "salmon", tx2gene = t2g, countsFromAbundance = "lengthScaledTPM")))
})
cat("Practical 2:", nrow(samples2), "samples;", paste(names(txi), sapply(txi, function(x) nrow(x$de$counts)), "genes", collapse = ", "), "\n")

# 1.2 Practical 3: protein values of each pipeline (final analysis of PC3_Proteomics.R)
read_matrix <- function(f) { m <- fread(file.path(DATA_DIR, PC3$dir, f)); x <- as.matrix(m[, -(1:2)]); rownames(x) <- m$Protein; x }
prot <- list(ttest = read_matrix(PC3$ttest_matrix), MSstats = read_matrix(PC3$msstats_matrix))
cond3 <- function(cols) factor(ifelse(startsWith(cols, PC3$test), PC3$test, PC3$reference), levels = c(PC3$reference, PC3$test))
cat("Practical 3:", paste(names(prot), sapply(prot, nrow), "proteins x", sapply(prot, ncol), "samples", collapse = "; "), "\n")
fwrite(data.table(Practical = c(rep("PC2", length(txi)), rep("PC3", length(prot))), Data = c(names(txi), names(prot)),
                  Features = c(sapply(txi, function(x) nrow(x$de$counts)), sapply(prot, nrow)),
                  Samples = c(rep(nrow(samples2), length(txi)), sapply(prot, ncol))), file.path(D1, "1_Inputs.csv"))


# Section 2. Differential statistics of every pipeline -------------------------
# Every pipeline gives, per gene or protein: log2FC, p-value, adjusted p-value and a
# signed z-score, z = sign(log2FC) x the normal quantile of p/2. The z-score puts all
# pipelines on the same scale; it ranks the genes for GSEA and is compared between pipelines.
# The functions below are reused, sample by sample, in the resampling (section 7).

D2 <- out_dir("2_Statistics")
zscore <- function(lfc, p) sign(lfc) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE)
de_pc2 <- function(source, method, keep_samples = samples2$sample) {
  sm <- samples2[samples2$sample %in% keep_samples, ]; sm$condition <- droplevels(sm$condition)
  x <- txi[[source]]
  if (method == "DESeq2") {
    t <- x$de; t$counts <- t$counts[, sm$sample]; t$abundance <- t$abundance[, sm$sample]; t$length <- t$length[, sm$sample]
    dds <- suppressMessages(DESeqDataSetFromTximport(t, colData = sm, design = ~ condition))
    dds <- suppressMessages(DESeq(dds[rowSums(counts(dds) >= 10) >= min(table(sm$condition)), ], quiet = TRUE))
    r <- results(dds, contrast = c("condition", PC2$test, PC2$reference))
    data.table(id = rownames(r), log2FC = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
  } else {
    design <- model.matrix(~ condition, sm)
    y <- DGEList(x$lv$counts[, sm$sample]); y <- calcNormFactors(y[filterByExpr(y, design), , keep.lib.sizes = FALSE])
    tt <- topTable(eBayes(lmFit(voom(y, design), design)), coef = 2, number = Inf, sort.by = "none")
    data.table(id = rownames(tt), log2FC = tt$logFC, pvalue = tt$P.Value, padj = tt$adj.P.Val)
  }
}
de_pc3 <- function(pipeline, drop = character(0)) {
  # t-test pipeline: Welch t-test on the final (imputed) matrix, as in PC3 2.9.
  # MSstats: two-group t-test with pooled variance on the MSstats protein values; for this
  # design (one run per biological sample) it is the linear model that MSstats fits.
  X <- prot[[pipeline]]; X <- X[, setdiff(colnames(X), drop), drop = FALSE]; g <- cond3(colnames(X))
  res <- t(apply(X, 1, function(y) {
    a <- y[g == PC3$test & is.finite(y)]; b <- y[g == PC3$reference & is.finite(y)]
    if (length(a) < 2 || length(b) < 2) return(c(NA, NA))
    r <- tryCatch(t.test(a, b, var.equal = pipeline == "MSstats"), error = function(e) NULL)
    if (is.null(r)) c(NA, NA) else c(mean(a) - mean(b), r$p.value)
  }))
  d <- data.table(id = rownames(X), log2FC = res[, 1], pvalue = res[, 2]); d[, padj := p.adjust(pvalue, "BH")]; d
}
stats <- list()
for (s in names(PC2$sources)) for (m in c("DESeq2", "limma")) stats[[paste0(s, "_", m)]] <- de_pc2(s, m)
for (k in names(prot)) stats[[k]] <- de_pc3(k)
for (n in names(stats)) {
  stats[[n]][, z := zscore(log2FC, pvalue)]
  fwrite(stats[[n]][order(pvalue)], file.path(D2, paste0("2_", n, "_statistics.csv")))
}

# 2.1 Lists of up- and down-regulated genes/proteins
#   PC2: adjusted p < FDR and |log2FC| >= LFC (as in Practical 2)
#   PC3: the four lists of relevant proteins of Practical 3 (t-test/MSstats x FC/Bayes)
lists <- list()
for (n in grep("^(Ensembl|Trinity)_", names(stats), value = TRUE)) {
  d <- stats[[n]][!is.na(padj) & padj < FDR & abs(log2FC) >= LFC]
  lists[[n]] <- list(up = d[log2FC > 0, id], down = d[log2FC < 0, id], universe = stats[[n]][!is.na(pvalue), id], pipeline = n)
}
for (n in names(PC3$lists)) {
  d <- fread(file.path(DATA_DIR, PC3$dir, PC3$lists[[n]]))
  k <- sub("_.*", "", n)
  lists[[n]] <- list(up = d[Direction == "up", Protein], down = d[Direction == "down", Protein],
                     universe = fread(file.path(DATA_DIR, PC3$dir, if (k == "ttest") PC3$ttest_table else PC3$msstats_table))[!is.na(pvalue) | is.infinite(log2FC), Protein],   # tested, plus absent in one condition
                     pipeline = k)
}
for (n in names(lists)) lists[[n]]$universe <- union(lists[[n]]$universe, c(lists[[n]]$up, lists[[n]]$down))   # proteins absent in one condition (MSstats) included
list_summary <- rbindlist(lapply(names(lists), function(n) data.table(List = n, Up = length(lists[[n]]$up), Down = length(lists[[n]]$down),
                                                                      Universe = length(lists[[n]]$universe))))
fwrite(list_summary, file.path(D2, "2.1_Lists_summary.csv"))
print(list_summary, row.names = FALSE)


# Section 3. Annotation: eggNOG-mapper for both organisms ----------------------
# One source for both proteomes, so that differences between Practicals 2 and 3 are not
# differences of annotation. eggNOG-mapper transfers annotation from orthologous groups:
# GO terms, KEGG orthologs and pathways, Pfam domains and COG functional categories.

D3 <- out_dir("3_Annotation")
read_emapper <- function(f) {
  x <- fread(file.path(DATA_DIR, f), skip = "#query", sep = "\t", quote = "", na.strings = c("-", ""), fill = TRUE)
  setnames(x, "#query", "query")
  x[!startsWith(query, "##")]
}
acc <- function(x) sub("^[^|]*\\|([^|]+)\\|.*$", "\\1", x)     # tr|W6QSR7|W6QSR7_ECTO5 -> W6QSR7
cog_names <- c(J = "Translation, ribosomal structure and biogenesis", A = "RNA processing and modification", K = "Transcription",
  L = "Replication, recombination and repair", B = "Chromatin structure and dynamics", D = "Cell cycle control, cell division, chromosome partitioning",
  Y = "Nuclear structure", V = "Defense mechanisms", T = "Signal transduction mechanisms", M = "Cell wall/membrane/envelope biogenesis",
  N = "Cell motility", Z = "Cytoskeleton", W = "Extracellular structures", U = "Intracellular trafficking, secretion, and vesicular transport",
  O = "Posttranslational modification, protein turnover, chaperones", X = "Mobilome: prophages, transposons", C = "Energy production and conversion",
  G = "Carbohydrate transport and metabolism", E = "Amino acid transport and metabolism", F = "Nucleotide transport and metabolism",
  H = "Coenzyme transport and metabolism", I = "Lipid transport and metabolism", P = "Inorganic ion transport and metabolism",
  Q = "Secondary metabolites biosynthesis, transport and catabolism", R = "General function prediction only", S = "Function unknown")
kegg_names <- if (file.exists(file.path(DATA_DIR, KEGG_MAP))) {   # file of Practical 3, or the names straight from KEGG
  unique(fread(file.path(DATA_DIR, KEGG_MAP))[, .(Term = Pathway, Name)])
} else tryCatch(fread("https://rest.kegg.jp/list/pathway", header = FALSE, sep = "\t", col.names = c("Term", "Name"))[, Term := sub("^path:", "", Term)],
                error = function(e) { cat("KEGG pathway names not available: map IDs are used as names.\n"); data.table(Term = character(), Name = character()) })
split_col <- function(d, col, pattern = ",") {   # one row per (gene, item) of a comma-separated column
  d <- d[!is.na(get(col)), .(id, v = as.character(get(col)))]
  d[, .(item = trimws(unlist(strsplit(v, pattern)))), by = id][item != ""]
}
build_sets <- function(em) {   # em: data.table with id + eggNOG columns -> list of (Term, id, Name) tables
  s <- list()
  go <- unique(split_col(em, "GOs"))
  info <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(go$item), columns = c("ONTOLOGY", "TERM"), keytype = "GOID"))
  go[, Ontology := info$ONTOLOGY[match(item, info$GOID)]][, Name := info$TERM[match(item, info$GOID)]]
  for (o in c("BP", "MF", "CC")) s[[paste0("GO_", o)]] <- unique(go[Ontology == o, .(Term = item, id, Name)])
  kg <- unique(split_col(em, "KEGG_Pathway"))
  kg[grepl("^[0-9]{5}$", item), item := paste0("map", item)]   # eggNOG-mapper 3 writes 00010, version 2 map00010
  kg <- kg[startsWith(item, "map")]
  s$KEGG <- merge(kg[, .(Term = item, id)], kegg_names, by = "Term", all.x = TRUE)[, Name := ifelse(is.na(Name), Term, Name)]
  s$Pfam <- unique(split_col(em, "PFAMs")[, item := sub("_[0-9]+_[0-9]+$", "", item)])[, .(Term = item, id, Name = item)]   # v3 adds _start_end
  if ("COG_category" %in% names(em)) {   # category letters (eggNOG-mapper 2); version 3 writes COG ids (COG0583) or "S" instead
    cg <- em[!is.na(COG_category) & grepl("^[A-Z]{1,5}$", COG_category), .(id, COG_category = as.character(COG_category))]
    cg <- cg[, .(item = unlist(strsplit(COG_category, ""))), by = id][item %in% names(cog_names) & item != "S"]   # S = function unknown
    if (uniqueN(cg$item) >= 5) s$COG <- unique(cg[, .(Term = item, id, Name = paste0(item, ": ", cog_names[item]))])
  }
  lapply(s[intersect(ANNOTATIONS, names(s))], unique)
}
# Chlamydomonas: proteins -> genes from the FASTA headers (gene:CHLRE_...), all isoforms of a gene pooled
hdr <- grep("^>", readLines(file.path(DATA_DIR, CHLAMY_FASTA)), value = TRUE)
p2g <- data.table(query = sub("^>(\\S+).*", "\\1", hdr), id = sub(".*\\bgene:(\\S+).*", "\\1", hdr))
em_c <- merge(read_emapper(EGGNOG_CHLAMY), p2g, by = "query")
em_e <- read_emapper(EGGNOG_EOLEO)[, id := acc(query)]
sets <- list(chlamy = build_sets(em_c), eoleo = build_sets(em_e))
organism <- function(pipeline) if (grepl("^(Ensembl|Trinity)", pipeline)) "chlamy" else "eoleo"
coverage <- rbindlist(lapply(names(lists), function(n) {
  u <- lists[[n]]$universe; s <- sets[[organism(lists[[n]]$pipeline)]]
  data.table(List = n, Universe = length(u), t(sapply(s, function(a) round(100 * mean(u %in% a$id), 1))),
             Any = round(100 * mean(u %in% unique(unlist(lapply(s, `[[`, "id")))), 1))
}), fill = TRUE)
fwrite(coverage, file.path(D3, "3_Annotation_coverage_percent.csv"))
cat("Annotated genes/proteins in each universe (%):\n"); print(coverage, row.names = FALSE)
fwrite(rbindlist(lapply(names(sets), function(o) rbindlist(lapply(names(sets[[o]]), function(a)
  data.table(Organism = o, Annotation = a, Terms = uniqueN(sets[[o]][[a]]$Term), Genes = uniqueN(sets[[o]][[a]]$id)))))), file.path(D3, "3_Annotation_sets.csv"))


# Section 4. ORA and GSEA of every pipeline, same settings ---------------------

D4 <- out_dir("4_Enrichment")
ora <- function(ids, universe, s) {   # hypergeometric test, universe = annotated genes tested in the pipeline
  u <- intersect(universe, s$id); x <- intersect(ids, u); if (length(x) < ORA_MIN) return(NULL)
  t <- s[id %in% u, .(K = uniqueN(id), k = uniqueN(id[id %in% x]), genes = paste(sort(unique(id[id %in% x])), collapse = "/")), by = .(Term, Name)]
  t <- t[K >= ORA_MIN & K <= MAX_SIZE]
  t[, `:=`(n = length(x), N = length(u))][, pvalue := phyper(k - 1, K, N - K, n, lower.tail = FALSE)][, padj := p.adjust(pvalue, "BH")]
  t[, Ratio := sprintf("%d/%d", k, n)][, Background := sprintf("%d/%d", K, N)]
  t[k >= ORA_MIN & padj < FDR][order(padj)]
}
gsea <- function(z, s, seed = SEED) {  # fgsea on the signed z-scores
  z <- z[is.finite(z) & !duplicated(names(z))]; z <- sort(z, decreasing = TRUE)
  pw <- split(s$id, s$Term); pw <- lapply(pw, function(v) intersect(unique(v), names(z)))
  pw <- pw[lengths(pw) >= MIN_SIZE & lengths(pw) <= MAX_SIZE]
  set.seed(seed)
  r <- suppressWarnings(fgsea::fgseaMultilevel(pw, z, minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0, BPPARAM = BiocParallel::SerialParam()))
  r <- as.data.table(r)[, .(Term = pathway, size, NES, pvalue = pval, padj)]
  merge(r, unique(s[, .(Term, Name)]), by = "Term", all.x = TRUE)
}
ora_all <- rbindlist(lapply(names(lists), function(ln) rbindlist(lapply(c("up", "down"), function(dirn) {
  s <- sets[[organism(lists[[ln]]$pipeline)]]
  rbindlist(lapply(names(s), function(a) { r <- ora(lists[[ln]][[dirn]], lists[[ln]]$universe, s[[a]]); if (!is.null(r) && nrow(r)) r[, `:=`(List = ln, Direction = dirn, Annotation = a)] }))
}))), fill = TRUE)
if (nrow(ora_all)) { setcolorder(ora_all, c("List", "Direction", "Annotation", "Term", "Name", "Ratio", "Background", "padj")); fwrite(ora_all[, -c("K", "k", "n", "N")], file.path(D4, "4.1_ORA_all.csv")) }
gsea_all <- rbindlist(lapply(names(stats), function(n) {
  s <- sets[[organism(n)]]; z <- setNames(stats[[n]]$z, stats[[n]]$id)
  rbindlist(lapply(names(s), function(a) gsea(z, s[[a]])[, `:=`(Pipeline = n, Annotation = a)]))
}))
fwrite(gsea_all[order(padj)], file.path(D4, "4.2_GSEA_all.csv"))
ora_sum <- if (nrow(ora_all)) dcast(ora_all, List + Direction ~ Annotation, fun.aggregate = length, value.var = "Term") else data.table()
gsea_sum <- dcast(gsea_all[padj < FDR], Pipeline ~ Annotation, fun.aggregate = length, value.var = "Term")
fwrite(ora_sum, file.path(D4, "4.1_ORA_summary.csv")); fwrite(gsea_sum, file.path(D4, "4.2_GSEA_summary.csv"))
cat("ORA: significant terms\n"); print(ora_sum); cat("GSEA: significant terms\n"); print(gsea_sum)


# Section 5. Descriptive comparison of the pipelines ---------------------------

D5 <- out_dir("5_Comparison")
# 5.1 Jaccard index of the significant ORA terms (both directions, all annotations) between lists
jac <- function(a, b) if (length(union(a, b)) == 0) NA else length(intersect(a, b)) / length(union(a, b))
for (o in c("chlamy", "eoleo")) {
  ln <- names(lists)[sapply(names(lists), function(n) organism(lists[[n]]$pipeline) == o)]
  terms <- lapply(ln, function(l) if (nrow(ora_all)) ora_all[List == l, paste(Direction, Annotation, Term)] else character(0)); names(terms) <- ln   # not "n": ora_all has a column n
  J <- outer(ln, ln, Vectorize(function(a, b) jac(terms[[a]], terms[[b]]))); dimnames(J) <- list(ln, ln)
  fwrite(data.table(List = ln, round(J, 2)), file.path(D5, paste0("5.1_Jaccard_ORA_", o, ".csv")))
  if (length(unique(J[is.finite(J)])) > 1)
    pheatmap(J, display_numbers = TRUE, number_format = "%.2f", cluster_rows = FALSE, cluster_cols = FALSE, na_col = "grey90",
             main = paste("Jaccard index of significant ORA terms -", o), filename = file.path(D5, paste0("5.1_Jaccard_ORA_", o, ".png")), width = 7, height = 6)
}
# 5.2 NES of every pair of pipelines (all terms tested in both)
nes_pairs <- rbindlist(lapply(names(PAIRS), function(pn) {
  a <- PAIRS[[pn]][1]; b <- PAIRS[[pn]][2]
  m <- merge(gsea_all[Pipeline == a, .(Annotation, Term, Name, NES_A = NES, padj_A = padj)],
             gsea_all[Pipeline == b, .(Annotation, Term, NES_B = NES, padj_B = padj)], by = c("Annotation", "Term"))
  m[, Pair := pn]
}))
nes_pairs[, Significant := ifelse(padj_A < FDR & padj_B < FDR, "both", ifelse(padj_A < FDR, "first only", ifelse(padj_B < FDR, "second only", "neither")))]
fwrite(nes_pairs, file.path(D5, "5.2_NES_pairs.csv"))
# NES of a term that is not enriched sits near +1 or -1, never near 0 (fgsea normalises each sign
# by its own null), so noise terms whose sign flips form the two off-diagonal clouds of the plot.
# For the terms enriched in at least one pipeline, the table gives how often both pipelines agree on the direction.
nes_cor <- nes_pairs[, .(Terms = .N, Correlation = round(cor(NES_A, NES_B), 2),
                         Enriched_any = sum(Significant != "neither"),
                         Same_direction_pct = round(100 * mean(sign(NES_A[Significant != "neither"]) == sign(NES_B[Significant != "neither"])), 1),
                         Both = sum(Significant == "both"),
                         First_only = sum(Significant == "first only"), Second_only = sum(Significant == "second only")), by = Pair]
fwrite(nes_cor, file.path(D5, "5.2_NES_correlation.csv")); print(nes_cor, row.names = FALSE)
p <- ggplot(nes_pairs, aes(NES_A, NES_B, colour = Significant)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(linetype = "dashed", colour = "grey50") + geom_point(size = 1, alpha = 0.6) + facet_wrap(~ Pair, ncol = 3) +
  scale_colour_manual(values = c(both = "#14213D", "first only" = "#B5541B", "second only" = "#2A6F97", neither = "grey75")) +
  theme_bw() + labs(x = "NES, first pipeline of the pair", y = "NES, second pipeline", colour = paste("GSEA adj. p <", FDR), title = "Same terms, two pipelines")
save_plot(p, D5, "5.2_NES_pairs.png", 13, 8)


# Section 6. Differential enrichment: where do two pipelines disagree? --------
# For the genes tested in both pipelines, d = z(first) - z(second). Ranking the genes by d
# and running GSEA finds the functions whose genes systematically get a higher (NES > 0)
# or lower (NES < 0) score in the first pipeline than in the second. Caution: genes are
# permuted, so correlation between the genes of a function makes this test optimistic;
# section 7 answers the same question by resampling the samples.

D6 <- out_dir("6_Differential")
diff_all <- rbindlist(lapply(names(PAIRS), function(pn) {
  a <- PAIRS[[pn]][1]; b <- PAIRS[[pn]][2]
  m <- merge(stats[[a]][is.finite(z), .(id, za = z)], stats[[b]][is.finite(z), .(id, zb = z)], by = "id")
  d <- setNames(m$za - m$zb, m$id); s <- sets[[organism(a)]]
  rbindlist(lapply(names(s), function(an) gsea(d, s[[an]])[, `:=`(Pair = pn, Annotation = an, Genes = nrow(m))]))
}))
fwrite(diff_all[order(padj)], file.path(D6, "6_Differential_enrichment.csv"))
diff_sum <- diff_all[, .(Genes = Genes[1], Terms_tested = .N, Discordant = sum(padj < FDR, na.rm = TRUE), Higher_in_first = sum(padj < FDR & NES > 0, na.rm = TRUE),
                         Higher_in_second = sum(padj < FDR & NES < 0, na.rm = TRUE)), by = Pair]
fwrite(diff_sum, file.path(D6, "6_Differential_summary.csv")); print(diff_sum, row.names = FALSE)


# Section 7. Resampling: leave one sample out ---------------------------------
# Each pipeline is run again without one sample at a time (jackknife). For each function the
# difference of NES between two pipelines is recomputed every time; its spread is the
# uncertainty that comes only from the choice of samples. A difference larger than that
# (|z| = |difference| / jackknife SE, t with n - 1 df, adjusted p < FDR) means that the pipeline
# changes the enrichment of that function beyond sampling noise. As in Practical 3, significance
# is not enough: a "pipeline effect" also needs |difference| >= DELTA_NES and the function
# enriched in at least one of the two pipelines.
# LONG STEP (5-15 min): results are saved and reused.
#   PC3: the protein values of each pipeline are kept (normalisation, imputation and MSstats
#   summarisation are not repeated); only the test is, without the sample left out.

D7 <- out_dir("7_Resampling")
pipes <- unique(unlist(PAIRS))
subsets <- function(pipe) if (organism(pipe) == "chlamy") samples2$sample else colnames(prot[[pipe]])
cache7 <- file.path(D7, "Resampling_NES.rds")
key7 <- paste(c("jk-v1", pipes, FDR, MIN_SIZE, MAX_SIZE, SEED, sapply(pipes, function(p) paste(subsets(p), collapse = ","))), collapse = ";")
jk <- if (file.exists(cache7)) readRDS(cache7) else list()
if (!identical(jk$key, key7)) jk <- list(key = key7)
for (p in pipes) {
  if (!is.null(jk[[p]])) next
  cat("Resampling", p, "...\n"); t0 <- Sys.time()
  s <- sets[[organism(p)]]
  jk[[p]] <- rbindlist(lapply(subsets(p), function(out) {
    st <- if (organism(p) == "chlamy") de_pc2(sub("_.*", "", p), sub(".*_", "", p), setdiff(samples2$sample, out)) else de_pc3(p, out)
    z <- setNames(zscore(st$log2FC, st$pvalue), st$id)
    rbindlist(lapply(names(s), function(a) gsea(z, s[[a]])[, .(Annotation = a, Term, NES)]))[, Left_out := out]
  }))
  saveRDS(jk, cache7)
  cat("  done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}
jk_all <- rbindlist(lapply(names(PAIRS), function(pn) {
  a <- PAIRS[[pn]][1]; b <- PAIRS[[pn]][2]
  full <- nes_pairs[Pair == pn, .(Annotation, Term, Name, Delta = NES_A - NES_B, NES_A, NES_B, Significant)]
  if (organism(a) == "chlamy") {   # PC2: the same samples are left out in both pipelines (paired)
    m <- merge(jk[[a]][, .(Annotation, Term, Left_out, A = NES)], jk[[b]][, .(Annotation, Term, Left_out, B = NES)], by = c("Annotation", "Term", "Left_out"))
  } else {                         # PC3: both pipelines have the same samples (same labels)
    m <- merge(jk[[a]][, .(Annotation, Term, Left_out, A = NES)], jk[[b]][, .(Annotation, Term, Left_out, B = NES)], by = c("Annotation", "Term", "Left_out"))
  }
  n <- length(subsets(a))
  se <- m[, .(n_sub = .N, SE = sqrt((.N - 1) / .N * sum(((A - B) - mean(A - B))^2)), Same_sign = NA_real_, d_sub = list(A - B)), by = .(Annotation, Term)][n_sub == n]
  r <- merge(full, se, by = c("Annotation", "Term"))
  r[, Same_sign := mapply(function(d, D) mean(sign(d) == sign(D)), d_sub, Delta)][, d_sub := NULL]
  r[, z := Delta / pmax(SE, 1e-6)][, pvalue := 2 * pt(-abs(z), df = n - 1)][, padj := p.adjust(pvalue, "BH"), by = Annotation]
  # Same logic as significance + relevance in Practical 3: the difference must be beyond sampling
  # noise, large (|delta NES| >= DELTA_NES) and about a function enriched in at least one pipeline.
  # A term going from NES +1 to -1 without being enriched in either is a real but irrelevant difference.
  r[, Pipeline_effect := padj < FDR & abs(Delta) >= DELTA_NES & Significant != "neither"]
  r[, Pair := pn]
}))
fwrite(jk_all[order(padj)], file.path(D7, "7_Resampling_NES_difference.csv"))
jk_sum <- jk_all[, .(Terms = .N, Beyond_noise = sum(padj < FDR), Pipeline_effect = sum(Pipeline_effect), Median_abs_Delta = round(median(abs(Delta)), 2),
                     Median_SE = round(median(SE), 2)), by = Pair]
fwrite(jk_sum, file.path(D7, "7_Resampling_summary.csv")); print(jk_sum, row.names = FALSE)
jk_all[, Class := ifelse(Pipeline_effect, "pipeline effect", ifelse(Significant != "neither", "enriched in one or both: within noise or small", "not enriched in either"))]
p <- ggplot(jk_all[order(Class != "not enriched in either")], aes(Delta, pmin(abs(z), 50), colour = Class)) + geom_point(size = 0.9, alpha = 0.7) + facet_wrap(~ Pair, ncol = 3) +
  scale_colour_manual(values = c("pipeline effect" = "#B5541B", "enriched in one or both: within noise or small" = "#14213D", "not enriched in either" = "grey80")) +
  geom_vline(xintercept = c(-DELTA_NES, DELTA_NES), linetype = 2, colour = "grey50") +
  theme_bw() + labs(x = "NES difference (first - second pipeline), all samples", y = "|difference| / jackknife SE (capped at 50)", colour = NULL,
                    title = "Is the pipeline effect larger than the effect of the samples?")
save_plot(p, D7, "7_Resampling_volcano.png", 13, 8)


# Section 8. Files for Cytoscape / ClueGO ---------------------------------------
# One text file per list and direction (one identifier per line), and the annotation of
# each organism as GMT files (one line per term: term, description, genes).

D8 <- out_dir("8_ClueGO")
for (n in names(lists)) for (dirn in c("up", "down"))
  writeLines(lists[[n]][[dirn]], file.path(D8, paste0("List_", n, "_", dirn, ".txt")))
for (o in names(sets)) for (a in names(sets[[o]])) {
  g <- sets[[o]][[a]][, .(genes = paste(unique(id), collapse = "\t")), by = .(Term, Name)]
  writeLines(paste(g$Term, g$Name, g$genes, sep = "\t"), file.path(D8, paste0("Annotation_", o, "_", a, ".gmt")))
}


# Section 9. Session information --------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\nDone. All results are in:", OUT, "\n")
