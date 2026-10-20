# =============================================================================
# PC4b_MOFA.R
# Systems Molecular Biology - Practical 4, part 2 (Multi-omics integration with MOFA) - 2026-2027
#
# MOFA (Multi-Omics Factor Analysis) finds a few hidden "factors" that explain the variation
# shared by several omics measured on the same samples, and tells for each factor which omics
# it affects and which features (genes, proteins, lipids) drive it.
#   Part A - simulated data: the true factors are known, so we can check what MOFA recovers.
#   Part B - real data: Chlamydomonas acclimated to low, moderate and high diurnal light,
#            sampled at five times of the day (3 x 5 x 3 = 45 samples), with the transcriptome,
#            the proteome (TMT) and the lipidome measured on the SAME samples.
#            Data: Dupuis, Ojeda et al. (2025) The Plant Cell 37: koaf086, Supplementary Data Sets.
#
# The models are trained with mofapy2 (Python) by the instructor (PC4b_train_MOFA.py) and given
# in ForStudents_PC4.zip: this script only prepares the data and reads the trained models,
# so you do not need Python.
#
# HOW TO USE
#   1. ForStudents_PC4.zip already created the folder BMSi/PC4_Enrichment with a MOFA subfolder.
#   2. Run the script section by section in RStudio (cursor in a section, Ctrl+Alt+T).
#   3. All results are written to the folder "results_MOFA" inside DATA_DIR.
#
# RESULTS (folder "results_MOFA")
#   1_Simulated      the simulated data and the true factors
#   2_Data           preparation of the real data: filters, TMT plexes, PCA before/after
#   3_Simulated_fit  what MOFA recovered from the simulated data
#   4_Variance       variance explained by each factor in each omics (real data)
#   5_Factors        factors against the design (light, time, replicate, plex) and physiology
#   6_Weights        features that drive each factor; enrichment of the weights (eggNOG)
# =============================================================================


# CONFIGURATION: the only part you may need to edit ---------------------------

DATA_DIR <- ""      # empty = automatic: Windows C:/Users/<you>/BMSi/PC4_Enrichment, macOS/Linux ~/BMSi/PC4_Enrichment
MOFA_DIR <- "MOFA"  # inside DATA_DIR: supplementary data and trained models
FILES <- list(metadata = "SupplementaryDataSet1_OmicsSampleMetadata.xlsx",
              omics    = "SupplementaryDataSet2_mRNA_protein_abundances.xlsx",
              lipids   = "SupplementaryDataSet7_Lipidome_MedianNormdPeakIntensity.xlsx",
              pigments = "SupplementaryDataSet8_Pigments_pmolPerMillionCells.xlsx")
# Annotation of Practical 4, part 1 (inside DATA_DIR)
EGGNOG_CHLAMY <- "chlamy.emapper.annotations"
CHLAMY_FASTA  <- "Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.pep.all.fa"
KEGG_MAP      <- "6_KEGG_KO_pathways.tsv"

# Data preparation (part B)
MIN_FPKM    <- 1      # mRNA: mean FPKM over the 45 samples
N_FEATURES  <- 2000   # most variable mRNAs and proteins kept for MOFA (lipids: all)
MAX_NA_PROT <- 4      # proteins with more missing values than this are left out

# Interpretation
FDR <- 0.05; MIN_SIZE <- 10; MAX_SIZE <- 500   # GSEA of the weights
TOP_N <- 15                                   # features listed per factor
SEED <- 2026


# Section 0. Packages ----------------------------------------------------------

if (getRversion() < "4.3.0") stop("This script needs R 4.3.0 or newer. You have R ", getRversion(), ".")
options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("data.table", "ggplot2", "pheatmap", "readxl", "fgsea", "GO.db", "AnnotationDbi", "BiocParallel")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
}
for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("Package check:\n"); for (p in pkgs) cat(sprintf("  %-14s %s\n", p, packageVersion(p)))
set.seed(SEED)

if (DATA_DIR == "") {
  home <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
  DATA_DIR <- file.path(home, "BMSi", "PC4_Enrichment")
}
DATA_DIR <- normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE)
MD <- file.path(DATA_DIR, MOFA_DIR)
OUT <- file.path(DATA_DIR, "results_MOFA")
out_dir <- function(name) { d <- file.path(OUT, name); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
save_plot <- function(p, folder, file, w = 8, h = 6) ggsave(file.path(folder, file), p, width = w, height = h, dpi = 300)
dir.create(file.path(MD, "input"), recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))
light_cols <- c(LL = "grey45", ML = "#E1A100", HL = "#C0392B")


# Section 1. Part A: simulated data with known factors --------------------------
# 45 samples with the design of the real experiment (3 light levels x 5 times x 3 replicates)
# and three omics. Four true factors, each with its own pattern of activity:
#   F1 "light"   : depends on the light level; affects mRNA, protein and lipids
#   F2 "time"    : a daily rhythm;              affects mRNA and protein
#   F3 "mRNA"    : sample-to-sample variation;  affects mRNA only
#   F4 "lipids"  : sample-to-sample variation;  affects lipids only
# In each omics only some features respond to a factor (sparse weights). 5% of the protein
# values and the lipids of three samples are missing, as happens with real data.

D1 <- out_dir("1_Simulated")
set.seed(SEED)   # the trained model was fitted to exactly these simulated data
meta_sim <- CJ(Light = c("LL", "ML", "HL"), Time = c("CTm02", "CTp02", "CTp06", "CTp10", "CTm10"), Rep = 0:2, sorted = FALSE)
meta_sim[, sample := sprintf("S%02d", .I)]
hours <- c(CTm02 = 22, CTp02 = 2, CTp06 = 6, CTp10 = 10, CTm10 = 14)
Z_true <- cbind(F1_light = c(LL = -1.2, ML = 0, HL = 1.2)[meta_sim$Light] + rnorm(45, 0, 0.25),
                F2_time  = sin(2 * pi * (hours[meta_sim$Time] - 4) / 24) * 1.2 + rnorm(45, 0, 0.25),
                F3_mRNA  = rnorm(45), F4_lipids = rnorm(45))
rownames(Z_true) <- meta_sim$sample
views_sim <- list(mRNA = 1000, protein = 600, lipids = 150)
active <- rbind(mRNA = c(1, 1, 1, 0), protein = c(1, 1, 0, 0), lipids = c(1, 0, 0, 1))   # which factor acts on which omics
colnames(active) <- colnames(Z_true)
noise_sd <- c(mRNA = 1, protein = 1.3, lipids = 0.8)
sim <- lapply(names(views_sim), function(v) {
  D <- views_sim[[v]]
  W <- sapply(seq_len(ncol(Z_true)), function(k) if (active[v, k] == 1) rbinom(D, 1, 0.2) * rnorm(D, 0, 1.5) else rep(0, D))
  Y <- Z_true %*% t(W) + matrix(rnorm(45 * D, 0, noise_sd[v]), 45, D)
  colnames(Y) <- sprintf("%s_%04d", v, seq_len(D)); rownames(W) <- colnames(Y); colnames(W) <- colnames(Z_true)
  if (v == "protein") Y[sample(length(Y), round(0.05 * length(Y)))] <- NA
  if (v == "lipids") Y[sample(45, 3), ] <- NA
  list(Y = Y, W = W)
})
names(sim) <- names(views_sim)
long_sim <- rbindlist(lapply(names(sim), function(v) {
  m <- as.data.table(sim[[v]]$Y, keep.rownames = "sample")
  melt(m, id.vars = "sample", variable.name = "feature", variable.factor = FALSE)[!is.na(value)][, view := v]
}))[, group := "group1"]
fwrite(long_sim[, .(sample, group, feature, view, value)], file.path(MD, "input", "sim_long.csv"))
fwrite(as.data.table(Z_true, keep.rownames = "sample")[meta_sim, on = "sample"], file.path(D1, "1_True_factors.csv"))
fwrite(rbindlist(lapply(names(sim), function(v) as.data.table(sim[[v]]$W, keep.rownames = "feature")[, view := v])), file.path(D1, "1_True_weights.csv"))
pheatmap(active, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.0f", legend = FALSE,
         color = c("grey95", "#14213D"), main = "Simulated design: which factor acts on which omics",
         filename = file.path(D1, "1_True_activity.png"), width = 6, height = 3)
cat("Simulated data:", paste(names(views_sim), unlist(views_sim), "features", collapse = "; "), "| 45 samples\n")


# Section 2. Part B: the real data ----------------------------------------------
# 2.1 Read the three omics and the design. 2.2 Filter and transform each omics.
# 2.3 TMT plexes: the 45 protein samples were run in three TMT plexes, and the plex dominates
#     the protein values (the authors report 98% of the variance in the first two PCs).
#     The plex is not in the metadata: it is recovered from the PCA. Each plex holds one replicate
#     of every light x time combination, so centring each protein within its plex removes the
#     plex without removing biology.
# 2.4 Write the input for MOFA: once with the plex corrected, once without (to compare).

D2 <- out_dir("2_Data")
need <- file.path(MD, unlist(FILES))
if (!all(file.exists(need))) stop("Missing in ", MD, ":\n  ", paste(basename(need[!file.exists(need)]), collapse = "\n  "))
meta <- as.data.table(read_excel(file.path(MD, FILES$metadata), sheet = "OmicsSampleMetadata"))
meta[, Light := factor(Population, levels = c("LL", "ML", "HL"))][, Time := factor(TimeCode, levels = c("CTm02", "CTp02", "CTp06", "CTp10", "CTm10"))]
meta[, Replicate := factor(ExperimentalReplicate)]
om <- as.data.table(read_excel(file.path(MD, FILES$omics), sheet = "mRNA_FPKMs_protein_MASICs", na = "NA"))
annot <- om[, .(id = Cre_GID, symbol = Cre_symbol, description = Description)]
to_matrix <- function(dt, cols, ids) { m <- as.matrix(dt[, ..cols]); storage.mode(m) <- "double"; rownames(m) <- ids; m }
mrna <- to_matrix(om, meta$mRNA_ColumnHeader, om$Cre_GID); colnames(mrna) <- meta$SampleID
prot <- to_matrix(om, meta$protein_ColumnHeader, om$Cre_GID); colnames(prot) <- meta$SampleID
prot <- prot[rowSums(!is.na(prot)) > 0, ]; mrna <- mrna[rowSums(!is.na(mrna)) > 0, ]
lip <- as.data.table(read_excel(file.path(MD, FILES$lipids), sheet = "Lipids_MedNormPeakIntensity"))
lipids <- to_matrix(lip, meta$SampleID, lip$LipidID)
pig <- as.data.table(read_excel(file.path(MD, FILES$pigments), sheet = "Pigments_pmolPerMillionCells"))
pigments <- t(to_matrix(pig, intersect(meta$SampleID, names(pig)), pig$Pigment))   # samples x pigments (replicates 0-2)
cat("Read:", nrow(mrna), "mRNAs,", nrow(prot), "proteins,", nrow(lipids), "lipids,", ncol(pigments), "pigments; 45 samples\n")

# 2.2 Transform and filter
mrna_l <- log2(mrna + 1)
mrna_l <- mrna_l[rowMeans(mrna) >= MIN_FPKM, ]
prot_l <- log2(prot); prot_l <- prot_l[rowSums(is.na(prot_l)) <= MAX_NA_PROT, ]

# 2.3 TMT plexes
pca <- function(m) { m <- m[complete.cases(m), ]; p <- prcomp(t(m - rowMeans(m))); list(x = p$x[, 1:2], var = round(100 * p$sdev[1:2]^2 / sum(p$sdev^2), 1)) }
p_raw <- pca(prot_l)
meta[, Plex := factor(paste0("plex", kmeans(p_raw$x, centers = 3, nstart = 25)$cluster[match(SampleID, rownames(p_raw$x))]))]
bal <- dcast(meta, Light + Time ~ Plex, fun.aggregate = length, value.var = "SampleID")
fwrite(bal, file.path(D2, "2.3_Plex_balance.csv"))
if (!all(as.matrix(bal[, -(1:2)]) == 1)) warning("The plexes are not balanced (one sample of each light x time per plex): centring within plex may remove biology.")
prot_c <- prot_l
for (pl in levels(meta$Plex)) { s <- meta[Plex == pl, SampleID]; prot_c[, s] <- prot_c[, s] - rowMeans(prot_c[, s], na.rm = TRUE) }
p_cor <- pca(prot_c)
pc_dt <- rbind(data.table(meta, PC1 = p_raw$x[meta$SampleID, 1], PC2 = p_raw$x[meta$SampleID, 2], Data = sprintf("As published (PC1 %s%%, PC2 %s%%)", p_raw$var[1], p_raw$var[2])),
               data.table(meta, PC1 = p_cor$x[meta$SampleID, 1], PC2 = p_cor$x[meta$SampleID, 2], Data = sprintf("Centred within plex (PC1 %s%%, PC2 %s%%)", p_cor$var[1], p_cor$var[2])))
p <- ggplot(pc_dt, aes(PC1, PC2, colour = Light, shape = Plex)) + geom_point(size = 2.5) + facet_wrap(~Data, scales = "free") +
  scale_colour_manual(values = light_cols) + labs(title = "Proteome PCA: the TMT plex hides the biology")
save_plot(p, D2, "2.3_Protein_PCA_plex.png", 11, 4.8)

# Keep the most variable features (variance after the corrections)
top_var <- function(m, n) m[order(apply(m, 1, var, na.rm = TRUE), decreasing = TRUE)[seq_len(min(n, nrow(m)))], ]
mrna_v <- top_var(mrna_l, N_FEATURES); prot_v <- top_var(prot_c, N_FEATURES); prot_raw_v <- prot_l[rownames(prot_v), ]
summary_dt <- data.table(Omics = c("mRNA", "protein", "lipids"),
                         Read = c(nrow(mrna), nrow(prot), nrow(lipids)), After_filters = c(nrow(mrna_l), nrow(prot_l), nrow(lipids)),
                         In_MOFA = c(nrow(mrna_v), nrow(prot_v), nrow(lipids)),
                         Missing_pct = round(100 * c(mean(is.na(mrna_v)), mean(is.na(prot_v)), mean(is.na(lipids))), 2))
fwrite(summary_dt, file.path(D2, "2.2_Features.csv")); print(summary_dt, row.names = FALSE)

# 2.4 Input for MOFA (long format: sample, group, feature, view, value)
to_long <- function(m, v) { d <- as.data.table(m, keep.rownames = "feature"); d <- melt(d, id.vars = "feature", variable.name = "sample", variable.factor = FALSE); d[!is.na(value)][, view := v] }
real <- rbind(to_long(mrna_v, "mRNA"), to_long(prot_v, "protein"), to_long(lipids, "lipids"))[, group := "group1"]
real_raw <- rbind(to_long(mrna_v, "mRNA"), to_long(prot_raw_v, "protein"), to_long(lipids, "lipids"))[, group := "group1"]
fwrite(real[, .(sample, group, feature, view, value)], file.path(MD, "input", "real_long.csv"))
fwrite(real_raw[, .(sample, group, feature, view, value)], file.path(MD, "input", "real_plex_uncorrected_long.csv"))
fwrite(cbind(meta[, .(SampleID, Light, LightIntensity, Time, ZTTime, Replicate, Plex, TemperatureC, CellDensity_Per_mL)],
             as.data.table(pigments[meta$SampleID, ])), file.path(D2, "2.4_Sample_metadata.csv"))


# Section 3. The trained models --------------------------------------------------
# Training is done once by the instructor with PC4b_train_MOFA.py (mofapy2): 10 factors to
# start with; factors that explain less than 1% of the variance in every omics are dropped
# during training; the omics are scaled to the same total variance, so that the largest omics
# does not dominate. Three models are provided:
#   sim             - the simulated data of section 1
#   real            - the real data, plex corrected
#   real_uncorrected- the real data, plex NOT corrected (to see what a batch effect does)

models <- c("sim", "real", "real_uncorrected")
read_model <- function(m) {
  d <- file.path(MD, "models", m)
  if (!file.exists(file.path(d, "factors.csv"))) stop("Trained model '", m, "' not found in ", d, ". It is in ForStudents_PC4.zip (instructor: run PC4b_train_MOFA.py).")
  list(Z = fread(file.path(d, "factors.csv")), W = fread(file.path(d, "weights.csv")), R2 = fread(file.path(d, "variance_explained.csv")))
}
fit <- setNames(lapply(models, read_model), models)
for (m in models) cat(sprintf("Model %-17s %d factors\n", m, ncol(fit[[m]]$Z) - 1))
fmat <- function(Z) { m <- as.matrix(Z[, -1]); rownames(m) <- Z$sample; m }


# Section 4. Part A: what did MOFA recover from the simulated data? -------------

D3 <- out_dir("3_Simulated_fit")
Zs <- fmat(fit$sim$Z)[rownames(Z_true), ]
C <- cor(Zs, Z_true)
pheatmap(round(C, 2), display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(-1, 1, length.out = 101),
         color = colorRampPalette(c("#2A6F97", "white", "#B5541B"))(100), main = "Correlation of MOFA factors with the true factors",
         filename = file.path(D3, "3.1_Factors_vs_truth.png"), width = 6, height = 4.5)
match_dt <- data.table(MOFA_factor = rownames(C), Best_true_factor = colnames(C)[apply(abs(C), 1, which.max)], r = round(apply(C, 1, function(x) x[which.max(abs(x))]), 2))
fwrite(match_dt, file.path(D3, "3.1_Factor_matching.csv")); print(match_dt, row.names = FALSE)
r2s <- dcast(fit$sim$R2, view ~ factor, value.var = "R2")
fwrite(r2s, file.path(D3, "3.2_Variance_explained.csv"))
r2_plot <- function(R2, title) ggplot(R2, aes(factor(factor, levels = unique(factor)), view, fill = R2)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(R2 >= 1, round(R2), "")), size = 3.5) + scale_fill_gradient(low = "grey97", high = "#14213D", name = "% variance") +
  labs(x = NULL, y = NULL, title = title) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(r2_plot(fit$sim$R2, "Simulated data: variance explained (%) by each factor in each omics"), D3, "3.2_Variance_explained.png", 8, 3.5)
# Weights: are the features that truly respond to a factor the ones with the largest MOFA weights?
Wt <- fread(file.path(D1, "1_True_weights.csv"))
wrec <- rbindlist(lapply(seq_len(nrow(match_dt)), function(i) {
  f <- match_dt$MOFA_factor[i]; tf <- match_dt$Best_true_factor[i]
  x <- merge(fit$sim$W[, c("feature", "view", f), with = FALSE], Wt[, c("feature", "view", tf), with = FALSE], by = c("feature", "view"))
  x[, .(MOFA_factor = f, True_factor = tf, Features = .N, Truly_active = sum(get(tf) != 0),
        AUC = if (any(get(tf) != 0) && any(get(tf) == 0)) round(wilcox.test(abs(get(f))[get(tf) != 0], abs(get(f))[get(tf) == 0])$statistic / (sum(get(tf) != 0) * sum(get(tf) == 0)), 2) else NA_real_), by = view]
}))
fwrite(wrec, file.path(D3, "3.3_Weights_recovery_AUC.csv")); print(wrec, row.names = FALSE)


# Section 5. Part B: variance explained and what each factor means ---------------

D4 <- out_dir("4_Variance")
for (m in c("real", "real_uncorrected")) {
  fwrite(dcast(fit[[m]]$R2, view ~ factor, value.var = "R2"), file.path(D4, paste0("4_Variance_explained_", m, ".csv")))
  save_plot(r2_plot(fit[[m]]$R2, paste0("Variance explained (%) - ", ifelse(m == "real", "plex corrected", "plex NOT corrected"))), D4, paste0("4_Variance_explained_", m, ".png"), 9, 3.5)
}
tot <- rbindlist(lapply(c("real", "real_uncorrected"), function(m) fit[[m]]$R2[, .(Total_R2 = round(sum(R2), 1)), by = view][, Model := m]))
fwrite(tot, file.path(D4, "4_Total_variance_explained.csv")); print(dcast(tot, view ~ Model, value.var = "Total_R2"), row.names = FALSE)

D5 <- out_dir("5_Factors")
md <- fread(file.path(D2, "2.4_Sample_metadata.csv"))
assoc <- function(Z) {
  Z <- fmat(Z); md2 <- md[match(rownames(Z), SampleID)]
  rbindlist(lapply(colnames(Z), function(f) {
    y <- Z[, f]
    r2 <- function(x) summary(lm(y ~ x))$r.squared
    data.table(Factor = f,
               Light = r2(factor(md2$Light)), Time = r2(factor(md2$Time)), Replicate = r2(factor(md2$Replicate)), Plex = r2(factor(md2$Plex)),
               Temperature = cor(y, md2$TemperatureC)^2, log_CellDensity = cor(y, log10(md2$CellDensity_Per_mL))^2,
               Chl_a = cor(y, md2$Chl_a, use = "complete.obs")^2, Chl_a_b_ratio = cor(y, md2$Chl_a / md2$Chl_b, use = "complete.obs")^2,
               Zeaxanthin = cor(y, md2$zeaxanthin, use = "complete.obs")^2, Lutein = cor(y, md2$lutein, use = "complete.obs")^2)
  }))
}
for (m in c("real", "real_uncorrected")) {
  a <- assoc(fit[[m]]$Z); fwrite(a[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], file.path(D5, paste0("5.1_Factor_associations_", m, ".csv")))
  am <- as.matrix(a[, -1]); rownames(am) <- a$Factor
  pheatmap(round(t(am), 2), display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("grey97", "#14213D"))(100),
           breaks = seq(0, 1, length.out = 101), main = paste("R2 of each factor with the design and physiology -", m),
           filename = file.path(D5, paste0("5.1_Factor_associations_", m, ".png")), width = 8, height = 5)
}
print(fread(file.path(D5, "5.1_Factor_associations_real.csv")), row.names = FALSE)
Zr <- fmat(fit$real$Z); zdt <- data.table(md[match(rownames(Zr), SampleID)], Zr)
zdt[, Time := factor(Time, levels = c("CTm02", "CTp02", "CTp06", "CTp10", "CTm10"))][, Light := factor(Light, levels = c("LL", "ML", "HL"))]
fl <- melt(zdt, id.vars = c("SampleID", "Light", "Time", "Replicate"), measure.vars = colnames(Zr), variable.name = "Factor", value.name = "Value")
p <- ggplot(fl, aes(Time, Value, colour = Light, group = Light)) + stat_summary(fun = mean, geom = "line", linewidth = 0.8) + geom_point(size = 1.4, alpha = 0.7) +
  facet_wrap(~Factor, scales = "free_y") + scale_colour_manual(values = light_cols) + labs(x = "Time (CTm02 = 2 h before light ... CTm10 = 2 h after dark)", y = "Factor value",
  title = "Factors of the real data across the day, by light level") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, D5, "5.2_Factors_over_time.png", 11, 7)
if (ncol(Zr) >= 2) {
  p <- ggplot(zdt, aes(Factor1, Factor2, colour = Light, shape = Time)) + geom_point(size = 2.8) + scale_colour_manual(values = light_cols) +
    labs(title = "Samples on the first two factors")
  save_plot(p, D5, "5.3_Factor1_vs_Factor2.png", 7, 5.5)
}


# Section 6. What drives each factor: weights and their enrichment ---------------
# The weight of a feature says how strongly, and in which direction, it follows the factor.
# Genes and proteins are ranked by their weight and analysed by GSEA with the eggNOG-mapper
# annotation of part 1. The real data use Chlamydomonas genome v6.1 identifiers (Cre01.g000050);
# the annotation of part 1 uses v5.5 genes (CHLRE_01g000050v5): most genes keep their number
# and are matched by renaming; new v6.1 genes (Cre..g80....) are left without annotation.
# Lipids are grouped by class (TG, DG, MGDG, DGDG, DGTS, PE, PG, PI, SQDG, PA).

D6 <- out_dir("6_Weights")
W <- fit$real$W; fcols <- setdiff(names(W), c("feature", "view"))
W <- merge(W, annot, by.x = "feature", by.y = "id", all.x = TRUE)
top <- rbindlist(lapply(fcols, function(f) W[order(-abs(get(f))), .SD[seq_len(min(TOP_N, .N))], by = view][, .(Factor = f, view, feature, symbol, description, weight = round(get(f), 3))]))
fwrite(top, file.path(D6, "6.1_Top_features_per_factor.csv"))
# Annotation (as in PC4_Enrichment.R)
read_emapper <- function(f) { x <- fread(file.path(DATA_DIR, f), skip = "#query", sep = "\t", quote = "", na.strings = c("-", ""), fill = TRUE); setnames(x, "#query", "query"); x[!startsWith(query, "##")] }
hdr <- grep("^>", readLines(file.path(DATA_DIR, CHLAMY_FASTA)), value = TRUE)
p2g <- data.table(query = sub("^>(\\S+).*", "\\1", hdr), id = sub(".*\\bgene:(\\S+).*", "\\1", hdr))
em <- merge(read_emapper(EGGNOG_CHLAMY), p2g, by = "query")[, id := sub("^CHLRE_(\\d+)g(\\d+)v5$", "Cre\\1.g\\2", id)]
split_col <- function(d, col) { d <- d[!is.na(get(col)), .(id, v = as.character(get(col)))]; d[, .(item = trimws(unlist(strsplit(v, ",")))), by = id][item != ""] }
go <- unique(split_col(em, "GOs"))
info <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(go$item), columns = c("ONTOLOGY", "TERM"), keytype = "GOID"))
go[, Ontology := info$ONTOLOGY[match(item, info$GOID)]][, Name := info$TERM[match(item, info$GOID)]]
kg <- unique(split_col(em, "KEGG_Pathway"))[grepl("^[0-9]{5}$", item), item := paste0("map", item)][startsWith(item, "map")]
kegg_names <- if (file.exists(file.path(DATA_DIR, KEGG_MAP))) unique(fread(file.path(DATA_DIR, KEGG_MAP))[, .(Term = Pathway, Name)]) else data.table(Term = character(), Name = character())
kg <- merge(kg[, .(Term = item, id)], kegg_names, by = "Term", all.x = TRUE)[, Name := ifelse(is.na(Name), Term, Name)]
sets <- list(GO_BP = unique(go[Ontology == "BP", .(Term = item, id, Name)]), KEGG = unique(kg))
cov <- W[view %in% c("mRNA", "protein"), .(Features = .N, Annotated = sum(feature %in% unique(c(sets$GO_BP$id, sets$KEGG$id)))), by = view]
fwrite(cov, file.path(D6, "6.2_Annotation_coverage.csv")); print(cov, row.names = FALSE)
gsea_w <- function(w, s) {
  w <- w[is.finite(w)]; pw <- lapply(split(s$id, s$Term), function(v) intersect(unique(v), names(w))); pw <- pw[lengths(pw) >= MIN_SIZE & lengths(pw) <= MAX_SIZE]
  if (!length(pw)) return(data.table())
  r <- suppressWarnings(fgsea::fgseaMultilevel(pw, sort(w, decreasing = TRUE), minSize = MIN_SIZE, maxSize = MAX_SIZE, eps = 0, BPPARAM = BiocParallel::SerialParam()))
  merge(as.data.table(r)[, .(Term = pathway, size, NES, padj)], unique(s[, .(Term, Name)]), by = "Term")
}
# GSEA of the weights of the factors that explain at least 2% of mRNA or protein variance
keep_f <- fit$real$R2[view %in% c("mRNA", "protein") & R2 >= 2, unique(factor)]
ws <- rbindlist(lapply(keep_f, function(f) rbindlist(lapply(c("mRNA", "protein"), function(v) {
  x <- W[view == v]; w <- setNames(x[[f]], x$feature)
  rbindlist(lapply(names(sets), function(a) gsea_w(w, sets[[a]])[, `:=`(Factor = f, view = v, Annotation = a)]))
}))))
if (nrow(ws)) {
  setcolorder(ws, c("Factor", "view", "Annotation", "Name", "size", "NES", "padj"))
  fwrite(ws[order(Factor, view, padj)], file.path(D6, "6.3_GSEA_weights.csv"))
  cat("GSEA of the weights: significant terms per factor and omics\n")
  print(dcast(ws[padj < FDR], Factor ~ view, fun.aggregate = length, value.var = "Term"), row.names = FALSE)
  best <- ws[padj < FDR][order(padj), head(.SD, 3), by = .(Factor, view)][, .(Factor, view, Name, NES = round(NES, 2), padj = signif(padj, 2))]
  print(best, row.names = FALSE)
}
# Lipid classes
lw <- W[view == "lipids"][, Class := sub("^([A-Za-z]+).*", "\\1", feature)]
lc <- melt(lw, id.vars = c("feature", "Class"), measure.vars = fcols, variable.name = "Factor", value.name = "weight")[, .(Species = .N, Mean_weight = round(mean(weight), 3)), by = .(Factor, Class)]
fwrite(lc, file.path(D6, "6.4_Lipid_classes_weights.csv"))
p <- ggplot(lc, aes(Factor, Class, fill = Mean_weight)) + geom_tile(colour = "white") + scale_fill_gradient2(low = "#2A6F97", mid = "white", high = "#B5541B", name = "Mean weight") +
  labs(x = NULL, y = NULL, title = "Lipid classes: mean weight on each factor") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, D6, "6.4_Lipid_classes_weights.png", 8, 4.5)


# Section 7. Session information ---------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
cat("\nDone. All results are in:", OUT, "\n")
