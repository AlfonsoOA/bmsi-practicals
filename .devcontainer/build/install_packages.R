# R packages of the practicals (PC2-PC4). CRAN packages come as binaries from the
# frozen Posit Package Manager snapshot of the rocker image; Bioconductor packages
# from the release given in BIOC_VERSION. The build stops if anything is missing.

options(Ncpus = max(1L, parallel::detectCores()))
cran <- c("BiocManager", "data.table", "ggplot2", "ggrepel", "pheatmap", "vegan", "ape", "UpSetR",
          "VennDiagram", "missForest", "imputeLCMD", "future", "future.apply", "Matrix", "readxl",
          "rstan", "BH", "StanHeaders", "RcppEigen", "RcppParallel", "pkgbuild", "zip", "rmarkdown")
bioc <- c("tximport", "DESeq2", "edgeR", "limma", "pcaMethods", "impute", "MSstats",
          "clusterProfiler", "enrichplot", "GO.db", "fgsea", "AnnotationDbi", "BiocParallel")

install.packages(cran)
BiocManager::install(version = Sys.getenv("BIOC_VERSION", "3.22"), ask = FALSE, update = FALSE)
BiocManager::install(bioc, ask = FALSE, update = FALSE)

missing <- setdiff(c(cran, bioc), rownames(installed.packages()))
if (length(missing) > 0) stop("Packages not installed: ", paste(missing, collapse = ", "))
for (p in c(cran, bioc)) suppressPackageStartupMessages(library(p, character.only = TRUE))
cat("All", length(c(cran, bioc)), "packages installed and loadable. R", as.character(getRversion()),
    "| Bioconductor", as.character(BiocManager::version()), "\n")
writeLines(capture.output(sessionInfo()), "/usr/local/share/bmsi_sessionInfo.txt")
