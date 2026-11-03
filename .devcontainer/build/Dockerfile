# Image for the Systems Molecular Biology computer practicals (UCO, 2026-27).
# Built by .github/workflows/build-image.yml and published to GHCR; students
# never build it: their codespace pulls the finished image.
#
# Contents: R 4.5 + CRAN/Bioconductor packages of PC2-PC4 (binary packages from
# Posit Public Package Manager, frozen at the date of the R release) and the
# conda environments of the terminal parts (PC1, PC2).
#
# Trinity and BUSCO are NOT installed (several GB, and Trinity needs far more than
# the 8 GB of a codespace): PC2_terminal.sh skips both when their precomputed
# results from ForStudents_PC2.zip are present. Their environments exist empty
# only so that the script does not try to create them.
FROM ghcr.io/rocker-org/devcontainer/r-ver:4.5

ARG BIOC_VERSION=3.22

# System tools: GNU time (memory checks), unzip, and libraries some packages need
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
 && apt-get install -y --no-install-recommends time unzip zip pigz \
      libglpk40 libxml2 libcurl4 libhdf5-dev libgsl27 \
 && rm -rf /var/lib/apt/lists/*

# R packages (one layer; fails the build if any package is missing).
# Packages that still ask for C++11 (MSstats) are compiled as C++17: current
# RcppArmadillo headers need at least C++14.
COPY install_packages.R /tmp/install_packages.R
RUN mkdir -p /root/.R && echo "CXX11STD = -std=gnu++17" > /root/.R/Makevars \
 && BIOC_VERSION=${BIOC_VERSION} Rscript /tmp/install_packages.R && rm -rf /tmp/* /root/.cache

# The Stan model of PC3 (section 4), compiled once here: compiling it in a codespace
# needs 2-3 GB of memory on top of the data and can exhaust the 8 GB. The workflow
# copies scripts/PC3_Proteomics.R into the build folder; the model code is taken from it.
COPY PC3_Proteomics.R /tmp/PC3_Proteomics.R
RUN mkdir -p /opt/bmsi/stan && Rscript -e ' \
      ex <- parse("/tmp/PC3_Proteomics.R"); \
      i <- which(vapply(ex, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("stan_code")), logical(1))); \
      stopifnot(length(i) == 1); eval(ex[[i]]); \
      f <- "/opt/bmsi/stan/protein_model.stan"; writeLines(stan_code, f); \
      invisible(rstan::stan_model(file = f, auto_write = TRUE)); stopifnot(file.exists(sub("stan$", "rds", f)))' \
 && chmod -R a+rX /opt/bmsi && rm -f /tmp/PC3_Proteomics.R

# Conda (Miniforge) and the environments of the terminal parts, same names and
# versions as in the practical scripts, so the scripts find them and skip installation
ENV CONDA_DIR=/opt/conda
RUN curl -fsSL -o /tmp/miniforge.sh \
      https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
 && bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && rm /tmp/miniforge.sh \
 && ${CONDA_DIR}/bin/conda config --system --set channel_priority strict \
 && ${CONDA_DIR}/bin/conda create -y -n bmsi_pc1 --override-channels -c conda-forge -c bioconda \
      python=3.11 fastqc=0.12.1 fastp=1.3.7 megahit=1.2.9 bwa=0.7.19 samtools=1.24 minimap2=2.31 \
      quast=5.3.0 seqkit=2.14.0 multiqc=1.35 sra-tools=3.4.1 blast=2.16.0 \
 && ${CONDA_DIR}/bin/conda create -y -n bmsi_pc2 --override-channels -c conda-forge -c bioconda \
      python=3.11 fastqc=0.12.1 fastp=1.3.7 multiqc=1.35 salmon=2.8.0 minimap2=2.31 seqkit=2.14.0 sra-tools=3.4.1 \
 && ${CONDA_DIR}/bin/conda create -y -n bmsi_pc2_trinity \
 && ${CONDA_DIR}/bin/conda create -y -n bmsi_pc2_busco \
 && ${CONDA_DIR}/bin/conda clean -afy \
 && chown -R 1000:1000 ${CONDA_DIR}

# conda available in every terminal of the user (bash), base not activated
RUN USERNAME=$(getent passwd 1000 | cut -d: -f1) \
 && su - ${USERNAME} -c "${CONDA_DIR}/bin/conda init bash && ${CONDA_DIR}/bin/conda config --set auto_activate_base false"
