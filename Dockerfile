# syntax=docker/dockerfile:1

FROM rocker/r-ver:4.6.1

# ------------------------------------------------------------------------------
# Environment
# ------------------------------------------------------------------------------

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    USE_BUNDLED_LIBUV=1 \
    RETICULATE_PYTHON="/root/.virtualenvs/pamchip-env/bin/python"

# ------------------------------------------------------------------------------
# 1. System Dependencies
# ------------------------------------------------------------------------------

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    unzip \
    git \
    patch \
    build-essential \
    cmake \
    pkg-config \
    python3 \
    python3-pip \
    python3-venv \
    libcairo2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libjpeg-dev \
    libtiff5-dev \
    libpng-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------------
# 2. Install Quarto
# ------------------------------------------------------------------------------

ARG QUARTO_VERSION="1.10.18"

RUN ARCH="$(dpkg --print-architecture)" && \
    case "$ARCH" in \
      amd64) QUARTO_ARCH="amd64" ;; \
      arm64) QUARTO_ARCH="arm64" ;; \
      *) echo "Unsupported architecture: $ARCH" && exit 1 ;; \
    esac && \
    curl \
      --retry 5 \
      --retry-connrefused \
      --retry-delay 2 \
      -L \
      -o /tmp/quarto.deb \
      "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${QUARTO_ARCH}.deb" && \
    dpkg -i /tmp/quarto.deb && \
    rm -f /tmp/quarto.deb

# ------------------------------------------------------------------------------
# 3. Install CRAN Package Dependencies
# ------------------------------------------------------------------------------

RUN Rscript - <<'RS'
options(
  repos = c(
    CRAN = "https://packagemanager.posit.co/cran/__linux__/bookworm/latest"
  )
)

pkgs <- c(
  "targets",
  "tarchetypes",
  "yaml",
  "dplyr",
  "purrr",
  "readr",
  "stringr",
  "tidyr",
  "httr2",
  "jsonlite",
  "reticulate",
  "quarto",
  "knitr",
  "rmarkdown",
  "kableExtra",
  "ggplot2",
  "here",
  "BiocManager"
)

install.packages(
  pkgs,
  Ncpus = 4
)

missing <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0L) {
  stop(
    "Failed to install CRAN packages: ",
    paste(missing, collapse = ", ")
  )
}
RS

# ------------------------------------------------------------------------------
# 4. Install Bioconductor Reporting / Enrichment Dependencies
# ------------------------------------------------------------------------------

RUN Rscript - <<'RS'
options(
  repos = BiocManager::repositories(
    site_repository = "https://packagemanager.posit.co/cran/__linux__/bookworm/latest"
  )
)

pkgs <- c(
  "ReactomePA",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot"
)

BiocManager::install(
  pkgs,
  ask = FALSE,
  update = FALSE,
  Ncpus = 4
)

missing <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing) > 0L) {
  stop(
    "Failed to install Bioconductor packages: ",
    paste(missing, collapse = ", ")
  )
}
RS

# ------------------------------------------------------------------------------
# 5. Verify R Reporting Environment
# ------------------------------------------------------------------------------

RUN Rscript - <<'RS'
required <- c(
  "targets",
  "tarchetypes",
  "yaml",
  "dplyr",
  "purrr",
  "readr",
  "stringr",
  "tidyr",
  "httr2",
  "jsonlite",
  "reticulate",
  "quarto",
  "knitr",
  "rmarkdown",
  "kableExtra",
  "ggplot2",
  "here",
  "BiocManager",
  "ReactomePA",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot"
)

cat("R:", R.version.string, "\n\n")

for (pkg in required) {

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required R package could not be loaded: ",
      pkg
    )
  }

  cat(
    pkg,
    ": ",
    as.character(
      utils::packageVersion(pkg)
    ),
    "\n",
    sep = ""
  )
}
RS

# ------------------------------------------------------------------------------
# 6. Create Reticulate Python Environment
# ------------------------------------------------------------------------------

RUN Rscript - <<'RS'
py_path <- reticulate::install_python(
  version = "3.10.14"
)

reticulate::virtualenv_create(
  "pamchip-env",
  python = py_path
)
RS

# ------------------------------------------------------------------------------
# 7. Install Python Dependencies for Kinase Library + GPS6
# ------------------------------------------------------------------------------

RUN /root/.virtualenvs/pamchip-env/bin/python -m pip install \
      --upgrade \
      pip \
      setuptools \
      wheel && \
    /root/.virtualenvs/pamchip-env/bin/python -m pip install \
      --no-cache-dir \
      "numpy<2" \
      pandas \
      joblib \
      scikit-learn \
      lightgbm \
      h5py \
      tensorflow \
      requests \
      kinase-library

# ------------------------------------------------------------------------------
# 8. Verify Python Environment During Image Build
# ------------------------------------------------------------------------------

RUN /root/.virtualenvs/pamchip-env/bin/python - <<'PY'
import sys

import numpy
import pandas
import joblib
import sklearn
import lightgbm
import h5py
import tensorflow as tf
import requests
import kinase_library

print("Python:", sys.version)
print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("joblib:", joblib.__version__)
print("sklearn:", sklearn.__version__)
print("lightgbm:", lightgbm.__version__)
print("h5py:", h5py.__version__)
print("tensorflow:", tf.__version__)
print("requests:", requests.__version__)
print("kinase_library: imported successfully")
PY

# ------------------------------------------------------------------------------
# 9. Verify Quarto
# ------------------------------------------------------------------------------

RUN quarto --version

# ------------------------------------------------------------------------------
# 10. Project Directory Setup
# ------------------------------------------------------------------------------

WORKDIR /project

COPY . /project

RUN mkdir -p \
    data/external/gps6 \
    data/external/kinase_library \
    data/external/signor \
    data/external/kinase_taxonomy \
    data/raw \
    results

# ------------------------------------------------------------------------------
# 11. Default Command
# ------------------------------------------------------------------------------

CMD ["Rscript", "-e", "targets::tar_make()"]