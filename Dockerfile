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
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libjpeg-dev \
    libtiff5-dev \
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
      -o quarto.deb \
      "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${QUARTO_ARCH}.deb" && \
    dpkg -i quarto.deb && \
    rm quarto.deb

# ------------------------------------------------------------------------------
# 3. Install R Package Dependencies
# ------------------------------------------------------------------------------

RUN Rscript -e ' \
  options( \
    repos = c( \
      CRAN = "https://packagemanager.posit.co/cran/__linux__/bookworm/latest" \
    ) \
  ); \
  pkgs <- c( \
    "targets", \
    "tarchetypes", \
    "yaml", \
    "dplyr", \
    "purrr", \
    "readr", \
    "stringr", \
    "tidyr", \
    "httr2", \
    "jsonlite", \
    "reticulate", \
    "quarto", \
    "knitr", \
    "rmarkdown", \
    "kableExtra", \
    "ggplot2", \
    "here" \
  ); \
  install.packages(pkgs, Ncpus = 4); \
  missing <- setdiff( \
    pkgs, \
    rownames(installed.packages()) \
  ); \
  if (length(missing) > 0) { \
    stop( \
      "Failed to install packages: ", \
      paste(missing, collapse = ", ") \
    ); \
  }'

# ------------------------------------------------------------------------------
# 4. Create Reticulate Python Environment
# ------------------------------------------------------------------------------

# Install a dedicated Python build through reticulate so the environment is
# compatible with R/reticulate and can also be called directly by GPS6.
RUN Rscript -e ' \
  py_path <- reticulate::install_python( \
    version = "3.10.14" \
  ); \
  reticulate::virtualenv_create( \
    "pamchip-env", \
    python = py_path \
  ) \
'

# ------------------------------------------------------------------------------
# 5. Install Python Dependencies for Kinase Library + GPS6
# ------------------------------------------------------------------------------

# The same virtualenv is used for:
#
#   - The Kinase Library through reticulate
#   - GPS6 through system2()
#
# GPS6 requires:
#
#   numpy
#   pandas
#   joblib
#   scikit-learn
#   lightgbm
#   h5py
#   tensorflow
#
# The Kinase Library additionally requires:
#
#   kinase-library
#   requests
#
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
# 6. Verify Python Environment During Image Build
# ------------------------------------------------------------------------------

# Fail the Docker build immediately if any required Python package cannot load.
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

print("Python:", sys.version)
print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("joblib:", joblib.__version__)
print("sklearn:", sklearn.__version__)
print("lightgbm:", lightgbm.__version__)
print("h5py:", h5py.__version__)
print("tensorflow:", tf.__version__)
print("requests:", requests.__version__)

try:
    import kinase_library
    print("kinase_library: imported successfully")
except ImportError:
    import kinase_library as kl
    print("kinase_library: imported successfully")
PY

# ------------------------------------------------------------------------------
# 7. Project Directory Setup
# ------------------------------------------------------------------------------

WORKDIR /project

COPY . /project

RUN mkdir -p \
    data/external/gps6 \
    data/external/kinase_library \
    data/external/signor \
    data/raw \
    results

# ------------------------------------------------------------------------------
# 8. Default Command
# ------------------------------------------------------------------------------

CMD ["Rscript", "-e", "targets::tar_make()"]
