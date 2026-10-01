# Syntax directive fixes the BuildKit platform warning
# syntax=docker/dockerfile:1

ARG TARGETPLATFORM=linux/amd64
FROM rocker/r-ver:4.6.1

# Set environment variables for non-interactive installations
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    USE_BUNDLED_LIBUV=1

# ------------------------------------------------------------------------------
# 1. System Dependencies & Python Compilation Build Tools
# ------------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    unzip \
    git \
    patch \
    build-essential \
    libffi-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    liblzma-dev \
    libgl1 \
    libglx-mesa0 \
    libglib2.0-0 \
    libcurl4-openssl-dev \
    libxml2-dev \
    libuv1-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------------
# 2. Install Quarto via .deb Package with Retries
# ------------------------------------------------------------------------------
ARG QUARTO_VERSION="1.4.550"
RUN curl --retry 5 --retry-connrefused --retry-delay 2 -L -o quarto.deb \
    "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
    && dpkg -i quarto.deb \
    && rm quarto.deb

# ------------------------------------------------------------------------------
# 3. Download and Install GPS 6.0 Standalone (Linux)
# ------------------------------------------------------------------------------
WORKDIR /opt/gps6
RUN curl --retry 5 --retry-connrefused --retry-delay 2 -O \
    https://gps.biocuckoo.cn/down/GPS_6.0_unix_20230422.tar.gz \
    && tar -zxvf GPS_6.0_unix_20230422.tar.gz \
    && rm GPS_6.0_unix_20230422.tar.gz \
    && chmod +x /opt/gps6/* || true

ENV PATH="/opt/gps6:$PATH"

# ------------------------------------------------------------------------------
# 4. Install R Package Dependencies & Reticulate Python Virtualenv
# ------------------------------------------------------------------------------
RUN Rscript -e ' \
  options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/bookworm/latest")); \
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
  installed <- install.packages(pkgs, Ncpus = 4); \
  missing <- setdiff(pkgs, rownames(installed.packages())); \
  if (length(missing) > 0) { \
    stop("Failed to install packages: ", paste(missing, collapse = ", ")); \
  }'

# Setup --enable-shared compatible Python via reticulate
ENV RETICULATE_PYTHON="/root/.virtualenvs/pamchip-env/bin/python"
RUN Rscript -e ' \
  py_path <- reticulate::install_python(version = "3.10.14"); \
  reticulate::virtualenv_create("pamchip-env", python = py_path); \
  reticulate::virtualenv_install("pamchip-env", packages = c("kinase-library", "requests", "pandas")) \
'

# ------------------------------------------------------------------------------
# 5. Project Directory Setup
# ------------------------------------------------------------------------------
WORKDIR /project

COPY . /project

RUN mkdir -p data/external/gps6 data/external/kinase_library data/external/signor data/raw

CMD ["Rscript", "-e", "targets::tar_make()"]
