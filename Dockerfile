# Use an official Python base image
ARG BASE_IMAGE=3.12-slim
FROM python:${BASE_IMAGE} AS base

# Set environment variables
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH="/opt/conda/bin:$PATH" \
    POETRY_HOME="/opt/poetry" \
    POETRY_VIRTUALENVS_CREATE=false

# Install required tools and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl bzip2 build-essential \
    git file pkg-config swig libcairo2-dev libnetpbm10-dev  \
    netpbm libpng-dev libjpeg-dev zlib1g-dev  \
    libbz2-dev libcfitsio-dev wcslib-dev && \
    rm -rf /var/lib/apt/lists/*

# Create a layer for Miniconda
FROM base AS miniconda
WORKDIR /tmp

# Download and install Miniconda
ARG MINICONDA_VERSION=py312_24.9.2-0
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    conda clean -afy

RUN conda install -y -c conda-forge astromatic-source-extractor astromatic-scamp astromatic-swarp astromatic-psfex astrometry gsl wcstools

# Create a layer for Poetry
FROM miniconda AS poetry

# Install Poetry
RUN curl -sSL https://install.python-poetry.org | python3 - && \
    ln -s $POETRY_HOME/bin/poetry /usr/local/bin/poetry

# Configure Poetry cache to use a shared volume
ENV POETRY_CACHE_DIR=/cache/poetry
VOLUME /cache/poetry

# Create a layer for rust
FROM poetry AS rust

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

ENV PATH="/root/.cargo/bin:${PATH}"

# Create a final layer for the application
FROM rust AS install
WORKDIR /app

# Copy only dependency files first to leverage Docker cache
COPY pyproject.toml poetry.lock ./

# Install dependencies
RUN poetry install --no-root

# Copy the rest of the application code
COPY . ./

RUN poetry install

FROM install AS final

ENV RAW_DATA_DIR=/data
ENV OUTPUT_DATA_DIR=/data

# Set the entry point for the container
ENTRYPOINT ["mirar-run"]
