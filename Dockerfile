# Use an official Python base image
ARG BASE_IMAGE=3.12-slim
FROM python:${BASE_IMAGE} as base

# Set environment variables
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH="/opt/conda/bin:$PATH" \
    POETRY_HOME="/opt/poetry" \
    POETRY_VIRTUALENVS_CREATE=false

# Install required tools and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl bzip2 build-essential && \
    rm -rf /var/lib/apt/lists/*

# Create a layer for Miniconda
FROM base as miniconda
WORKDIR /tmp

# Download and install Miniconda
ARG MINICONDA_VERSION=py312_24.9.2-0
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    conda clean -afy

RUN conda install -c conda-forge astromatic-source-extractor astromatic-scamp astromatic-swarp astromatic-psfex astrometry gsl wcstools

# Create a layer for Poetry
FROM miniconda as poetry

# Install Poetry
RUN curl -sSL https://install.python-poetry.org | python3 - && \
    ln -s $POETRY_HOME/bin/poetry /usr/local/bin/poetry

# Configure Poetry cache to use a shared volume
ENV POETRY_CACHE_DIR=/cache/poetry
VOLUME /cache/poetry

# Create a layer for rust
FROM poetry as rust

RUN echo "here"

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

ENV PATH="/root/.cargo/bin:${PATH}"

## Create a layer for postgres
#FROM rust as postgres
#
## Set up database
#RUN apt-get update
#RUN apt-get install -y postgresql postgresql-common
##RUN apt-get install -y postgresql-server-dev-14
#RUN #service postgresql status
#RUN ls -a /var/run/postgresql/
#RUN postgres
#RUN service postgresql start
#RUN service postgresql status
#RUN psql -U postgres -c "CREATE USER runner WITH PASSWORD 'runner_password'; GRANT ALL PRIVILEGES ON DATABASE postgres TO runner; ALTER USER runner WITH SUPERUSER;"
#
## Create a layer for q3c
#FROM postgres as q3c
#
#RUN git clone https://github.com/segasai/q3c.git
#RUN make -C q3c
#RUN make -C q3c install

# Create a final layer for the application
FROM rust as final
WORKDIR /app

# Copy only dependency files first to leverage Docker cache
COPY pyproject.toml poetry.lock ./

# Install dependencies
RUN poetry install --no-root

# Copy the rest of the application code
COPY . .

# Set the entry point for the container
ENTRYPOINT ["bash"]
