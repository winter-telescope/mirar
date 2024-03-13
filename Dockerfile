FROM debian:latest
WORKDIR /usr/src/mirar

RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY . .

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH

#
#
#RUN /usr/local/bin/python -m pip install --upgrade pip
#RUN pip install -e .
#RUN ln -s /usr/bin/SWarp /usr/bin/swarp
#RUN ln -s /usr/bin/source-extractor /usr/bin/sex
