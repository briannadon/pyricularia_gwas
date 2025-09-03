# Start from the rubric image and add bioinformatics tools
FROM ubuntu:latest

# Switch to root to install packages
USER root

COPY ./submission/ /tmp/files/

# Install system dependencies first
RUN apt-get update && apt-get install -y wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install mambaforge using a specific version URL
RUN wget https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Mambaforge-24.3.0-0-Linux-x86_64.sh -O mambaforge.sh && \
    bash mambaforge.sh -b -p /opt/conda && \
    rm mambaforge.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Configure conda channels (mambaforge already has conda-forge)
RUN conda config --add channels bioconda && \
    conda config --set channel_priority strict

# Create separate environments for each tool to avoid version conflicts
# Create all environments in a single RUN command to better handle dependencies
RUN mamba create -n blast -y blast && \
    mamba create -n samtools -y samtools && \
    mamba create -n bwa -y bwa && \
    mamba create -n bowtie2 -y bowtie2 && \
    mamba create -n bedtools -y bedtools && \
    mamba create -n fastqc -y fastqc && \
    mamba create -n seqtk -y seqtk && \
    mamba create -n star -y -c bioconda star && \
    mamba create -n subread -y -c bioconda subread && \
    mamba create -n multiqc -y -c bioconda multiqc && \
    mamba clean --all -y

RUN  mamba create -n nextflow -y -c bioconda nextflow

# Create R environment with base R
# Install R base first, then add packages incrementally to avoid solver issues
RUN mamba create -n R -y -c conda-forge r-base=4.3.* && \
    mamba clean --all -y

# Install R packages in smaller groups to avoid dependency conflicts
RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             source /opt/conda/etc/profile.d/mamba.sh && \
             conda activate R && \
             mamba install -y -c conda-forge r-tidyverse r-devtools r-biocmanager zlib" && \
    mamba clean --all -y

# Create Python environment with common bioinformatics and data science packages
RUN mamba create -n python -y -c conda-forge python=3.11 && \
    mamba clean --all -y

# Install canonical Python packages for bioinformatics and data science
RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             source /opt/conda/etc/profile.d/mamba.sh && \
             conda activate python && \
             mamba install -y -c conda-forge \
                 numpy \
                 pandas \
                 scipy \
                 matplotlib \
                 seaborn \
                 scikit-learn \
                 jupyterlab \
                 ipython \
                 biopython \
                 pysam \
                 pyvcf3 \
                 requests \
                 plotly \
                 statsmodels \
                 networkx \
                 h5py" && \
    mamba clean --all -y


# Optional: Install Bioconductor packages (comment out if build time is too long)
RUN echo '#!/usr/bin/env Rscript\n\
if (!require("BiocManager", quietly = TRUE))\n\
     install.packages("BiocManager", repos="https://cloud.r-project.org")\n\
 BiocManager::install(version = "3.20", ask = FALSE, update = TRUE)\n\
# # Install commonly used Bioconductor packages\n\
 BiocManager::install(c(\n\
     "DESeq2",\n\
     "edgeR",\n\
     "limma",\n\
     "GenomicRanges",\n\
     "Biostrings"\n\
 ), ask = FALSE, update = FALSE)' > /tmp/install_bioc.R && \
     bash -c "source /opt/conda/etc/profile.d/conda.sh && \
              source /opt/conda/etc/profile.d/mamba.sh && \
              conda activate R && \
              Rscript /tmp/install_bioc.R" && \
     rm /tmp/install_bioc.R


# Create workspace directory if it doesn't exist
RUN mkdir -p /workspace

# Ensure the sandboxing user has access to conda
RUN if id -u 999 >/dev/null 2>&1; then \
        chown -R 999:999 /opt/conda && \
        chmod -R 755 /opt/conda; \
    fi

# Initialize conda/mamba for bash shell
RUN /opt/conda/bin/conda init bash && \
    /opt/conda/bin/mamba init bash && \
    echo "source /opt/conda/etc/profile.d/conda.sh" >> /etc/bash.bashrc && \
    echo "source /opt/conda/etc/profile.d/mamba.sh" >> /etc/bash.bashrc

# Switch back to the sandboxing user (if it exists in the base image)
# USER 999

# Ensure PATH includes conda for the user
ENV PATH="/opt/conda/bin:$PATH"

# Initialize conda/mamba for the user
RUN if [ -f /home/user/.bashrc ]; then \
        /opt/conda/bin/conda init bash && \
        /opt/conda/bin/mamba init bash; \
    fi

# The rubric image already has the startup command configured
