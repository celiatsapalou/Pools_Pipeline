# Use an official Miniconda image as a parent image
FROM condaforge/mambaforge

# Set the working directory in the container
WORKDIR /workspace

# Install zsh
RUN mamba install -c conda-forge zsh -y

# Change default shell to zsh
SHELL ["/bin/zsh", "-c"]

# Install mamba from the conda-forge channel
RUN conda install -c conda-forge mamba

# Create the Conda environment using mamba for faster installs
RUN mamba create --name snakemake-env python=3.9 -y

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "snakemake-env", "/bin/bash", "-c"]

# Install Snakemake using mamba for better performance
RUN mamba install -c bioconda -c conda-forge snakemake -y

# Ensure mamba is activated upon start
RUN echo "source activate snakemake-env" >> ~/.bashrc
