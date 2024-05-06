# Use an official Miniconda image as a parent image
FROM condaforge/mambaforge

# Set the working directory in the container
WORKDIR /workspace

# Create the Conda environment using mamba for faster installs
RUN mamba create --name snakemake-env -c bioconda -c conda-forge snakemake zsh -y && \
    mamba init zsh && \
    mamba clean --all -y && \
    echo "source /workspace/mambaforge/bin/activate" >> ~/.bashrc && \
    echo "conda activate snakemake-env" >>  ~/.bashrc
