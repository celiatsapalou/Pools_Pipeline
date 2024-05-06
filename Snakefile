rule all:
    input:
        expand("results/{sample}/{sample}.{cell}.bam", sample=["sample1", "sample2"], cell=["cell1", "cell2"])

# Create initial mock FASTQ files
rule create_mock_fastq:
    output:
        fastq="mock_data/{sample}.{cell}.fastq"
    shell:
        "echo 'Fake sequence data for {wildcards.sample} {wildcards.cell}' > {output}"

# Convert FASTQ to BAM - store intermediate BAM files separately
rule fastq_to_bam:
    input:
        fastq="mock_data/{sample}.{cell}.fastq"
    output:
        bam="intermediate/{sample}.{cell}.bam"
    shell:
        "echo 'Converting {input.fastq} to BAM' && touch {output}"

# Generate demultiplexing table
rule demultiplexing:
    output:
        table="demultiplex_table.csv"
    shell:
        """
        echo 'cell,Pool,1KG_identified_sample,z-score_value,SNP_nb,Trustable' > {output}
        echo 'cell1,pool1,sample1,5.95255829,19,TRUE' >> {output}
        echo 'cell2,pool1,sample2,5.986672013,15,TRUE' >> {output}
        """

# Checkpoint for parsing demultiplexing table
checkpoint parse_demultiplex_table:
    input:
        "demultiplex_table.csv"
    output:
        "parsed_table.txt"
    script:
        "scripts/parse_table.py"

# Function to gather outputs based on checkpoint
def gather_samples(wildcards):
    checkpoint_output = checkpoints.parse_demultiplex_table.get(**wildcards).output[0]
    import pandas as pd
    table = pd.read_csv(checkpoint_output, sep=',')
    return expand("intermediate/{{sample}}.{{cell}}.bam", zip, sample=table['1KG_identified_sample'], cell=table['cell'])

# Organize samples into final directories
rule organize_samples:
    input:
        gather_samples
    output:
        bam="results/{sample}/{sample}.{cell}.bam"
    shell:
        """
        mkdir -p results/{wildcards.sample}
        mv {input} {output}
        """
