import pandas as pd
import random


input_structure = {
    "poolA": ["cell1A", "cell2A", "cell3A", "cell4A"],
    "poolB": ["cell1B", "cell2B", "cell3B", "cell4B"],
}
pool_composition = {
    "poolA": ["sampleV", "sampleW"],
    "poolB": ["sampleX", "sampleY"],
}


rule all:
    input:
        [
            expand(
                "test-run/{pool}/arbigent/{sample}.txt",
                pool=pool,
                sample=pool_composition[pool],
            )
            for pool in pool_composition.keys()
        ],


# Create initial mock FASTQ files
rule create_mock_fastq:
    output:
        fastq="test-run/{pool}/fastq/{cell}.fastq",
    shell:
        "echo 'Fake sequence data for {wildcards.pool} {wildcards.cell}' > {output}"


# Convert FASTQ to BAM - store intermediate BAM files separately
rule fastq_to_bam:
    input:
        fastq="test-run/{pool}/fastq/{cell}.fastq",
    output:
        bam="test-run/{pool}/bam/{cell}.bam",
    shell:
        "echo 'Converting {input.fastq} to BAM' && touch {output}"


# Generate demultiplexing table
checkpoint demultiplexing:
    input:
        lambda wc: expand(
            "test-run/{pool}/bam/{cell}.bam",
            pool=wc.pool,
            cell=input_structure[wc.pool],
        ),
    output:
        table="test-run/{pool}/demultiplexing/demultiplexing_table.csv",
    run:
        random.seed(3)  # for reproducibility

        # Prepare the output DataFrame list
        data = []

        # Iterate over each cell in the pool
        for cell in input_structure[wildcards.pool]:
            # Select a random sample from the corresponding pool
            sample = random.choice(pool_composition[wildcards.pool])

            # Generate a random z-score value, normally distributed between 5 and 6, with 4 decimals
            z_score = round(random.uniform(5, 6), 4)

            # Generate a random randint SNP_nb between 50 and 100
            SNP_nb = random.randint(50, 100)

            # Prepare a dictionary for the row
            row = {
                "cell": cell,
                "Pool": wildcards.pool,
                "1KG_identified_sample": sample,
                "z-score_value": z_score,
                "SNP_nb": SNP_nb,
                "Trustable": "TRUE",
            }

            # Append row to the data list
            data.append(row)

            # Create a DataFrame from the list of dictionaries
        df = pd.DataFrame(data)

        # Write the DataFrame to file with append mode off (write all at once)
        df.to_csv(output.table, mode="w", index=False, header=True)


# # Function to gather outputs based on checkpoint
# def gather_samples(wildcards):
#     # checkpoint_output = checkpoints.parse_demultiplex_table.get(**wildcards).output[0]
#     checkpoint_output = checkpoints.demultiplexing.get(**wildcards).output.table
#     import pandas as pd

#     table = pd.read_csv(checkpoint_output, sep=",")
#     table = table.loc[
#         (table["Pool"] == wildcards.pool)
#         & (table["1KG_identified_sample"] == wildcards.sample)
#     ]

#     return expand(
#         "test-run/{pool}/bam/{cell}.bam",
#         pool=wildcards.pool,
#         cell=table["cell"].unique().tolist(),
#     )


# Organize samples into final directories
rule organize_samples:
    input:
        bam="test-run/{pool}/bam/{cell}.bam",
    output:
        bam_samplewise="test-run/{pool}/bam_samplewise/{sample}/{cell}.bam",
    params:
        bam_samplewise_dir=lambda wc, input, output: "/".join(output.bam_samplewise.split("/")[:-1]),
    shell:
        """
        echo {params.bam_samplewise_dir}
        echo {output.bam_samplewise}
        echo {input}
        mkdir -p {params.bam_samplewise_dir}
        cp {input} {output.bam_samplewise}
        """


def return_cells_samplewise(wildcards):
    checkpoint_output = checkpoints.demultiplexing.get(**wildcards).output.table
    table = pd.read_csv(checkpoint_output, sep=",")
    table = table.loc[
        (table["Pool"] == wildcards.pool)
        & (table["1KG_identified_sample"] == wildcards.sample)
    ]

    return expand(
        "test-run/{pool}/bam_samplewise/{sample}/{cell}.bam",
        pool=wildcards.pool,
        sample=wildcards.sample,
        cell=table["cell"].unique().tolist(),
    )


rule arbigent:
    input:
        return_cells_samplewise,
        demultiplexing="test-run/{pool}/demultiplexing/demultiplexing_table.csv"
    output:
        "test-run/{pool}/arbigent/{sample}.txt",
    shell:
        # using awk, filter the demultiplexing table to only include the sample of interest
        "awk -F',' '{{if ($3 == \"{wildcards.sample}\") print $0}}' {input.demultiplexing} > {output}"
