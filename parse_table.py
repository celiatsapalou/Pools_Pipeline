import pandas as pd

# Load the demultiplex table using a comma as the delimiter
df = pd.read_csv(snakemake.input[0], delimiter=',')

# Remove any surrounding whitespace from column names and values
df.columns = df.columns.str.strip()
df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

# Print columns to check correct parsing
print("Columns found:", df.columns)

# Print the first few rows to verify data integrity
print(df.head())

# Generate output paths
df['output_path'] = df.apply(lambda row: f"results/{row['1KG_identified_sample']}/{row['1KG_identified_sample']}.{row['cell']}.bam", axis=1)

# Save the output paths to the specified output file
df['output_path'].to_csv(snakemake.output[0], header=False, index=False)
