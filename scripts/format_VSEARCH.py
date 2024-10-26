import pandas as pd
from tqdm import tqdm

def load_metadata(metadata_file):
    metadata_df = pd.read_csv(metadata_file, sep='\t', low_memory=False)
    metadata_df = metadata_df.set_index('COInrID')
    return metadata_df

def load_sequences(sequence_file):
    sequence_df = pd.read_csv(sequence_file, sep='\t', header=None, names=['COInrID', 'seqID', 'sequence'])
    sequence_df = sequence_df.set_index('COInrID')
    return sequence_df

def generate_output(metadata_df, sequence_df, output_file):
    with open(output_file, 'w') as out_file:
        # Use tqdm to add a progress bar for the sequence processing loop
        for coid, sequence_row in tqdm(sequence_df.iterrows(), total=len(sequence_df), desc="Processing sequences"):
            if coid in metadata_df.index:
                meta_row = metadata_df.loc[coid]
                
                # Extract taxonomic information and check for missing high-level taxonomy
                kingdom = meta_row.get('kingdom')
                phylum = meta_row.get('phylum')
                class_ = meta_row.get('class')
                family = meta_row.get('family') if pd.notna(meta_row.get('family')) else meta_row.get('subfamily')
                
                # Skip record if any high-level taxonomy is missing
                if pd.isna(kingdom) or pd.isna(phylum) or pd.isna(class_) or pd.isna(family):
                    continue
                
                # Optional taxonomic levels
                genus = meta_row.get('genus')
                species = meta_row.get('species')
                
                sequence = sequence_row['sequence']
                
                # Construct the header line conditionally based on available taxonomy
                header = f">{coid};tax=k:{kingdom},d:{phylum},c:{class_},f:{family}"
                if pd.notna(genus):
                    header += f",g:{genus}"
                if pd.notna(species):
                    header += f",s:{species}"
                
                # Write the formatted record to the output file
                out_file.write(f"{header}\n{sequence}\n")

def main(metadata_file, sequence_file, output_file):
    # Load the metadata and sequence files
    metadata_df = load_metadata(metadata_file)
    sequence_df = load_sequences(sequence_file)
    
    # Generate the output file with a progress bar
    generate_output(metadata_df, sequence_df, output_file)
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    sequence_file = 'COInr.tsv'
    metadata_file = 'COInr_metadata.tsv'
    output_file = 'COI.COInr.24-10-26.fasta'

    main(metadata_file, sequence_file, output_file)
