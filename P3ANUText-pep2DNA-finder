import json
import csv
import os
from collections import Counter
import tkinter as tk
from tkinter import filedialog

def select_json_file():
    #Opens a GUI file dialog to select the input JSON file.
    root = tk.Tk()
    root.withdraw() # Hide the main root window
    
    file_path = filedialog.askopenfilename(
        title="Select NGS Data JSON File",
        filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")]
    )
    return file_path

def get_peptide_input():
    #Prompts user in console for the target peptide sequence.
    while True:
        peptide = input("Enter the peptide sequence of interest (e.g., SHSSCHHR): ").strip().upper()
        if peptide:
            return peptide
        print("Sequence cannot be empty. Please try again.")

def find_dna_subsequences(json_data, target_peptide):
    #Finds instances of the target peptide in the protein sequences, 
    #calculates the corresponding codon offsets, and extracts the exact DNA sub-sequences.
    extracted_dna_sequences = []
    
    for entry_id, entry_data in json_data.items():
        protein_seq = entry_data.get("proteinSequence", "")
        dna_seq = entry_data.get("sequences", "")
        
        # Check if the target peptide exists in this entry
        if target_peptide in protein_seq:
            start_aa_idx = protein_seq.find(target_peptide)
            end_aa_idx = start_aa_idx + len(target_peptide)
            
            # Map amino acid indices to DNA codon indices (1 AA = 3 nucleotides)
            start_dna_idx = start_aa_idx * 3
            end_dna_idx = end_aa_idx * 3
            
            # Ensure the indices are valid within the DNA sequence length
            if end_dna_idx <= len(dna_seq):
                target_dna = dna_seq[start_dna_idx:end_dna_idx]
                extracted_dna_sequences.append(target_dna)
                
    return extracted_dna_sequences

def write_results_to_csv(dna_counts, target_peptide):
    #Outputs the DNA sequence variants and their absolute counts to a CSV file.
    output_filename = f"dna_variants_{target_peptide}.csv"
    
    # Sort variants by frequency (highest first)
    sorted_variants = dna_counts.most_common()
    
    try:
        with open(output_filename, mode='w', newline='', encoding='utf-8') as csv_file:
            writer = csv.writer(csv_file)
            # Write header
            writer.writerow(["Target Peptide", "DNA Variant Sequence", "Frequency Count"])
            
            # Write data rows
            for dna_seq, count in sorted_variants:
                writer.writerow([target_peptide, dna_seq, count])
                
        print(f"\nSuccess! Results written to: {os.path.abspath(output_filename)}")
        print(f"Found {len(dna_counts)} unique DNA combination(s) out of {sum(dna_counts.values())} total matches.")
    except Exception as e:
        print(f"Error writing CSV file: {e}")

def main():
    print("NGS Peptide-to-DNA Variant Extractor")
    
    #Select the file
    print("Please select your JSON data file using the popup window...")
    json_path = select_json_file()
    if not json_path:
        print("Selection cancelled. Exiting program.")
        return
        
    #Load the JSON data
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            json_data = json.load(f)
    except Exception as e:
        print(f"Failed to parse JSON file: {e}")
        return

    #Get user target peptide sequence
    target_peptide = get_peptide_input()
    
    #Extract corresponding DNA sub-sequences
    matching_dna = find_dna_subsequences(json_data, target_peptide)
    
    if not matching_dna:
        print(f"\nThe peptide sequence '{target_peptide}' was not found in any dataset records.")
        return
        
    #Count distinct combinations
    dna_counts = Counter(matching_dna)
    
    #Generate CSV Report
    write_results_to_csv(dna_counts, target_peptide)

if __name__ == "__main__":
    main()
