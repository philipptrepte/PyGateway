#!/usr/bin/env python3

import argparse
import csv
import os
import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqIO.XdnaIO import XdnaIterator
from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord

def read_xdna_file(xdna_file):
    """Reads a .xdna file and transforms all sequences to upper case."""
    with open(xdna_file, "rb") as file:
        sequence = next(XdnaIterator(file))
        sequence.seq = sequence.seq.upper()
        name = str(xdna_file)
        return sequence, name
    
def read_dna_file(dna_file):
    """Reads a .dna file and transforms all sequences to upper case."""
    data = snapgene_file_to_dict(filepath=dna_file)
    sequence = Seq(data["seq"])
    sequence.seq = sequence.upper()
    name = str(dna_file)
    return sequence, name

def read_gb_file(gb_file):
    """Reads a .gb file and transforms all sequences to upper case."""
    with open(gb_file, "r") as file:
        sequence = SeqIO.read(file, "genbank")
        sequence.seq = sequence.seq.upper()
        name = str(gb_file)
        return sequence, name

def create_expression_vector(entry_vector_file, destination_vector_file, output_folder, features_file):
    print(f"Creating expression vector from {entry_vector_file} and {destination_vector_file}")
    if entry_vector_file.endswith(".xdna"):
        entry_vector, entry_vector_name = read_xdna_file(entry_vector_file)
    elif entry_vector_file.endswith(".dna"):
        entry_vector, entry_vector_name = read_dna_file(entry_vector_file)
    elif entry_vector_file.endswith(".gb"):
        entry_vector, entry_vector_name = read_gb_file(entry_vector_file)
    else:
        raise ValueError("Invalid file extension for entry vector. Only .xdna, .dna, and .gb files are supported.")

    # Read the entry and destination vectors
    entry_vector_name = entry_vector_name.split("/")[-1]
    entry_name = entry_vector_name.split("pDONR")[-1]
    entry_name = entry_vector_name.split("pdonr")[-1]
    entry_name = entry_vector_name.split("PDONR")[-1]
    entry_name = entry_name.split("pENTR")[-1]
    entry_name = entry_name.split("pentr")[-1]
    entry_name = entry_name.split("PENTR")[-1]
    entry_name = entry_name.split("201")[-1]
    entry_name = entry_name.split("221")[-1]
    entry_name = entry_name.split("223")[-1]
    entry_name = entry_name.split("_")[1]
    entry_name = entry_name.split(".xdna")[0]
    entry_name = entry_name.split(".dna")[0]
    entry_name = entry_name.split(".gb")[0]

    if destination_vector_file.endswith(".xdna"):
        destination_vector, destination_name = read_xdna_file(destination_vector_file)
        destination_name = destination_name.split("/")[-1]
        destination_name = destination_name.split(".xdna")[0]
    elif destination_vector_file.endswith(".dna"):
        destination_vector, destination_name = read_dna_file(destination_vector_file)
        destination_name = destination_name.split("/")[-1]
        destination_name = destination_name.split(".dna")[0]
    elif destination_vector_file.endswith(".gb"):
        destination_vector, destination_name = read_gb_file(destination_vector_file)
        destination_name = destination_name.split("/")[-1]
        destination_name = destination_name.split(".gb")[0]
    else:
        raise ValueError("Invalid file extension for destination vector. Only .xdna, .dna, and .gb files are supported.")
        
    new_name = destination_name.replace("GW", entry_name)
    
    gateway_dict = {"att1_shared" : "TTTGTACAAAAAAG", 
                    "att2_shared" : "CTTTCTTGTACAAAGT"}
    
    # Detect att1 and att2 sites in the entry vector
    att1_entry_start = entry_vector.seq.find(gateway_dict["att1_shared"])
    att1_entry_end = att1_entry_start + len(gateway_dict["att1_shared"]) + 3
    if(att1_entry_start == -1):
        att1_entry_end = 0
        raise ValueError("No att1 site found in the entry vector.")    

    att2_entry_start = entry_vector.seq.find(gateway_dict["att2_shared"]) -3
    if(att2_entry_start == -4):
        att2_entry_start = 0
        raise ValueError("No att2 site found in the entry vector.")
        
    sequence_between_entry = entry_vector.seq[att1_entry_end:att2_entry_start]
    if sequence_between_entry is None:
        raise ValueError("No sequence between att1 and att2 found in the entry vector.")
    

    # Detect att1 and att2 sites in the destination vector
    att1_destination_start = destination_vector.seq.find(gateway_dict["att1_shared"])
    att1_destination_end = att1_destination_start + len(gateway_dict["att1_shared"]) + 3
    if(att1_destination_start == -1):
        att1_destination_end = 0
        raise ValueError("No att1 site found in the destination vector.")
    
    att2_destination_start = destination_vector.seq.find(gateway_dict["att2_shared"]) -3
    if(att2_destination_start == -4):
        att2_destination_start = 0
        raise ValueError("No att2 site found in the destination vector.")
    
    sequence_between_destination = destination_vector.seq[att1_destination_end:att2_destination_start]
    if sequence_between_destination is None:
        raise ValueError("No sequence between att1 and att2 found in the destination vector.")
    
    # Replace sequence_between_destination with sequence_between_entry
    new_sequence = destination_vector.seq[:att1_destination_end] + sequence_between_entry + destination_vector.seq[att2_destination_start:]

    # Create a new SeqRecord for the expression vector
    expression_vector = SeqRecord(
        new_sequence,
        id=new_name,
        name=new_name,
        description=f"Expression vector containing '{entry_name}' shuttled into '{destination_name}'",
        annotations={"molecule_type": "DNA", "topology": "circular"}
    )
    
    # Add features to expression vector
    # Read the features from CSV file
    print(f"Adding features from '{features_file}'")
    with open(features_file, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            Name, Feature, Type, *_ = row
            Name = str(Name)
            Feature = str(Feature).upper()
            Type = str(Type) if Type else "misc_feature"

            # Modify the feature sequence to include alternative bases
            ambiguous_bases = {
                "M": "[AC]", "S": "[GC]", "W": "[AT]", "Y": "[CT]",
                "R": "[AG]", "K": "[GT]", "V": "[ACG]", "H": "[ACT]",
                "D": "[AGT]", "B": "[CGT]", "N": "[ACGT]"
            }
            modified_feature = Feature
            for base, replacement in ambiguous_bases.items():
                modified_feature = modified_feature.replace(base, replacement)

            # Generate the reverse complement of the feature
            reverse_complement_feature = str(Seq(Feature).reverse_complement())
            modified_reverse_complement_feature = reverse_complement_feature
            for base, replacement in ambiguous_bases.items():
                modified_reverse_complement_feature = modified_reverse_complement_feature.replace(base, replacement)

            # Perform the search
            try:
                match = re.search(modified_feature, str(expression_vector.seq))
                if not match:
                    match = re.search(modified_reverse_complement_feature, str(expression_vector.seq))
            except re.error as e:
                print(f"Regex error for feature '{Name}': {e}")
                continue

            if match:
                start = match.start()
                end = match.end()
            else:
                start = -1
                end = -1

            if start != -1:
                if not Type:  # Ensure Type is not empty
                    raise ValueError(f"Feature '{Name}' does not have a valid type.")
                seq_feature = SeqFeature(
                    FeatureLocation(start, end),
                    type=Type,
                    qualifiers={"label": Name}
                )
                expression_vector.features.append(seq_feature)
    
    # Append the entry_name as a feature with type = CDS
    start = expression_vector.seq.find(sequence_between_entry)+2
    end = start + len(sequence_between_entry)-6
    if start != -1:
        entry_feature = SeqFeature(
            FeatureLocation(start, end),
            type="CDS",
            qualifiers={"label": entry_name}
        )
        expression_vector.features.append(entry_feature)
    
    # Save the expression vector map to a file
    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    output = os.path.join(output_folder, new_name)
    SeqIO.write(expression_vector, output, "genbank")
    print(f"Expression vector map saved to {output}")

def process_folder(entry_vector_folder, dest_vector_folder, features_file, output_folder):
    """Process all files in the given folders"""
    for entry_vector_file in os.listdir(entry_vector_folder):
        entry_vector_path = os.path.join(entry_vector_folder, entry_vector_file)
        for destination_vector_file in os.listdir(dest_vector_folder):
            destination_vector_path = os.path.join(dest_vector_folder, destination_vector_file)
            try:
                create_expression_vector(entry_vector_path, destination_vector_path, features_file, output_folder)
            except ValueError as e:
                print(f"Error processing {entry_vector_file} and {destination_vector_file}: {e}. Skipping to the next file.")
                continue

def main():
    parser = argparse.ArgumentParser(description="Create expression vector by performing a virtual LR reaction between an entry and destination vector")
    parser.add_argument("-e", "--entry_vector_folder", type=str, help="The folder containing all entry vector files", default="import/entry_vectors/")
    parser.add_argument("-d", "--dest_vector_folder", type=str, help="The folder containing all destination vector files", default="import/destination_vectors/")
    parser.add_argument("-f", "--features_file", type=str, help="The CSV file containing the features to be added to the expression vector", default="import/features/all_features.csv")
    parser.add_argument("-o", "--output_folder", type=str, help="The output folder", default="output/")
    args = parser.parse_args()
    process_folder(args.entry_vector_folder, args.dest_vector_folder, args.output_folder, args.features_file)

if __name__ == "__main__":
    main()