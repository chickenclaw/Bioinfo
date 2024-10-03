# this program can be execute on python 3.11.9 & biopython API 1.78
# using conda to manage env and activate the biopython env to execute
# enter [python extractProteinCodingGene -h] in terminal will show the usage
# function: extract all nucleotides sequence of CDS from gb file, and produce fasta file
# hope the program can help you process data more quickly than manual method. Have fun~
# made by Pin Li 2024/6/17

import os
import glob
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def extract_cds_sequences(filepath):
    cds_sequences = []
    species_name = None
    try:
        with open(filepath, 'r') as file:
            for record in SeqIO.parse(file, 'genbank'):
                print("Now be reading: " + record.description)
                if not species_name:
                    species_name = record.annotations.get('organism', 'unknown_species').replace(' ', '_')
                for feature in record.features:
                    if feature.type == "CDS":
                        try:
                            cds_seq = feature.extract(record.seq)
                            gene = feature.qualifiers.get('gene', [''])[0]
                            print("")
                            print(feature.qualifiers.get("gene")[0])
                            print(cds_seq)
                            cds_record = SeqRecord(cds_seq, id=gene, description='')
                            cds_sequences.append(cds_record)
                        except Exception as e:
                            print(f"Error processing CDS feature in {filepath}: {e}")
    except IOError as e:
        print(f"Error reading file {filepath}: {e}")
    return cds_sequences, species_name

def write_to_fasta(sequences, output_filename, id_name):
    with open(output_filename, 'w') as output_handle:
        SeqIO.write(sequences, output_handle, 'fasta')

def set_output_filename(filepath, suffix='.fasta'):
    return os.path.splitext(filepath)[0] + suffix

def mode1_output(file_patterns):
    for file_pattern in file_patterns:
        if '*' in file_pattern:  # 處理 * 情況
            for filepath in glob.glob(file_pattern):
                cds_sequences, species_name = extract_cds_sequences(filepath)
                concatenated_seq = SeqRecord(Seq(''.join([str(seq.seq) for seq in cds_sequences])), id=species_name, description='')
                output_filename = set_output_filename(filepath)
                write_to_fasta([concatenated_seq], output_filename, species_name)
        else:  # 處理單一檔案或多個檔案情況
            cds_sequences, species_name = extract_cds_sequences(file_pattern)
            concatenated_seq = SeqRecord(Seq(''.join([str(seq.seq) for seq in cds_sequences])), id=species_name, description='')
            output_filename = set_output_filename(file_pattern)
            write_to_fasta([concatenated_seq], output_filename, species_name)

def mode2_output(file_patterns):
    for file_pattern in file_patterns:
        if '*' in file_pattern:  # 處理 * 情況
            for filepath in glob.glob(file_pattern):
                process_single_file_mode2(filepath)
        else:  # 處理單一檔案或多個檔案情況
            process_single_file_mode2(file_pattern)

def process_single_file_mode2(filepath):
    cds_sequences, species_name = extract_cds_sequences(filepath)
    folder_name = set_output_filename(filepath, '')
    os.makedirs(folder_name, exist_ok=True)
    
    gene_count = {}
    
    for seq_record in cds_sequences:
        gene_name = seq_record.id
        if gene_name in gene_count:
            gene_count[gene_name] += 1
            gene_name = f"{gene_name}_{gene_count[gene_name]}"
        else:
            gene_count[gene_name] = 1
            gene_name = f"{gene_name}_1"
        
        output_filename = os.path.join(folder_name, gene_name + '.fasta')
        write_to_fasta([seq_record], output_filename, os.path.basename(filepath))

def main():
    parser = argparse.ArgumentParser(description='Extract CDS sequences from GenBank files and output to individual FASTA files.',
                                     epilog= '-mode 1 or -mode 2 is necessary args')
    parser.add_argument('files', nargs='+', 
                        help='One or more .gb files or patterns to process (e.g., *.gb or file1.gb file2.gb)')
    parser.add_argument('-mode', choices=['1', '2'], required=True, 
                        help='Output mode: 1 (concatenate sequences) or 2 (separate gene files)')
    args = parser.parse_args()

    if args.mode == '1':
        mode1_output(args.files)
    elif args.mode == '2':
        mode2_output(args.files)

if __name__ == '__main__':
    main()
