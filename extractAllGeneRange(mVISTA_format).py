# this program can be execute on python 3.11.9 & biopython API 1.78
# using conda to manage env and activate the biopython env to execute
# plz switch to the folder storing gb file
# enter [python extractAllExon -h] in terminal will show the usage
# function: extract all exon range and strand to fit mVISTA annotation file format
# made by Pin Li 2024/7/19

import os
import glob
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def extract_exon_sequences(filepath):
    gene_rangeAndName_list = []
    times = 0
    species_name = None
    try:
        with open(filepath, 'r') as file:
            for record in SeqIO.parse(file, 'genbank'):
                print("Now be reading: " + record.description)
                if not species_name:
                    species_name = record.annotations.get('organism', 'unknown_species').replace(' ', '_')
                for feature in record.features:
                    #print(feature)
                    if feature.type == "CDS" or feature.type == "tRNA" or feature.type == "rRNA": # 選擇要把哪些類型的 data 抓入
                        times += 1
                        try:
                            gene_name = feature.qualifiers.get("gene")[0]
                            gene_range = feature.location
                            nameRangeTypeTuple = (gene_name, gene_range, feature.type)
                            gene_rangeAndName_list.append(nameRangeTypeTuple)
                            print(gene_name)
                            print(gene_range)
                            print("")
                        except Exception as e: 
                            print(f"Error processing CDS feature in {filepath}: {e}")
                print(f"total parsed gene: {times}\n")
    except IOError as e: # 讀不到檔案
        print(f"Error reading file {filepath}: {e}")
    return gene_rangeAndName_list

def write_to_mVISTA(gene_rangeAndName_list, output_filename):
    symbol = None
    with open(output_filename, 'w') as output:
        for gene_name, gene_range, gene_type in gene_rangeAndName_list:
            match gene_range.strand:
                case 1: # + 正股表示為 >
                    symbol = ">"
                case -1: # - 反股表示為 <
                    symbol = "<"
            if gene_type == "CDS":
                suffix = "exon"
            else: 
                suffix = "utr"
            output.write(f"{symbol} {gene_range.start + 1} {gene_range.end} {gene_name}\n") # gene 層
            for exon in gene_range.parts:
                output.write(f"{exon.start + 1} {exon.end} {suffix}\n") # gene 內的各個 exon 
            output.write("\n")

def set_output_filename(filepath, suffix='.fasta'):
    return os.path.splitext(filepath)[0] + suffix

def output_mVISTA(file_patterns):
        for file_pattern in file_patterns:
            if '*' in file_pattern:  # 處理 * 情況
                for filepath in glob.glob(file_pattern):
                    print(filepath)
                    gene_rangeAndName_list= extract_exon_sequences(filepath)
                    output_filename = set_output_filename(filepath,'.txt')
                    write_to_mVISTA(gene_rangeAndName_list, output_filename)
            else:  # 處理單一檔案或多個檔案情況
                print(file_pattern)
                gene_rangeAndName_list = extract_exon_sequences(file_pattern)
                output_filename = set_output_filename(file_pattern,'.txt')
                write_to_mVISTA(gene_rangeAndName_list, output_filename)

def main():
    parser = argparse.ArgumentParser(description='Extract CDS sequences from GenBank files and output to individual FASTA files.')
    parser.add_argument('files', nargs='+', 
                        help='One or more .gb files or patterns to process (e.g., *.gb or file1.gb file2.gb)')
    args = parser.parse_args()
    output_mVISTA(args.files)
if __name__ == '__main__':
    main()
