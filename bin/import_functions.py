import gzip
from Bio import SeqIO
import time
import pandas as pd
import os

def import_fasta_qc(gz_file):
    '''
    Imports sequences in fa.gz format
    Prints basic info about the file
    '''
    with gzip.open(gz_file, "rt") as sequences:
        # check number and lengths of amplicons in the file
        amplicons = 0
        lengths = []
        counts = {}
        
        all_sequences = []

        for record in SeqIO.parse(sequences, "fasta"):
            amplicons += 1
            lengths.append(len(record.seq))
            all_sequences.append(str(record.seq))

        #print("Reading the input file")
        print(f"There is a total of {amplicons} sequences in the fasta file")

        # fill the counts dictionary with each unique value in lengths and the sum
        for l in lengths:
            if l in counts:
                counts[l] += 1
            else:
                counts[l] = 1
        # prints out the items in the dictionary
        for number, frequency in counts.items():
            print(f"There are {frequency} sequences of length {number}")
        time.sleep(3)

        # retrieve file name and save in a dictionary
        file_name = gz_file.split("/")[-1].split(".")[0]
        dict_sequences = {file_name: all_sequences}

        return(dict_sequences)
    

def import_ref(reference):
    '''
    Imports reference genome 
    Performs some basic checks: nucleotide content and unique sequence
    '''
    with gzip.open(reference, "rt") as ref_genome:
        print("Reading the genome reference file...")
        time.sleep(2)
        for record in SeqIO.parse(ref_genome, "fasta"):
            # check nucleotide content
            sequence = str(record.seq).upper()
            non_canonical = set(sequence) - set("ATGCN") # finds non canonical bases in the reference genome
            prohibited = set(sequence) - set("ATGCNRYSWKMBDHV") # finds bases that are not expected in reference
            if non_canonical:
                print(f'Non canonical bases found in reference genome: {non_canonical}')
                time.sleep(3)
            if prohibited:
                print(f'Prohibited bases found in reference genome: {prohibited}')
                time.exit(1)

            # check unique sequence

    return(sequence)


def import_coordinates(coordinates):
    '''
    Imports coordinates in bed format
    '''
    # read bed
    bed_coordinates = pd.read_csv(coordinates, sep = '\t', header = None)
    # assign column names
    bed_coordinates.columns = ["chr", "start", "end", "amplicon"]
    print("importing amplicon coordinates...")
    time.sleep(3)
    print("...done")
    return(bed_coordinates)


def import_genes(gene_path):
    '''
    Load gene files from input path
    Add a gene label and combine all genes in a single df
    '''
    # load genes coordinates from bed files
    all_files = os.listdir(gene_path)
    gene_files = [gene for gene in all_files if 'gene' in gene and gene.endswith("bed")]

    # create a dataframe to combine all genes
    gene_df = pd.DataFrame()

    for gene in gene_files:
        # get gene basename
        gene_basename = gene.split(".")[0]
        gene_coordinates = pd.read_csv(gene_path + gene, sep = '\t', header = None)
        # assign column names
        gene_coordinates.columns = ["ensembl_ID", "chr", "start", "end", "width"]
        # add label
        gene_coordinates["gene"] = gene_basename

        # append to common dataframe
        gene_df = pd.concat([gene_df, gene_coordinates], ignore_index=True)

    return(gene_df)




