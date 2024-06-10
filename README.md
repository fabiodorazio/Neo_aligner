
**Neo Aligner**


## Introduction

This is a simple aligner for short sequences using the Align module from Biopython.
The aim of this program is to identify single nucleotide mutations in provided amplicon sequences by mapping them to a subset of the reference genome.
The program generates summary tables and plots in the Output directory

## Quick Start

1. Install [`Conda`](https://conda.io/miniconda.html)

2. Install the required packages from the provided environment.yml file with the following command:

```bash
conda env create -f environment.yml
```

3. cd into the bin directory

4. Start the aligner with the following command:

```bash
python3 main.py ../assets/ ../assets/genome.fa.gz ../assets/amplicon_coordinates.bed -n Mut_Signatures_Combined_Amplicons --mutational-signatures --combine-amplicons --gene-overlaps
```

## Pipeline Summary

The following 3 arguments are mandatory for neo_aligner:

1. "input_amp"    >  "Input amplicon sequences in fasta format"
2. "reference"    >  "Input reference genome in fasta format"
3. "coordinates"  >  "A bed file containing the coordinates of each amplicon"

The following 6 arguments are optional:

1. "-n","name"                  >  "Output file name without the format: can be the name of the analysis, individual, date, etc..."
2. "-qc", "qcmode"              >  "If provided, runs the program in qc-mode only and does not perform the alignment"
3. "-o", est="output_dir"       >  "Output directory path"
4. "--combine-amplicons"        >  "If provided, combines amplicons based on coordinates and basenames"
5. "--mutational-signatures"    >  "If provided, estimates the type of mutational signature"
6. "--gene-overlaps"            >  "Provide gene coordinates to identify the ensembl ID of the gene where the mutation is found"

## Description

Neo_aligner starts by importing the required files provided by the user. It performs some qc checks on the input files:

1. For each amplicon, returns the number of sequences and the length of each
2. Reads the reference genome and looks for canonical and non-canonical bases. The program is terminated if non-expected bases are found

To perform the alignment, Neo_aligner uses the Align module from Biopython. Due to memory limitations of this module, the program subsets the reference genome based on the amplicon coordinates provided, so to align the amplicon only to its expected target.
Note: this will not work if there are other types of mutations either than single nucleotide substitutions

For each alignment performed, the coordinate of the mutated based found is registered, alongside with a matrix that stores the position and the sequence number of each mutation.
The alignment is run for each amplicon provided in the input directory and the reference and alternative alleles are retrieved

Neo_aligner aims to perform a simple mutational signature convertion. This is particularly important for cancer samples where the origin is uncertain.

If target gene/exons coordinates are provided, each mutation can be mapped back and labelled with the gene and exons where it belongs.

Finally, the mutation frequency is calculated for substitutions or signatures and graphics are generated