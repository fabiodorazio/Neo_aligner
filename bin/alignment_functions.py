import pandas as pd
import numpy as np
from Bio import Align
import import_functions as impF
import os


def assign_coordinates(amplicon, coordinates):
    '''
    Retrieves the coordinates of the amplicon
    Adjusts the coordinates provided so to match the lengths of the amplicons
    Defines the new start position based on the lengths of the amplicons
    '''
    # extract amplicon name and sequences from dict
    amp_name = list(amplicon.keys())[0]
    amp_sequence = list(amplicon.values())[0]
    # get file name and map it the amplicon to the coordinates
    amp_row = coordinates.loc[coordinates["amplicon"] == amp_name]
    # get start and end coordinates
    start = amp_row["start"]
    end = amp_row["end"]

    # amplicon length
    amp_length = int(end - start)

    # obtain max length from the list of sequences
    l = max(len(seq) for seq in amp_sequence)

    # adjust start coordinate
    # if the coordinates are shorter than the max amplicon, shift the start downstream
    if amp_length < l:
        difference = int(amp_length - l)
        new_start = int(start + difference)
    # if the coordinates are longer than the max amplicon, shift the start upstream
    elif amp_length > l:
        difference = int(amp_length - l)
        new_start = int(start + difference)
    else:
        new_start = start
        
    # return: start, end, max length, 
    output_dict = {
        "amplicon": amp_name,
        "start": new_start,
        "end": end,
        "max_length": l
    }

    output_df = pd.DataFrame(output_dict)

    return(output_df)


def aligner(seqA, seqB, coord_start):
    # instance of aligner class
    aligner = Align.PairwiseAligner()
    # perform alignment
    alignments = aligner.align(seqA, seqB)
    # retrieve alignments
    best_alignment = alignments[0]
    seqA_aligned, seqB_aligned = best_alignment.target, best_alignment.query
    
    # define start coordinate
    coordinate = int(coord_start["start"])

    # count mutations
    mutations = 0
    coordinate_list = []
    # iterate over each sequence after combining them in a single iterable
    for ref, alt in zip(seqA_aligned, seqB_aligned):
        coordinate += 1
        if ref != alt and ref != "-" and alt != "-":
            mutations += 1
            # append coordinate
            coordinate_list.append(coordinate)
            # append ref base
            coordinate_list.append(ref)
            # append alt base
            coordinate_list.append(alt)
     
    return(coordinate_list)


def run_alignment(amplicons_info, genome, coords):
    '''
    Retrieves coordinates and amplicon sequences
    Runs the aligner in a pairwise manner between each amplicon and the subset reference genome for memory efficiency
    Gathers all the mutations in a list of lists 
    Creates a matrix for heatmaps
    '''
    # retrieve start and end coordinates
    start = int(coords["start"])
    end = int(coords["end"])

    # retrieve amplicon sequences
    amplicons_sequences = list(amplicons_info.values())[0]
    
    # run alignmemnt and store coordinates and bases of mismatches in lists
    coordinates_list = []
    ref_list = []
    alt_list = []

    # create a matrix to store binarised mutation positions
    mut_matrix = pd.DataFrame(0, index = range(1,len(amplicons_sequences)+1), columns = range(start,end))

    # add counter for tracking
    n = 1
    for sequence in amplicons_sequences:
        print(f"Aligning sequence {n} to reference {start}:{end}...")
        n +=1

        # run the aligner on each sequence
        l = aligner(sequence, genome[start:end], coords)
        print(f"Found {l} mismatch/es")
        #time.sleep(0.01)

        if l:
            coordinates_list.append(l[0])
            ref_list.append(l[1])
            alt_list.append(l[2])
            # replace 0s with 1s when a mutation occurs
            mut_matrix.loc[n,l[0]] = 1

    # combine output in a df
    output_dict = {
        "coordinates":coordinates_list,
        "reference_allele":ref_list,
        "alternative_allele":alt_list 
    }

    output_df = pd.DataFrame(output_dict)
    
    return(output_df, mut_matrix)


def gather_mutations(genome, coords, input_path):
    '''
    Runs the aligner for each provided amplicon file
    Combines the mutations form different sequences for the same amplicon
    Counts the total number of mutations found and total number of sequences
    Each iteration is labelled so to allow to analyse each amplicon separately during downstream analyses
    '''
    # list amplicon files
    all_files = os.listdir(input_path)
    list_of_amplicons = [amp for amp in all_files if 'amplicon' in amp and amp.endswith("fa.gz")]

    # label each iteration
    lab_iteration = 1

    # create a final dataframe for concatenation
    final_df = pd.DataFrame()
    # create a final dictionary for plotting heatmaps
    matrix_dict = {}

    # iterate over list of files
    for amp in list_of_amplicons:
        lab_iteration += 1
        # read file
        amplicons_sequences = impF.import_fasta_qc(input_path + amp)

        # retrieve total number of sequences
        tot_seq = len(list(amplicons_sequences.values())[0])
        # retrieve amplicon base name
        amplicon_name = amp.split(".")[0]
        # assign coordinates to amplicon
        assigned_coordinates = assign_coordinates(amplicons_sequences, coords)
        # run aligner
        coordinates_list, matrix = run_alignment(amplicons_sequences, genome, assigned_coordinates)

        # add matrix to dictionary for later plotting
        matrix_dict[amp] = matrix
        # generate outputs: amp name, coordinate, mut count, total sequences, n of iterations and nature of substitution

        # nature of substitution
        coordinates_list["substitution"] = coordinates_list["reference_allele"] + '>' + coordinates_list["alternative_allele"]
        df = coordinates_list.groupby(coordinates_list.columns.tolist()).size().reset_index(name = "count")
        df["iteration"] = lab_iteration
        df["amplicon"] = amplicon_name
        df["sequence_total"] = tot_seq

        final_df = pd.concat([final_df, df], ignore_index=True)

    return(final_df, matrix_dict)




