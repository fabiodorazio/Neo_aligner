import os
import time
import pandas as pd
import import_functions as impF
import alignment_functions as algF


def generate_signature(mut_table):
    '''
    Mutational substitutions are generated from ref and alt columns and combined in thr 6 canonical classes where the pyrimidine is mutated
    '''
    print("Converting substitutions in mutational signatures...")
    time.sleep(2)
    # nature of substitution
    mut_table["substitution"] = mut_table["reference_allele"] + '>' + mut_table["alternative_allele"]
    
    # create a dictionary to link the substitutions
    repl_dic = {'G>T': 'C>A', 'G>C':'C>G', 'G>A':'C>T', 'A>T':'T>A', 'A>G':'T>C', 'A>C':'T>G'}
    # replace mutations with mut signatures
    mut_table["signature"] = mut_table["substitution"].replace(repl_dic)
     
    return(mut_table)


def mutation_frequency(mut_table, output_path, combine_amplicons: bool = True):
    '''
    Retrieve total number of sequences per amplicon
    Retrieve total number of mutations
    Combine amp1 and amp2
    Calculate frequency: the number of mutations should be divided by the number of amplicons. Overlapping sequences will be counted twice
    Number of amplicons supporting the mutation
    '''
    # Combine amplicons
    if combine_amplicons:
        # create a new df to combine amplicons with overlapping sequences
        mut_table_name = mut_table
        # check if there are overlapping mutations on different amplicons
        overlapping_amp = mut_table_name.groupby("coordinates").filter(lambda x: x["amplicon"].nunique() > 1)
        # extract amplicon suffix for new name
        suffixes = []
        for name in overlapping_amp["amplicon"]:
            suffix = name.split("_")[-1]
            suffixes.append(suffix)
        # concatenate all elements in the list
        amplicon_combined = '_'.join(suffixes)
        # replace amplicon names
        mut_table_name["orig_amplicon"] = mut_table_name["amplicon"]
        mut_table_name.loc[overlapping_amp.index, "amplicon"] = "amplicon_" + amplicon_combined

        # combine amplicons with the same name and sum the total number of sequences
        mut_table_name["combined_total"] = mut_table_name.groupby("amplicon")["sequence_total"].transform("sum")
        # sum the number of mutations when grouping by coordinate and type of mutation
        mut_table_name["combined_count"] = mut_table_name.groupby(["coordinates", "substitution"])["count"].transform("sum")

        # calculate mutation frequency
        mut_table_name["frequency"] = mut_table_name["combined_count"]/mut_table_name["combined_total"]
    

    # Do not combine amplicons, instead calculate frequency based on iteration: each amplicon input is independent
    else:
        # create a new df to combine amplicons with overlapping sequences
        mut_table_name = mut_table
        mut_table_name["frequency"] = mut_table_name["count"]/mut_table_name["sequence_total"]

    # save
    mut_table_name.to_csv(output_path, index=False)

    return(mut_table_name)
    

def mutational_signature(mut_table):
    '''
    Calculates the percentage of mutational signature across all amplicons: it assumes all amplicons come from the same individual
    '''
    # calculate total number of mutations detected in the assay
    tot_mutations = mut_table["count"].sum()

    # group mut counts based on the signature type
    mut_table_sign = mut_table.groupby("signature")["count"].sum()

    # calculate percentage of all mutations identified
    mut_table_sign_perc = (mut_table_sign/tot_mutations)*100
    
    return(mut_table_sign_perc)


def gene_overlap(gene, mut_table, output_path):
    '''
    Overlap mutations to gene by coordinates
    Identify the exon where the mutation falls
    Some mutations may fall outside the exon coordinates. 
    To assign them to the gene, the min start and max end of the genes are obtained
    (The above would need strand info for other inputs)
    '''
    # create list to store ensembl ids
    Ensembl_list = []
    gene_list = []

    # iterate through rows in the mutation table to extract the mutation coordinates
    for index, row in mut_table.iterrows():
        # get coordinates of the mutations
        mutation = row["coordinates"]
        # find within which range the mutation falls and extract the ensembl id and the gene
        gene_overlap = gene[gene.apply(lambda x: (x["start"] <= mutation <= x["end"]), axis=1)]
        # add gene and ensembl to list if they are found, NA otherwise
        if not gene_overlap.empty:
            Ensembl_list.append(gene_overlap.iloc[0]["ensembl_ID"])
            #gene_list.append(gene_overlap.iloc[0]["gene"])
        else:
            Ensembl_list.append("NA")
            #gene_list.append("NA")

        # get min start and max end for each gene (strand not included)  
        gene_coordinates = gene.groupby("gene").agg({"start": "min", "end": "max"}).reset_index()
        # 
        assign_gene = gene_coordinates[gene_coordinates.apply(lambda x: (x["start"] <= mutation <= x["end"]), axis=1)]
        if not assign_gene.empty:
            gene_list.append(assign_gene.iloc[0]["gene"])
        else:
            gene_list.append("NA")

    # add elements to dataframe
    mut_table["ensembl_ID"] = Ensembl_list
    mut_table["gene"] = gene_list
    
    # save
    mut_table.to_csv(output_path, index=False)

    return(mut_table)




