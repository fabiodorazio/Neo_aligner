import os
import argparse
import time
import sys
import utils
import import_functions as impF
import alignment_functions as algF
import summary_functions as sumF
import graphics_functions as grfF

def get_args():
    parser = argparse.ArgumentParser(prog="neo_aligner")
    parser.add_argument("input_amp", help = "Input amplicon sequences in fasta format")
    parser.add_argument("reference", help ="Input reference genome in fasta format")
    parser.add_argument("coordinates", help = "A bed file containing the coordinates of each amplicon")
    parser.add_argument("-n", dest="name", default=None, help="Output file name without the format: can be the name of the analysis, individual, date, etc...")
    parser.add_argument("-qc", dest="qcmode", action="store_true", help="If provided, runs the program in qc-mode only and does not perform the alignment")
    parser.add_argument("-o", dest="output_dir", default="../Outputs", help="Output directory path")
    parser.add_argument("--combine-amplicons", dest="comb", action="store_true", help="If provided, combines amplicons based on coordinates and basenames")
    parser.add_argument("--mutational-signatures", dest="mut", action="store_true", help="If provided, estimates the type of mutational signature")
    parser.add_argument("--gene-overlaps", dest="gene", action="store_true", help="Provide gene coordinates to identify the ensembl ID of the gene where the mutation is found")

    args = parser.parse_args()

    input = args.input_amp
    ref = args.reference
    coord = args.coordinates
    output = args.output_dir
    name = args.name
    qc = args.qcmode
    comb = args.comb
    mut = args.mut
    gene = args.gene

    return (input, ref, coord, output, name, qc, comb, mut, gene)


if __name__ == "__main__":
    INPUT, REF, COORD, OUTPUT, NAME, QC, COMB, MUT, GENE = get_args()
    print("Starting the aligner...")
    time.sleep(2)

    if NAME is None:
        NAME = utils.get_basename(INPUT)

    # checks if output exists or create one
    output_dir = utils.check_output_dir(OUTPUT)
    # creates output path for the aligner
    output_analysis = output_dir + NAME + ".csv"
    # sets location for plot output, checks if exists or creates one
    output_plot_dir = utils.check_output_dir(output_dir + '/Plot_outputs/')
    output_plot_dir_basename = output_plot_dir + NAME

    ########## IMPORT FILES ##########
    reference = impF.import_ref(REF)
    coordinates = impF.import_coordinates(COORD)
    genes = impF.import_genes(INPUT)

    ########## RUN ALIGNER ##########
    # run aligner in QC only mode
    if QC:
        all_files = os.listdir(INPUT)
        list_of_amplicons = [amp for amp in all_files if 'amplicon' in amp and amp.endswith("fa.gz")]
        for amp in list_of_amplicons:
            print(f"Reading {amp}...")
            amplicons_sequences = impF.import_fasta_qc(INPUT + amp)
        sys.exit(1)

    # >> run aligner 
    aligned_sequences, matrix_heatmap = algF.gather_mutations(reference, coordinates, INPUT)
    
    ########## PERFORM CALCULATIONS ##########
    # mutational signature
    if MUT:
        mut_table = sumF.generate_signature(aligned_sequences)
    else:
        mut_table = aligned_sequences

    final_mut_table = sumF.mutation_frequency(mut_table, output_analysis, COMB)

    # gene overlaps
    if GENE:
        gene = impF.import_genes(INPUT)
        # overlap with provided genes
        output_gene_analysis = output_dir + NAME + "_Genes.csv"
        gene2 = sumF.gene_overlap(gene, final_mut_table, output_gene_analysis)


    ########## GENERATE PLOTS ##########
    grfF.mutation_frequency_substitution(final_mut_table, output_plot_dir_basename)
    grfF.mutation_freq_heatmap(matrix_heatmap, output_plot_dir_basename)

    if MUT:
        grfF.mutation_frequency_signature(final_mut_table, output_plot_dir_basename)
