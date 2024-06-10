import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme(color_codes=True)


def mutation_freq_heatmap(matrix_dict, output_plot_path):
    '''
    Iterate over items in the input dictionary 
    Plot heatmaps of matrices 
    Save
    '''
    # plot
    for name, matrix in matrix_dict.items():
        sns.heatmap(matrix)
        plt.savefig(f"{output_plot_path}_{name}_Mut_Heatmap.png", bbox_inches = "tight")
        plt.close()


def mutation_frequency_substitution(mut_table, output_plot_path):
    '''
    Generates a barplot showing the mutation frequency of each substitution found
    '''
    # mutation frequency by substitution
    sub = mut_table.groupby("substitution")["count"].sum()
    total = mut_table["count"].sum()

    f1 = sub/total
    f1 = f1.reset_index()
    f1.columns = ["substitution", "frequency"]
    # plot
    sns.barplot(f1, x='substitution', y='frequency')
    plt.savefig(f'{output_plot_path}_Mut_Sub_Freq.png')
    plt.close()


def mutation_frequency_signature(mut_table, output_plot_path):
    '''
    Generates a barplot showing the mutation frequency of each signature found
    '''
    # mutation frequency by substitution
    sign = mut_table.groupby("signature")["count"].sum()
    tot = mut_table["count"].sum()

    f1 = sign/tot
    f1 = f1.reset_index()
    f1.columns = ["signature", "frequency"]
    # plot
    sns.barplot(f1, x="signature", y="frequency")
    plt.savefig(f"{output_plot_path}_Mut_Sign_Freq.png")
    plt.close()

