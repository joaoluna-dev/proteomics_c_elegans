#Importing modules
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import seaborn as sns
import os

#Plotting script

def plot_enrichment_data(enriched_data, plot, plotsdir, titlename):
    if plot == "scatterplot":

        # Filtering data
        try:
            enriched_data['P-Value'] = enriched_data['P-Value'].astype(float) # Converting P-values to float
            enriched_data['FDR'] = enriched_data['FDR'].astype(float) # Converting FDR values to float
            enriched_data = enriched_data[enriched_data['P-Value'] <= 0.05] # Filtering P-values <= 0.05
            enriched_data = enriched_data[enriched_data['FDR'] <= 0.05] # Filtering FDR values <= 0.05
        except Exception as e:
            print(f"Erro ao filtrar os dados por P-value e FDR: {e}")

        # Getting top 20 most represented terms and his values
        sorted_object = enriched_data.sort_values(by='Count', ascending=False)
        top_20_sorted_object = sorted_object.head(20)
        terms = top_20_sorted_object['Term']
        pvalue = top_20_sorted_object['P-Value']
        proteins_count = top_20_sorted_object['Count']

        # Creating base plot
        plt.figure(figsize=(17, 10))
        enriched_plot = sns.scatterplot(x=proteins_count,
                                        y=terms,
                                        size=proteins_count,
                                        hue=pvalue,
                                        sizes=(100,1000),
                                        palette='coolwarm_r',
                                        legend=False)
        enriched_plot.set_title(f"Top 20 terms from {titlename}")
        enriched_plot.set_xlabel("Count")
        enriched_plot.set_ylabel("Terms")

        # Creating adjusted p-values label
        normalized_pvalues = Normalize(vmin=pvalue.min(), vmax=pvalue.max())
        sm = ScalarMappable(cmap='coolwarm_r', norm=normalized_pvalues)
        sm.set_array([])
        colorbar = plt.colorbar(sm, ax=enriched_plot, orientation='vertical', pad=0.05)
        colorbar.set_label('P-value', rotation=270, labelpad=15)

        # Resizing elements and saving plot to a file
        plt.tight_layout()
        plt.savefig(os.path.join(plotsdir, f"{titlename}.tiff"), dpi=300, format="tiff", bbox_inches='tight')
        plt.close()
        print(f"Scatterplot {titlename} plotado.")

    elif plot == "horizontal barplot":

        # Filtering data
        try:
            enriched_data['P-Value'] = enriched_data['P-Value'].astype(float)  # Converting P-values to float
            enriched_data['FDR'] = enriched_data['FDR'].astype(float)  # Converting FDR values to float
            enriched_data = enriched_data[enriched_data['P-Value'] <= 0.05]  # Filtering P-values <= 0.05
            enriched_data = enriched_data[enriched_data['FDR'] <= 0.05]  # Filtering FDR values <= 0.05
        except Exception as e:
            print(f"Erro ao filtrar os dados por P-value e FDR: {e}")

        # Getting top 20 most represented terms and his values
        sorted_object = enriched_data.sort_values(by='Count', ascending=False)
        top_20_sorted_object = sorted_object.head(20)
        terms = top_20_sorted_object['Term']
        pvalue = top_20_sorted_object['P-Value']
        proteins_count = top_20_sorted_object['Count']

        # Creating base plot
        plt.figure(figsize=(17, 10))
        enriched_plot = sns.barplot(x=proteins_count,
                                    y=terms,
                                    orient="h",
                                    hue=pvalue,
                                    dodge=False,
                                    palette='coolwarm_r')
        enriched_plot.get_legend().remove()
        enriched_plot.set_title(f"Top 20 terms from {titlename}")
        enriched_plot.set_xlabel("Count")
        enriched_plot.set_ylabel("Terms")

        # creating log adjusted p-values label
        normalized_pvalues = Normalize(vmin=pvalue.min(), vmax=pvalue.max())
        sm = ScalarMappable(cmap='coolwarm_r', norm=normalized_pvalues)
        sm.set_array([])
        colorbar = plt.colorbar(sm, ax=enriched_plot, orientation='vertical', pad=0.05)
        colorbar.set_label('P-value', rotation=270, labelpad=15)

        # Resizing elements and saving plot to a file
        plt.tight_layout()
        plt.savefig(os.path.join(plotsdir, f"{titlename}.tiff"), dpi=300, format="tiff", bbox_inches='tight')
        plt.close()
        print(f"Barplot horizontal {titlename} plotado.")



try:
    enrichment_data_dir = sys.argv[1]
    plot_type = sys.argv[2]
    enrichment_plot_dir = sys.argv[3]
except IndexError:
    print("Uso: python plot.py enrichment_data_dir plot_type enrichment_plot_dir")
    sys.exit(1)
except Exception as e:
    print(f"Erro durante a obtenção dos argumentos: {e}")
    sys.exit(1)


# Calling plot_enrichment_data for each .csv file located in enrichment_data_dir
for file in os.listdir(enrichment_data_dir):
    if ".csv" in file:
        try:
            df = pd.read_csv(os.path.join(enrichment_data_dir, file))
            df_titlename = file.split(".")[0]
            plot_enrichment_data(
                enriched_data=df,
                plot=plot_type,
                plotsdir=enrichment_plot_dir,
                titlename=df_titlename
            )
        except Exception as e:
            print(f"Erro durante a plotagem dos dados: {e}")
            sys.exit(1)



