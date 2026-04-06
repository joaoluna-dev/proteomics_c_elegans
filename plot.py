# Importando módulos
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import seaborn as sns
import os
import re

# Script de plotagem

# Dicionário de mapeamento de sufixos de arquivo para nomes padronizados
TITLE_MAP = {
    "GO_BP": "GO: Biological Process",
    "GO_CC": "GO: Cellular Component",
    "GO_MF": "GO: Molecular Function",
    "Kegg_Pathways": "KEGG Pathways",
    "Reactome_Pathways": "Reactome Pathways",
}


def get_display_title(filename_base):
    """Obtém o título padronizado a partir do nome do arquivo."""
    for suffix, display_name in TITLE_MAP.items():
        if filename_base.endswith(suffix):
            return display_name
    return filename_base

def plot_subplot_enrichment(ax, enriched_data, plot, titlename):
    """Plota dados de enriquecimento em um eixo de subplot."""
    # Filtrando dados
    try:
        enriched_data = enriched_data.copy()
        enriched_data['P-Value'] = enriched_data['P-Value'].astype(float)
        enriched_data['FDR'] = enriched_data['FDR'].astype(float)
        enriched_data = enriched_data[enriched_data['P-Value'] <= 0.05]
        enriched_data = enriched_data[enriched_data['FDR'] <= 0.05]
    except Exception as e:
        print(f"Erro ao filtrar os dados por P-value e FDR: {e}")

    # Obtendo os 10 termos mais representativos e seus valores
    sorted_object = enriched_data.sort_values(by='Count', ascending=False)
    top_10_sorted_object = sorted_object.head(10)
    terms = top_10_sorted_object['Term']
    pvalue = top_10_sorted_object['P-Value']
    proteins_count = top_10_sorted_object['Count']

    # Configurando o subplot como uma caixinha independente
    ax.set_frame_on(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
        spine.set_edgecolor('black')

    # Adicionando fundo branco ao subplot para separação visual
    ax.set_facecolor('white')

    if len(terms) == 0:
        ax.set_xlabel("Count")
        ax.set_ylabel("Terms")
        display_title = get_display_title(titlename)
        ax.text(1.12, 0.5, display_title, transform=ax.transAxes, rotation=90,
                va='center', ha='left', fontsize=12, fontweight='bold',
                bbox=dict(boxstyle='square,pad=0.4', fc='lightgray', ec='none'))
        return

    if plot == "scatterplot":
        enriched_plot = sns.scatterplot(x=proteins_count,
                                        y=terms,
                                        size=proteins_count,
                                        hue=pvalue,
                                        sizes=(100, 1000),
                                        palette='coolwarm_r',
                                        legend=False,
                                        ax=ax)
    elif plot == "barplot":
        enriched_plot = sns.barplot(x=proteins_count,
                                    y=terms,
                                    orient="h",
                                    hue=pvalue,
                                    dodge=False,
                                    palette='coolwarm_r',
                                    ax=ax)
        enriched_plot.get_legend().remove()
    else:
        raise ValueError(f"Tipo de plot desconhecido: {plot}")

    enriched_plot.set_xlabel("Count")
    enriched_plot.set_ylabel("Terms")

    # Criando legenda de p-values individual para este subplot
    normalized_pvalues = Normalize(vmin=pvalue.min(), vmax=pvalue.max())
    sm = ScalarMappable(cmap='coolwarm_r', norm=normalized_pvalues)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.01, fraction=0.04)
    cbar.set_label('P-value', rotation=270, labelpad=10)

    # Adicionando título na lateral direita com fundo cinza (estilo margin_titles do FacetGrid)
    display_title = get_display_title(titlename)
    ax.text(1.12, 0.5, display_title, transform=ax.transAxes, rotation=90,
            va='center', ha='left', fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='square,pad=0.4', fc='lightgray', ec='none'))


try:
    enrichment_data_dir = "data/output/enrichment_data"
    if not os.path.exists(enrichment_data_dir):
        enrichment_data_dir = input("Insira o caminho da pasta com os arquivos CSV de enriquecimento \n"
                                    "Certifique-se de que não há outros arquivos .CSV na mesma pasta: ")
    plot_type = sys.argv[1]
    enrichment_plot_dir = "data/output/enrichment_plots"
    if not os.path.exists(enrichment_plot_dir):
        os.mkdir(enrichment_plot_dir)
except IndexError:
    print("Uso: python plot.py plot_type")
    sys.exit(1)
except Exception as e:
    print(f"Erro durante a obtenção dos argumentos: {e}")
    sys.exit(1)


# Lendo e categorizando os arquivos CSV dinamicamente
go_data = {}
pathway_data = {}

for file in os.listdir(enrichment_data_dir):
    if not file.endswith(".csv"):
        continue
    try:
        df = pd.read_csv(os.path.join(enrichment_data_dir, file))
        file_base = file.split(".")[0]

        # Categorizando arquivos com base nos padrões de nomenclatura
        if re.search(r'_GO_(BP|MF|CC)$', file_base):
            go_data[file_base] = df
        elif re.search(r'_(Kegg|Reactome)_', file_base):
            pathway_data[file_base] = df
    except Exception as e:
        print(f"Erro ao ler o arquivo {file}: {e}")


# Plotando subplots de GO (BP, CC, MF) agrupados
if go_data:
    go_order = ['BP', 'CC', 'MF']
    go_keys = sorted(go_data.keys(), key=lambda x: go_order.index(x.split('_')[-1]) if x.split('_')[-1] in go_order else 999)
    n_go = len(go_keys)

    fig, axes = plt.subplots(n_go, 1, figsize=(17, 7 * n_go))
    if n_go == 1:
        axes = [axes]

    all_go_cbs = []
    for idx, key in enumerate(go_keys):
        ax = axes[idx]
        plot_subplot_enrichment(ax, go_data[key], plot_type, key)

    fig.subplots_adjust(right=0.82, hspace=0.15, left=0.12, top=0.95, bottom=0.05)
    fig.suptitle("GO Enrichment Analysis", x=0.5, y=0.98, fontsize=16, fontweight='bold')
    plt.savefig(os.path.join(enrichment_plot_dir, "GO_Enrichment.tiff"), dpi=300, format="tiff", bbox_inches='tight')
    plt.close()
    print(f"Plot de enriquecimento GO salvo com {n_go} subplots.")


# Plotando subplots de vias (Kegg, Reactome) agrupados
if pathway_data:
    pathway_order = ['Kegg', 'Reactome']
    pathway_keys = sorted(pathway_data.keys(), key=lambda x: next((i for i, p in enumerate(pathway_order) if p in x), 999))
    n_pathways = len(pathway_keys)

    fig, axes = plt.subplots(n_pathways, 1, figsize=(17, 7 * n_pathways))
    if n_pathways == 1:
        axes = [axes]

    for idx, key in enumerate(pathway_keys):
        ax = axes[idx]
        plot_subplot_enrichment(ax, pathway_data[key], plot_type, key)

    fig.subplots_adjust(right=0.82, hspace=0.15, left=0.12, top=0.95, bottom=0.05)
    fig.suptitle("Pathway Enrichment Analysis", x=0.5, y=0.98, fontsize=16, fontweight='bold')
    plt.savefig(os.path.join(enrichment_plot_dir, "Pathway_Enrichment.tiff"), dpi=300, format="tiff", bbox_inches='tight')
    plt.close()
    print(f"Plot de enriquecimento de vias salvo com {n_pathways} subplots.")



