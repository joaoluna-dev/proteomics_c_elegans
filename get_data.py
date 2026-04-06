import os
import sys
import pandas as pd
import omicscope as omics
import matplotlib.pyplot as plt
import math


def read_proteomics_file(rawfilepath, method, control, foldchange):
    """Lê o arquivo de dados de proteômica e retorna um objeto OmicScope."""
    # Remove espaços e aspas do caminho
    rawfilepath = rawfilepath.strip().strip('"')

    if not os.path.exists(rawfilepath):
        raise FileNotFoundError(f"Arquivo com os dados de proteômica não encontrado: {rawfilepath}. Verifique o caminho inserido e tente novamente.")
    elif os.path.getsize(rawfilepath) == 0:
        raise Exception(f"O arquivo com os dados da proteômica fornecido está vazio: {rawfilepath}. Verifique o conteúdo do mesmo e tente novamente.")

    try:
        data = omics.OmicScope(rawfilepath, Method=method, ControlGroup=control, FoldChange_cutoff=foldchange, PValue_cutoff=0.05)
        return data
    except Exception as e:
        raise Exception(f"Erro ao ler o arquivo de proteômica: {str(e)}")


def plot_data(rawfiledata, analysis_type, titlename, plotsdir, analysis, proteins):

    if analysis_type == "conditions_barplot":
        try:
            print(f"Gerando {analysis}...")
            rawfiledata.bar_protein(save=os.path.join(plotsdir, f"{titlename}_{analysis}"), dpi=300, vector=False)
        except Exception as e:
            print(f"Erro ao gerar o barplot de condições: {str(e)}")
            sys.exit(1)

    elif analysis_type == "conditions_boxplot":
        try:
            print(f"Gerando {analysis}...")
            rawfiledata.boxplot_protein(save=os.path.join(plotsdir, f"{titlename}_{analysis}"), dpi=300, vector=False)
        except Exception as e:
            print(f"Erro ao gerar o boxplot de condições: {str(e)}")
            sys.exit(1)

    elif analysis_type == "dynamic_range":
        try:
            print(f"Gerando {analysis}...")
            rawfiledata.DynamicRange(save=os.path.join(plotsdir, f"{titlename}_{analysis}"), dpi=300, vector=False)
        except Exception as e:
            print(f"Erro ao gerar o dynamic range: {str(e)}")
            sys.exit(1)


    elif analysis_type == "ma_plot":
        try:
            print(f"Gerando {analysis}...")
            rawfiledata.MAplot(*proteins, save=os.path.join(plotsdir, f"{titlename}_{analysis}"), dpi=300, vector=False)
        except Exception as e:
            print(f"Erro ao gerar o MA plot: {str(e)}")
            sys.exit(1)

    elif analysis_type == "correlation_heatmap":
        try:
            print(f"Gerando {analysis}...")
            rawfiledata.correlation(linewidth=0, save=os.path.join(plotsdir, f"{titlename}_{analysis}"), dpi=300, vector=False)
        except Exception as e:
            print(f"Erro ao gerar o correlation heatmap: {str(e)}")
            sys.exit(1)

    elif analysis_type == "expression_heatmap":
        try:
            print(f"Gerando {analysis}...")
            rawfiledata.heatmap(linewidth=0, save=os.path.join(plotsdir, f"{titlename}_{analysis}"), dpi=300, vector=False)
        except Exception as e:
            print(f"Erro ao gerar o heatmap: {str(e)}")
            sys.exit(1)

    elif analysis_type == "id_barplot":
        try:
            print(f"Gerando {analysis}...")
            rawfiledata.bar_ident(save=os.path.join(plotsdir, f"{titlename}_{analysis}"), dpi=300, vector=False)
        except Exception as e:
            print(f"Erro ao gerar o barplot de identificação: {str(e)}")
            sys.exit(1)

    elif analysis_type == "pca":
        try:
            print(f"Gerando {analysis}...")
            rawfiledata.pca(save=os.path.join(plotsdir, f"{titlename}_{analysis}"), dpi=300, vector=False)
        except Exception as e:
            print(f"Erro ao gerar o pca: {str(e)}")
            sys.exit(1)

    elif analysis_type == "volcano":
        try:
            print(f"Gerando {analysis}...")
            rawfiledata.volcano(*proteins, save=os.path.join(plotsdir, f"{titlename}_{analysis}"), dpi=300, vector=False)
        except Exception as e:
            print(f"Erro ao gerar o volcano: {str(e)}")
            sys.exit(1)


try:
    # Uso: python get_data.py caminho/para/a/planilha.xlsx metodo grupo_controle caminho/para/a/saida fc
    # Obtendo dados das análises a partir do usuario
    raw_file_path = sys.argv[1]
    sample_name = raw_file_path.split("/")[-1].split(".")[0]

    # Obtendo método de proteômica
    method = sys.argv[2]

    # Obtendo grupo controle
    control_group = sys.argv[3]

    # Obtendo outputh path
    output_path = sys.argv[4]
    plots_dir = os.path.join(output_path, "plots")

    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)

    # Obtendo foldchange
    fc = float(sys.argv[5])

    #Realizando a transformação do FC em Log2FC
    try:
        if fc <= 0:
            raise ValueError("FC deve ser maior que zero.")
        log2fc_cutoff = math.log2(fc)
    except ValueError as e:
        print(f"Erro: Valor de FC inválido. {str(e)}")
        sys.exit(1)

    # Lendo os dados e criando os arquivos com as planilhas
    # Obtendo os dados da proteômica a partir da planilha original
    print("Lendo o arquivo de dados...")
    raw_file_data = read_proteomics_file(raw_file_path, method, control_group, log2fc_cutoff)
    print("Arquivo de dados lido com sucesso.")

    # Tabela de parâmetros
    print("Obtendo parâmetros do estudo")
    print(pd.DataFrame(raw_file_data.Params))

    # Criando tabela com as DEPs
    print("Criando tabela de DEPs...")
    deps_dataframe = pd.DataFrame(raw_file_data.deps)
    deps_dataframe.to_excel(os.path.join(output_path, f"{sample_name}_deps.xlsx"))
except Exception as e:
    print(f"Erro durante a obtenção das DEPs: {e}")

analysis_types = ["volcano",
                  "ma_plot",
                  "correlation_heatmap",
                  "expression_heatmap",
                  "pca"]

# Filtragem das DEPs a partir do log2FC, ordenamento das UP e Down reguladas e obtenção das top 5
try:
    DEPs = []
    #obtendo as top 5 proteínas mais upreguladas
    upregulated_proteins_dataframe = (deps_dataframe[deps_dataframe['log2(fc)'] >= log2fc_cutoff]
                                      .sort_values(by=['log2(fc)'], ascending=False).head(5))
    upregulated_proteins_list = upregulated_proteins_dataframe['gene_name'].tolist()
    DEPs = DEPs + upregulated_proteins_list

    #obtendo os top5 proteínas mais downreguladas
    #ascending=True para obter os valores mais negativos no topo do dataframe, durante a filtragem
    downregulated_proteins_dataframe = (deps_dataframe[deps_dataframe['log2(fc)'] <= log2fc_cutoff]
                                        .sort_values(by=['log2(fc)'], ascending=True).head(5))
    downregulated_proteins_list = downregulated_proteins_dataframe['gene_name'].tolist()
    DEPs = DEPs + downregulated_proteins_list
except Exception as e:
    print(f"Erro ao obter as proteinas diferencialmente reguladas: {e}")

for condition in analysis_types:
    plot_data(
        rawfiledata=raw_file_data,
        analysis_type=condition,
        titlename=sample_name,
        plotsdir=plots_dir,
        analysis=condition,
        proteins=DEPs
    )
print("Plotagem dos dados concluída.")
