import os
import json

#script para criação dos diretórios necessários para a ferramenta

#verifica se o diretório data existe

if not os.path.exists('data'):
    print("Criando diretórios de dados...", end="\r")
    os.mkdir('data')
    os.mkdir('data/input')
    os.mkdir('data/output')
    os.mkdir('data/temp')

if os.path.exists('data'):
    if not os.path.exists('data/input'):
        print("Criando diretório de input...", end="\r")
        os.mkdir('data/input')
    if not os.path.exists('data/output'):
        print("Criando diretório de output...", end="\r")
        os.mkdir('data/output')

if not os.path.exists('config.json'):
    config_data = {
        "FC": 1.2,
        "Method": "Progenesis",
        "Control": "Nome_controle"
    }

    with open('config.json', 'w') as outfile:
        json.dump(config_data, outfile)

    print("Arquivo de configuração config.json criado com sucesso! Verfique-o para adaptar as"
          "suas preferências de análise.")

print("Configuração concluída. Insira a(s) planilha(s) no diretório de input e execute o get_data.py para obter as DEPs do estudo de interesse")
print("Execução do get_data.py:")
print("python3 get_data.py caminho/para/a/planilha.xlsx")
