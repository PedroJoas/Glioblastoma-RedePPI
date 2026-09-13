import requests 
from collections import deque
import pandas as pd
from time import sleep

string_api_url = "https://version-12-0.string-db.org/api"
output_format = "tsv"
method = "network"

##
## Construct URL
##

request_url = "/".join([string_api_url, output_format, method])

##
## Set parameters
##

fila = deque(['TERT', 'PTEN', 'EGFR'])
visitados = set(fila)
on_off = 0 

while fila:
    gene = fila.popleft()
    print(f'--------- {gene} ---------', end='\r')
    params = {
        "identifiers" : gene, # your protein
        "species" : 9606 # NCBI/STRING taxon identifier 
    }

    ##
    ## Call STRING
    ##

    response = requests.post(request_url, data=params)
    with open('dados.csv', 'a') as f:
        f.write(response.text.strip().replace('	', ',')+'\n')

    sleep(1)

    if on_off == 0:
        request_url = request_url.replace('tsv', 'tsv-no-header')
        on_off = 1

    dados = pd.read_csv('dados.csv')
    novos = list(dados['preferredName_A'].unique())
    print(f'Quantidade de genes único: {len(novos)}', end='\r')
    
    if len(novos) >= 400:
        break
    for novo in novos:
        if novo not in visitados:
            fila.append(novo)
        visitados.add(novo)

    