import requests
from collections import deque
import pandas as pd
from time import sleep

string_api_url = "https://version-12-0.string-db.org/api"
output_format = "tsv-no-header"
method = "network"
request_url = "/".join([string_api_url, output_format, method])

# Parâmetros
SEMENTES = ['TERT', 'PTEN', 'EGFR']
SCORE_MIN = 0.7          # filtro de score
MAX_GRAU = 4          # profundidade máxima a partir das sementes
MAX_NOS = 3000           # limite de segurança (opcional)

# BFS com controle de profundidade
fila = deque([(s, 0) for s in SEMENTES])   # (gene, profundidade)
visitados = set(SEMENTES)
interacoes = []                             # acumula arestas

while fila and len(visitados) < MAX_NOS:
    gene, grau = fila.popleft()
    
    if grau > MAX_GRAU:
        continue
    
    print(f'--------- {gene} (grau {grau}) ---------', end='\r')
    
    params = {"identifiers": gene, "species": 9606}
    response = requests.post(request_url, data=params)
    sleep(1)
    
    # Processa resposta
    linhas = [l.split('\t') for l in response.text.strip().split('\n') if l]
    if not linhas:
        continue
    
    df = pd.DataFrame(linhas)
    # STRING tsv-no-header: colunas típicas incluem A, B, score...
    # Ajuste os índices conforme o formato real
    # Supondo: 0=stringId_A, 1=stringId_B, 2=preferredName_A, 3=preferredName_B, ..., score na última coluna
    df = df.rename(columns={2: 'A', 3: 'B', df.columns[-1]: 'score'})
    df['score'] = pd.to_numeric(df['score'], errors='coerce')
    df = df[df['score'] >= SCORE_MIN]
    
    for _, row in df.iterrows():
        a, b = row['A'], row['B']
        interacoes.append((a, b, row['score']))
        
        # Expande BFS
        if grau < MAX_GRAU:
            for parceiro in (a, b):
                if parceiro not in visitados:
                    visitados.add(parceiro)
                    fila.append((parceiro, grau + 1))

# Salva arestas únicas
df_inter = pd.DataFrame(interacoes, columns=['A', 'B', 'score']).drop_duplicates()
df_inter.to_csv('interacoes_filtradas.csv', index=False)

print(f'\nNós: {len(visitados)} | Arestas: {len(df_inter)}')