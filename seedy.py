import heapq
import pandas as pd

dados = pd.read_csv('dados.csv')
dados = dados.drop_duplicates()
sementes = ['TERT', 'PTEN', 'EGFR']
cluster = {k: 1 for k in sementes}
threshold = 0.8

fila_prioridade = []   
scores_na_fila = {}

def atualizar_fila(proteina, score):
    if proteina in cluster:
        return
    if proteina in scores_na_fila:
        if scores_na_fila[proteina] >= score:
            return  
        scores_na_fila[proteina] = score
        heapq.heappush(fila_prioridade, (-score, proteina))
    else:
        scores_na_fila[proteina] = score
        heapq.heappush(fila_prioridade, (-score, proteina))

for gene in sementes:
    mask = ((dados['preferredName_A'] == gene) | 
            (dados['preferredName_B'] == gene)) & \
           (dados['score'] >= threshold)
    
    sub = dados[mask]
    for _, row in sub.iterrows():
        parceiro = row['preferredName_B'] if row['preferredName_A'] == gene else row['preferredName_A']
        score = row['score']
        atualizar_fila(parceiro, score)

while fila_prioridade:
    neg_score, proteina = heapq.heappop(fila_prioridade)
    maxScore = -neg_score

    if scores_na_fila.get(proteina) != maxScore:
        continue

    del scores_na_fila[proteina]
    cluster[proteina] = maxScore

    mask = ((dados['preferredName_A'] == proteina) |
            (dados['preferredName_B'] == proteina)) & \
           (dados['score'] >= threshold)

    for _, row in dados[mask].iterrows():
        vizinho = row['preferredName_B'] if row['preferredName_A'] == proteina else row['preferredName_A']
        thisScore = maxScore * row['score']

        if thisScore < threshold:
            continue

        atualizar_fila(vizinho, thisScore)

print(f'Total de proteínas no cluster: {len(cluster)}')

arestas_finais = dados[
    (dados['preferredName_A'].isin(cluster.keys())) &
    (dados['preferredName_B'].isin(cluster.keys())) &
    (dados['score'] >= threshold)
].copy()

arestas_finais['score_A'] = arestas_finais['preferredName_A'].map(cluster)
arestas_finais['score_B'] = arestas_finais['preferredName_B'].map(cluster)

arestas_finais.to_csv('rede_ppi_arestas.csv', sep=',', index=False)