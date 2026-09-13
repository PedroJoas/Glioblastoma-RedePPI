import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

dados = pd.read_csv('rede_ppi_arestas.csv')
dados_limpo = dados.drop_duplicates()

print(dados.shape)
print(dados_limpo.shape)

G = nx.Graph()
print('Adicionando genes')
for i, row in dados_limpo.iterrows():
    gene_a = row['preferredName_A']
    gene_b = row['preferredName_B']
    score = row['score']

    if not G.has_edge(gene_a, gene_b) and gene_a != gene_b:
        G.add_edge(gene_a, gene_b, score=score)

print('Salvando edgelist')
nx.write_edgelist(G, "graph.edgelist")

pos = nx.kamada_kawai_layout(G)
nx.draw(G, pos=pos)
print('Salvando grafo')
plt.savefig("graph_plot.png", format="PNG", dpi=300)

plt.show()