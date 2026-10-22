import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

dados = pd.read_csv('interacoes_filtradas.csv')

print(dados.shape)

G = nx.Graph()
print('Adicionando genes')
for i, row in dados.iterrows():
    gene_a = row['A']
    gene_b = row['B']
    score = row['score']

    if not G.has_edge(gene_a, gene_b) and gene_a != gene_b and score >= 0.8:
        G.add_edge(gene_a, gene_b, score=score)

print('Salvando edgelist')
nx.write_edgelist(G, "graph.edgelist")

nx.draw(G)
print('Salvando grafo')
plt.savefig("graph_plot.png", format="PNG", dpi=300)

plt.show()