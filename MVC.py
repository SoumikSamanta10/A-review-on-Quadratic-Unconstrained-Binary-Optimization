# import structure for Q matrix
from collections import defaultdict

# Import networkx for graph tools
import networkx as nx

# Import dwave_networkx for d-wave graph tools/functions
import dwave_networkx as dnx

# Import matplotlib.pyplot to draw graphs on screen
import matplotlib
matplotlib.use("agg")    # must select backend before importing pyplot
import matplotlib.pyplot as plt

# Set the solver we're going to use
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

sampler = EmbeddingComposite(DWaveSampler())

# Create empty graph
G = nx.Graph()

# Add edges to graph - this also adds the nodes
G.add_edges_from([(1,2),(1,7),(2,3),(2,7),(3,4),(3,6),(4,6),(4,5),(4,7),(5,7),(6,7)])

#Formulate our problem
Q = defaultdict(int)
P = 8       #Penalty term

# Update Q matrix for every edge in the graph
# We try to minimize (1-yP)x_{i} + P*x_{i}*x{j}
# We subtract P for a node everytime it comes as part of an edge
# and add P to all the matrix elements representing edges
for i in G.nodes:
    Q[(i,i)] = 1

for i, j in G.edges:
    #print("{}   {}".format(i,j))
    Q[(i,i)]-= P
    Q[(j,j)]-= P
    Q[(i,j)]+= P/2
    Q[(j,i)]+= P/2

#print("The matrix after assignment is: {}".format(Q))

chainstrength = 8
numruns = 10

sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type': 'chimera'}))
S = sampler.sample_qubo(Q,num_reads=numruns)

# ------- Print results to user -------
print('-' * 60)
print('{:>15s}{:>15s}{:^15s}'.format('Not Cover','Cover','Energy'))
print('-' * 60)

for sample, E in S.data(fields=['sample','energy']):
    S0 = [k for k,v in sample.items() if v == 0]        #first set (0 nodes)
    S1 = [k for k,v in sample.items() if v == 1]        #second set (1 nodes)
    print('{:>15s}{:>15s}{:^15s}'.format(str(S0),str(S1),str(E)))


# THE NEXT PART IS DISPLAYING THE BEST RESULT GRAPHICALLY
# ------- Display results to user -------
# Grab best result
# Note: "best" result is the result with the lowest energy
# Note2: the look up table (lut) is a dictionary, where the key is the node index
#   and the value is the set label. For example, lut[5] = 1, indicates that
#   node 5 is in set 1 (S1).
lut = S.first.sample

# Interpret best result in terms of nodes and edges
S0 = [node for node in G.nodes if not lut[node]]            #storing the first set
S1 = [node for node in G.nodes if lut[node]]                #storing the second set
cut_edges = [(u, v) for u, v in G.edges if lut[u]!=lut[v]]  #storing edges that were cut (between S1 and S0)  
uncut_edges = [(u, v) for u, v in G.edges if lut[u]==lut[v]]#storing edges that weren't cut (within each set)

# Print the solution for the user
print('Minimum vertex cover found is', len(S1))
print(S1)

# Saving the initial graph
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, nodelist=G.nodes, node_color='r')
nx.draw_networkx_edges(G, pos, edgelist=G.edges, style='solid', width=3)
nx.draw_networkx_labels(G, pos)

filename = "pipelines_original.png"
plt.savefig(filename, dpi=600)
print("\nOriginal plot has been saved to {}". format(filename))

plt.clf()   #clear old graph for new one

# Display best result
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, nodelist=S0, node_color='r')
nx.draw_networkx_nodes(G, pos, nodelist=S1, node_color='c')
nx.draw_networkx_edges(G, pos, edgelist=cut_edges, style='solid', alpha=0.5, width=3)
nx.draw_networkx_edges(G, pos, edgelist=uncut_edges, style='solid', width=3)
nx.draw_networkx_labels(G, pos)

# Saving the best result as a png
filename = "pipelines.png"
plt.savefig(filename, dpi=600)
print("\nYour plot is saved to {}".format(filename))

print(S)
