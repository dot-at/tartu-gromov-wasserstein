import networkx as nx
import numpy as np
from numba import jit

@jit()
def rand(n):
    a = np.random.random(n)
    my_sum = 0
    for i in a:
        my_sum += i;
    for i in range(n):
        a[i] /= my_sum
    return a

@jit()
def export(name,g,Distance,Distribution="uniform"):
    f = open(name,"w")
    tmp = "n " + str(len(g)) + "\n"
    f.writelines(tmp)
    if Distribution == "uniform":
        tmp = "u\n"
        f.writelines(tmp)
    else:
        tmp = "g "
        ran = rand(len(g))
        for i in ran:
            tmp += str(i) + " "
    tmp += "\n"
    f.writelines(tmp)
    for i in range(1,len(g)):
        tmp = "d "
        for j in range(i):
            tmp += str(Distance[g[i]][g[j]]) + " "
        tmp += "\n"
        f.writelines(tmp)
#######################################################
#
# Uncomment the line you need and input the desired parameters
#
#######################################################
# g = nx.fast_gnp_random_graph(n, p, seed=None, directed=False)
# g = nx.gnp_random_graph(n, p, seed=None, directed=False)
# g = nx.dense_gnm_random_graph(n, m, seed=None)
# g = nx.gnm_random_graph(n, m, seed=None, directed=False)
# g = nx.erdos_renyi_graph(n, p, seed=None, directed=False)
# g = nx.binomial_graph(n, p, seed=None, directed=False)
# g = nx.newman_watts_strogatz_graph(n, k, p, seed=None)
# g = nx.watts_strogatz_graph(n, k, p, seed=None)
# g = nx.connected_watts_strogatz_graph(n, k, p, tries=100, seed=None)
# g = nx.random_regular_graph(d, n, seed=None)
# g = nx.barabasi_albert_graph(n, m, seed=None)
# g = nx.powerlaw_cluster_graph(n, m, p, seed=None)
# g = nx.random_lobster(n, p1, p2, seed=None)
# g = nx.random_powerlaw_tree(n, gamma=3, seed=None, tries=100)

D = nx.all_pairs_shortest_path_length(g)
C = nx.connected_components(g)
#################################
#
# To output the results in another folder change path variable
#
#################################
path = ""
name_index = 0
#################################
#
# To discard small graphs change treshold
#
#################################
treshold = 0
for i in C:
    if len(i) > treshold:
        export(path + str(name_index),sorted(i),D,"general")
        name_index += 1
