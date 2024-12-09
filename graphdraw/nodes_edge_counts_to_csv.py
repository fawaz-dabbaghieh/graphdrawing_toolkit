import os
import sys
import pickle
from collections import defaultdict
from BubbleGun.Graph import Graph


def int_list():
	return [0,0]


def add_len_to_node(graph, n):
	graph.nodes[n].seq_len = int(graph.nodes[n].optional_info.split("\t")[0].split(":")[-1])


def count_nodes_edges(gaf, edges, nodes, index):
	#  index is 0 or 1, 0 for cancer and 1 for healthy
	#  nodes here is a dict of node IDs and an int counter for number of
	#  edges is a dict with (n1, n2) as keys and n1 > n2
	#  reads in that gaf that touched that node

	with open(gaf, "r") as infile:
		for l in infile:
			l = l.strip().split()
			seq_name = l[0]
			path = l[5]  # the alignment path in the graph
			path = path.replace("<", ",").replace(">", ",")[1:].split(",")
			# node counts
			for s in path:
				nodes[s][index] += 1
			# edge counts
			for i in range(len(path) - 1):
				if path[i] > path[i+1]:
					edges[(path[i], path[i+1])][index] += 1
				else:
					edges[(path[i+1], path[i])][index] += 1


if len(sys.argv) < 5:
	print("You need to give the graph cancer.gaf healthy.gaf and an output prefix")
	sys.exit()


print("loading graph")
graph = Graph(sys.argv[1])
for n in graph.nodes.keys():
	add_len_to_node(graph, n)


nodes = defaultdict(int_list)
edges = defaultdict(int_list)
# cancer counts
print("counting cancer reads")
count_nodes_edges(sys.argv[2], edges, nodes, 0)
print("counting healthy reads")
count_nodes_edges(sys.argv[3], edges, nodes, 1)


prefix = sys.argv[4]
print(f"outputting a table for node counts")
with open(f"{prefix}_node_counts.csv", "w") as outfile:
	outfile.write("node_id,cancer_counts,healthy_counts,node_length\n")
	for n, counts in nodes.items():
		outfile.write(f"{n},{counts[0]},{counts[1]},{graph.nodes[n].seq_len}\n")

print(f"outputting a table for edge counts")
with open(f"{prefix}_edge_counts.csv", "w") as outfile:
	outfile.write("edge_id,cancer_counts,healthy_counts\n")
	for e, counts in edges.items():
		outfile.write(f"{e[0]}_{e[1]},{counts[0]},{counts[1]}\n")
