import os
import sys
import pdb
from BubbleGun.Graph import Graph


def write_gfa(graph, output_file="output_file.gfa"):
	"""
	this one will write a graph but with new node IDs
	:param graph: the graph object
	:param output_file: path to output file
	:return: writes a gfa file
	"""
	nodes = graph.nodes

	if os.path.exists(output_file):
		print(f"file {output_file} already exists, rewriting")
		
	f = open(output_file, "w+")

	for n1 in graph.nodes.values():

		line = str("\t".join(("S", str(n1.id), n1.seq)))
		line += "\t" + n1.optional_info  # adding optional info

		f.write(line + "\n")

		# writing edges
		edges = []
		# overlap = str(graph.k - 1) + "M\n"

		for n in n1.start:
			overlap = str(n[2]) + "M\n"
			# I am checking if the are nodes I want to write
			# I think I can remove this later as I implemented the .remove_node
			# to the Graph class that safely removes a node and all its edges
			# So there shouldn't be any edges to removed
			if n[1] == 0:
				edge = str("\t".join(("L", str(n1.id), "-", graph.nodes[n[0]].id, "+", overlap)))
				edges.append(edge)
			else:
				edge = str("\t".join(("L", str(n1.id), "-", graph.nodes[n[0]].id, "-", overlap)))
				edges.append(edge)

		for n in n1.end:
			overlap = str(n[2]) + "M\n"

			if n[1] == 0:
				edge = str("\t".join(("L", str(n1.id), "+", graph.nodes[n[0]].id, "+", overlap)))
				edges.append(edge)
			else:
				edge = str("\t".join(("L", str(n1.id), "+", graph.nodes[n[0]].id, "-", overlap)))
				edges.append(edge)

		for e in edges:
			f.write(e)

	f.close()


if len(sys.argv) < 4:
	print("you need to give input gfa, node suffix, output gfa")
	sys.exit()

in_gfa = sys.argv[1]
suffix = sys.argv[2]
if not suffix.startswith("_"):
	suffix = "_" + suffix
out_gfa = sys.argv[3]

graph = Graph(in_gfa)
for n in graph.nodes.values():
	n.id = n.id + suffix

write_gfa(graph, out_gfa)
