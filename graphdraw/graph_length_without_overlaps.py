import sys
from BubbleGun.Graph import Graph


def return_tuple(s1, s2, overlap):
	if s1>s2:
		return (s1, s2, overlap)
	return (s2, s1, overlap)


if len(sys.argv) < 2:
	print("Please give an input GFA")

graph = Graph(sys.argv[1])

edges_counted = set()
total_overlaps = 0
total_length = 0
n_edges = 0

for n in graph.nodes.values():
	total_length += n.seq_len

	for nn in n.start:
		n_edges += 1
		overlap = nn[2]
		pair = return_tuple(n.id, nn[0], overlap)
		if pair in edges_counted:
			continue
		else:
			edges_counted.add(pair)
			total_overlaps += overlap

	for nn in n.end:
		n_edges += 1
		overlap = nn[2]
		pair = return_tuple(n.id, nn[0], overlap)
		if pair in edges_counted:
			continue
		else:
			edges_counted.add(pair)
			total_overlaps += overlap


print(f"number of nodes: {len(graph)}")
print(f"number of edges: {int(n_edges/2)}")
print(f"total graph length: {total_length}")
print(f"total graph length without overlaps: {total_length - total_overlaps}")
