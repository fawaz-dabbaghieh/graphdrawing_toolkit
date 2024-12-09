from graphdraw.graph_io import read_gfa, write_gfa
from collections import deque
from graphdraw.Node import Node
import sys
import os


class Graph:
	"""
	Graph object containing the important information about the graph
	"""

	__slots__ = ['nodes', 'edge_counts']

	def __init__(self, graph_file=None, edge_count=False):
		if graph_file is not None:
			if not os.path.exists(graph_file):
				print("Error! Check log file.")
				logging.error("graph file {} does not exist".format(graph_file))
				sys.exit()
			# loading nodes from file
			self.nodes = read_gfa(gfa_file_path=graph_file)
			if edge_count:
				self.edge_counts = self.get_edges_counts(graph_file)
		else:
			self.nodes = dict()
			self.edge_counts = dict()

	def __len__(self):
		"""
		overloading the length function
		"""
		return len(self.nodes)

	def __str__(self):
		"""
		overloading the string function for printing
		"""
		return "The graph has {} Nodes".format(len(self.nodes))

	def __contains__(self, key):
		"""
		overloading the in operator to check if node exists in graph
		"""
		return key in self.nodes

	def __getitem__(self, key):
		"""
		overloading the bracket operator
		"""
		try:
			return self.nodes[key]
		except KeyError:
			return None

	def __setitem__(self, key, value):
		"""
		overloading setting an item in nodes
		"""
		if isinstance(value, Node):
			self.nodes[key] = value
		else:
			raise ValueError("the object given to set should be a Node object")

	def __delitem__(self, key):
		"""
		overloading deleting item
		"""
		del self.nodes[key]

# todo make remove start and remove end separate so I can use the same functions
#   for removing one edge
	def remove_node(self, n_id):
		"""
		remove a node and its corresponding edges
		"""
		starts = [x for x in self.nodes[n_id].start]
		for n_start in starts:
			overlap = n_start[2]
			if n_start[1] == 1:
				self.nodes[n_start[0]].end.remove((n_id, 0, overlap))
			else:
				self.nodes[n_start[0]].start.remove((n_id, 0, overlap))

		ends = [x for x in self.nodes[n_id].end]
		for n_end in ends:
			overlap = n_end[2]
			if n_end[1] == 1:
				self.nodes[n_end[0]].end.remove((n_id, 1, overlap))
			else:
				self.nodes[n_end[0]].start.remove((n_id, 1, overlap))

		del self.nodes[n_id]

	def remove_lonely_nodes(self):
		"""
		remove singular nodes with no neighbors
		"""
		nodes_to_remove = [n.id for n in self.nodes.values() if len(n.neighbors()) == 0]
		for i in nodes_to_remove:
			self.remove_node(i)


	def copy_edges(self, node1, node2, side):
		"""
		copy edges from one side and one node to another
		"""
		if side == 0:
			for e in self.nodes[node1].start:
				self.nodes[node2].start.append(e)
		elif side == 1:
			for e in self.nodes[node1].end:
				self[node2].end.append(e)
		else:
			print("The side has to be either 0 or 1 for start or end")


	def write_graph(self, set_of_nodes=None,
					output_file="output_graph.gfa",
					append=False):
		"""writes a graph file as GFA

		list_of_nodes can be a list of node ids to write
		ignore_nodes is a list of node ids to not write out
		if append is set to true then output file should be an existing
		graph file to append to
		modified to output a modified graph file
		"""
		if not output_file.endswith(".gfa"):
			output_file += ".gfa"
		# print("I am here")
		write_gfa(self, set_of_nodes=set_of_nodes, output_file=output_file,
				  append=append)


	def remove_edge(self, edge):
		n1, side1, n2, side2, overlap = edge
		if side1 == 0:
			self.nodes[n1].remove_from_start(n2, side2, overlap)
		else:
			self.nodes[n1].remove_from_end(n2, side2, overlap)

		if side2 == 0:
			self.nodes[n2].remove_from_start(n1, side1, overlap)
		else:
			self.nodes[n2].remove_from_end(n1, side1, overlap)

	# write BFS here, no need for doing the stupid two direction stuff
	# fix it in both GFASubgraph and BubbleGun
	def reset_visited(self):
		for n in self.nodes.values():
			n.visited = False

	def bfs(self, start_id, size):
		self.reset_visited()

		if len(self.nodes[start_id].neighbors()) == 0:
			return {start_id}

		queue = deque()
		queue.append(start_id)
		self.nodes[start_id].visited = True
		neighborhood = set()
		while (len(neighborhood) <= size) and len(queue) > 0:
			s = queue.popleft()
			neighborhood.add(s)
			self.nodes[s].visited = True
			for n in self.nodes[s].neighbors():
				if not self.nodes[n].visited:
					queue.append(n)
					self.nodes[n].visited = True
		return neighborhood
