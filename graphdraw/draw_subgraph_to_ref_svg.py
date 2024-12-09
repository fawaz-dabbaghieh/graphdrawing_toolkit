import sys
import os
import logging
import random
import pysam
from graphdraw.color_graph_with_chromosomes import Alignment
from graphdraw.Graph import Graph
from graphdraw.SVG import *
from graphdraw.extract_subgraph_from_ref_coord import extract_subgraph

"""
I need to use this to find the q point for a quadratic curve
dx = x1 - x2
midpoint = ((x1 + x2) / 2, (y1 + y2) / 2)
slope = -dx / (y1 - y2)
x = sqrt(altitude*altitude - dx*dx) / slope + midpoint.x
y = slope * (x - midpoint.x) + midpoint.y

To calculate the middle point for the quadratic curve of an edge I can do the following
I have the two points I am connecting (x1, y1) (x2, y2)
I get the mid point (x1+x2/2, y1+y2/2)
calculate the norm
xn, yn = -(y2-y1), x2-x1
unit = sqrt(xn^2 + yn^2)
the unit norm is (xn/unit, yn/unit) = ()

then to get the new point from the mid point I just do
mid_point + h*
"""


def sys_exit():
	print("Error Happened! Please check the log file")
	sys.exit()


def get_nodes_info(in_bam, coordinates):
	"""
	Takes a bam file and uses the coordinates to fetch the alignments and filters them based on top score
	Uses then these alignments for drawing
	"""
	try:
		if in_bam.endswith(".bam"):
			alignment_file = pysam.AlignmentFile(in_bam, "rb")
		else:
			alignment_file = pysam.AlignmentFile(in_bam, "r")
	except FileNotFoundError:
		logging.error(f"The file {in_bam} does not exist")
		sys_exit()
	except ValueError:
		logging.error(f"The file {in_bam} is not a bam file")
		sys_exit()

	chrom, coord = coordinates.split(":")
	start, end = coord.split("-")
	all_alignments = dict()
	# I probably need to change this filtering later and do a better one
	for a in alignment_file.fetch(chrom, int(start), int(end)):
		if a.flag in {0, 16}:
			if a.qname in all_alignments:
				all_alignments[a.qname] = compare(all_alignments[a.qname], a)
			else:
				all_alignments[a.qname] = Alignment(a)
	nodes = dict()
	for a in all_alignments.values():
		nodes[a.q_name] = {"start":int(a.ref_start), "end":int(a.ref_end), "orientation":a.orientation, "chromosome":a.ref_name, "mapq":a.mapq}
	return nodes


def get_ref_info(nodes, chromosome, legend):
	"""
	looks at all the alignments provided and takes the smallest start and largest end as the reference to draw against
	"""
	# the smallest start
	start = min([x["start"] for x in nodes.values()])
	end = max([x["end"] for x in nodes.values()])
	assert end > start

	# 10 px for padding
	reference = {"start": start, "end": end, "length": end - start, "start_px": legend, "color": "Black", "thick": 3, "chromosome": chromosome}
	return reference


def order_nodes(nodes, start=True):
	"""
	Orders nodes by start or end
	"""
	if start:
		ordered_keys = sorted(nodes, key=lambda k: nodes[k]["start"])
	else:
		ordered_keys = sorted(nodes, key=lambda k: nodes[k]["end"])
	return list(ordered_keys)


def calculate_opacity(mapq):
	"""
	maps the mapq value of the alignment to the opacity
	"""
	assert mapq <= 60
	assert mapq >= 0
	rr = 80/60
	return ((mapq)*rr + 20)/100


def check_for_clash(nodes, wiggle_room):
	"""
	I will brute force the check because the number of nodes is small enough
	"""
	row_intervals = dict()
	for n in nodes.keys():
		if nodes[n]['y'] not in row_intervals:
			row_intervals[nodes[n]['y']] = [(nodes[n]['x1'], nodes[n]['x2'], n)]
		else:
			row_intervals[nodes[n]['y']].append((nodes[n]['x1'], nodes[n]['x2'], n))

	# brute force for each node in a row
	plus_or_minus = [0, 1]
	for row in row_intervals.values():
		for x1, x2, n in row:
			for nx1, nx2, nn in row:
				if n != nn:
					if (nx1 <= x1 <= nx2) or (nx1 <= x2 <= nx2):
						if random.choice(plus_or_minus) == 0:
							nodes[n]['y'] = nodes[n]['y'] + wiggle_room
						else:
							nodes[n]['y'] = nodes[n]['y'] - wiggle_room


def calculate_nodes_coord(nodes, ordered_nodes, reference, output_svg):
	"""
	calculates the x and y coordinates of each alignment, to know where to draw each node
	compared to the reference

	also resolves clashes in case one node will be drawn over another one 
	"""
	# figure out length, start and end of each node based on their alignment to the reference in pixels
	x_coord = []
	for idx, n in enumerate(ordered_nodes):
		length = nodes[n]["end"] - nodes[n]["start"]
		len_perc = length/reference["length"]  # how much does this node cover of the reference
		px_len = reference['end_px'] - reference['start_px']
		n_len_in_px = round(len_perc*px_len)  # length of the node in pixels compared to the reference
		n_start = round(((nodes[n]["start"] - reference["start"])/reference["length"])*px_len)  # where to start
		nodes[n]["x1"] = n_start + output_svg.legend
		nodes[n]["x2"] = n_start + n_len_in_px + output_svg.legend
		# x_coord.append((n_start, n_len_in_px))  # where to start and how long
	#### TODO: i need to cap the changes in y between 25-100% of the height
	# and always check the new y if it is in that constraint, otherwise, either loop from the bottom
	# or just push it a little down

	# I draw the graph after keeping 20% on top of the image for the reference
	# and 10 pixels at the bottom
	y_upper_border = int(output_svg.height * 0.20)
	y_lower_border = output_svg.height - 10
	y_range = y_lower_border - y_upper_border
	n_rows = 5  # how many rows to spread the nodes over

	allowed_y = [0] * n_rows  # I am only allowing 4 y coordinates
	nodes_in_row = []  # I need this later to alternate the text by top and bottom of nodes to prevent clashes
	for _ in range(n_rows):
		nodes_in_row.append(list())

	for i in range(0, n_rows):
		allowed_y[i] = int(y_range/n_rows)*i + y_upper_border

	# drawing nodes
	previous = 0
	for idx, n in enumerate(ordered_nodes):
		nodes[n]['y'] = allowed_y[idx%n_rows]
		# nodes[n]['y'] = random.choice(allowed_y)
		nodes_in_row[allowed_y.index(nodes[n]['y'])].append(n)
		# pdb.set_trace()
		# nodes_in_row[idx%n_rows].append(n)

		style = dict()
		style["stroke"] = "Black"
		style["stroke-width"] = 3
		# random_mapq = random.choice(range(0, 60))
		# style["stroke-opacity"] = str(calculate_opacity(random_mapq))

		style["stroke-opacity"] = calculate_opacity(int(nodes[n]['mapq'])*0.4)  # opacity based on mapping quality
		nodes[n]['style'] = style

	# checking for clashes and changing the coordinates of the node
	# just chose to run it 3 time to get rid of most of the clashes, some might still persist
	# pdb.set_trace()
	check_for_clash(nodes, (y_range/n_rows)//4)
	check_for_clash(nodes, (y_range/n_rows)//3)
	check_for_clash(nodes, (y_range/n_rows)//3)

	# add nodes to SVG
	for n in nodes:
		# pdb.set_trace()
		rotation_value = 0
		rotate_node = (nodes[n]['x2'] - nodes[n]['x1'])*rotation_value
		if nodes[n]['orientation'] == "-":
			nodes[n]['style']['stroke'] = "Red"
			nodes[n]['y1'] = nodes[n]['y'] + rotate_node
			nodes[n]['y2'] = nodes[n]['y'] - rotate_node
		else:
			nodes[n]['y1'] = nodes[n]['y'] - rotate_node
			nodes[n]['y2'] = nodes[n]['y'] + rotate_node
		output_svg.add_line(nodes[n]['x1'], nodes[n]['y1'], nodes[n]['x2'], nodes[n]['y2'], nodes[n]['style'])

	"""
	generate style for the text as well and pass it to add_txt
	"""
	for row in nodes_in_row:
		for idx, n in enumerate(row):
			mid_point = nodes[n]['x1'] + (nodes[n]['x2'] - nodes[n]['x1'])//2
			if idx%2 == 0:  # even
				output_svg.add_txt(mid_point, nodes[n]['y'] + 1, n, "Black", 3)
			else:
				output_svg.add_txt(mid_point, nodes[n]['y'] + 1, n, "Black", 3)


def draw_subgraph(in_bam, in_graph, coordinates, svg_out, svg_width, svg_height, svg_legend):
	"""
	1- I think I need to go through the TSV and get all the nodes that aligned in the coordinates given
	then plot these nodes in that coordinates and connect any available edges
	problem here is that these nodes might not have edges

	2- Second way of doing this is by extracting a sub-graph, then looking at the TSV table for the alignment of each 
	node, then we plot
	problem here is that the alignment might not be in the same area or even on different chromosomes because of repeats and so on
	"""

	"""
	TODO integrate the graph subgraph extraction here so I don't have to run it separately before

	1- First a function that filters the BAM
	2- Then using the nodes that come back from this, run BFS
	3- now I have the ndoes dict and the graph and the rest can run as is
	"""

	nodes = get_nodes_info(in_bam, coordinates)
	chrom, coord = coordinates.split(":")
	start, end = coord.split("-")
	graph = Graph(in_graph)
	nodes_to_remove = []
	for n in graph.nodes.keys():
		if n not in nodes:
			nodes_to_remove.append(n)
	for n in nodes_to_remove:
		graph.remove_node(n)
	# pdb.set_trace()
	# extract_subgraph(in_graph, list(nodes.keys()), in_bam, chrom, int(start), int(end), 0, "tmp.gfa")
	# graph = Graph("tmp.gfa")
	assert set(nodes.keys()) == set(graph.nodes.keys())
	# graph = Graph(in_graph)
	out_file = svg_out

	chromosome = list(set([x['chromosome'] for x in nodes.values()]))
	assert len(chromosome) == 1

	reference = get_ref_info(nodes, chromosome[0], svg_legend)
	output_svg = SVG(svg_height, svg_width, svg_legend)
	reference['end_px'] = output_svg.width - svg_legend  # offset on the end by 10 pixels
	output_svg.add_ref(reference)
	ordered_nodes = order_nodes(nodes)

	# calculates the x and y for the start and end of each node and adds this information to nodes
	# also draws the lines into output_svg
	calculate_nodes_coord(nodes, ordered_nodes, reference, output_svg)

	# draw edges
	style = dict()
	style["stroke"] = "Gray"
	style["stroke-width"] = 0.2
	style["stroke-opacity"] = 1
	drawn_edges = set()
	# for n_id, n in graph.nodes.items():
	for n_id in nodes.keys():
		n = graph.nodes[n_id]
		# I need to consider alignment orientation here
		for nn in n.start:
			nn_id = nn[0]  # neighbors id
			if n_id > nn_id:
				pair = (n_id, 0, nn_id, nn[1])
			else:
				pair = (nn_id, nn[1], n_id, 0)
			if pair not in drawn_edges:
				drawn_edges.add(pair)
			else:
				continue

			n1_len = nodes[n_id]['x2'] - nodes[n_id]['x1']
			n2_len = nodes[nn_id]['x2'] - nodes[nn_id]['x1']

			"""
			move all this big mess into a separate function
			calculates all probabilities then maybe another function for the bottom part
			that assigns which curve to add
			I also need to avoid drawing edges twice, I can simply keep a set of drawn edges and avoid the ones already there
			"""
			start_start = [nodes[n_id]['x1'], nodes[n_id]['y1'], nodes[nn_id]['x1'], nodes[nn_id]['y1'], n1_len, n2_len, 0, 0, style]
			start_end = [nodes[n_id]['x1'], nodes[n_id]['y1'], nodes[nn_id]['x2'], nodes[nn_id]['y2'], n1_len, n2_len, 0, 1, style]
			end_start = [nodes[n_id]['x2'], nodes[n_id]['y2'], nodes[nn_id]['x1'], nodes[nn_id]['y1'], n1_len, n2_len, 1, 0, style]
			end_end = [nodes[n_id]['x2'], nodes[n_id]['y2'], nodes[nn_id]['x2'], nodes[nn_id]['y2'], n1_len, n2_len, 1, 1, style]

			if nn[1] == 0:  # nn connects to n from the start
				# start to start connection
				if nodes[n_id]['orientation'] == "+" and nodes[nn_id]['orientation'] == "+": # stays start-start
					output_svg.add_curve(*start_start)
				elif nodes[n_id]['orientation'] == "-" and nodes[nn_id]['orientation'] == "+": # becomes end-start
					output_svg.add_curve(*end_start)
				elif nodes[n_id]['orientation'] == "+" and nodes[nn_id]['orientation'] == "-": # becomes start-end
					output_svg.add_curve(*start_end)
				else:  # - - so it becomes end-end
					output_svg.add_curve(*end_end)

			else:  # nn connects to n from the end
				# start to end
				if nodes[n_id]['orientation'] == "+" and nodes[nn_id]['orientation'] == "+": # stays start-end
					output_svg.add_curve(*start_end)
				elif nodes[n_id]['orientation'] == "-" and nodes[nn_id]['orientation'] == "+": # becomes end-end
					output_svg.add_curve(*end_end)
				elif nodes[n_id]['orientation'] == "+" and nodes[nn_id]['orientation'] == "-": # becomes start-start
					output_svg.add_curve(*start_start)
				else:  # - - so it becomes end-start
					output_svg.add_curve(*end_start)

		for nn in n.end:
			nn_id = nn[0]  # neighbors id
			if n_id > nn_id:
				pair = (n_id, 1, nn_id, nn[1])
			else:
				pair = (nn_id, nn[1], n_id, 1)
			if pair not in drawn_edges:
				drawn_edges.add(pair)
			else:
				continue

			n1_len = nodes[n_id]['x1'] - nodes[n_id]['x1']
			n2_len = nodes[nn_id]['x2'] - nodes[nn_id]['x1']

			start_start = [nodes[n_id]['x1'], nodes[n_id]['y1'], nodes[nn_id]['x1'], nodes[nn_id]['y1'], n1_len,
				n2_len, 0, 0, style]
			start_end = [nodes[n_id]['x1'], nodes[n_id]['y1'], nodes[nn_id]['x2'], nodes[nn_id]['y2'], n1_len, n2_len, 0, 1, style]
			end_start = [nodes[n_id]['x2'], nodes[n_id]['y2'], nodes[nn_id]['x1'], nodes[nn_id]['y1'], n1_len, n2_len, 1, 0, style]
			end_end = [nodes[n_id]['x2'], nodes[n_id]['y2'], nodes[nn_id]['x2'], nodes[nn_id]['y2'], n1_len, n2_len, 1, 1, style]

			if nn[1] == 0:
				# from end to start

				if nodes[n_id]['orientation'] == "+" and nodes[nn_id]['orientation'] == "+": # stays end-start
					output_svg.add_curve(*end_start)
				elif nodes[n_id]['orientation'] == "-" and nodes[nn_id]['orientation'] == "+": # becomes start-start
					output_svg.add_curve(*start_start)
				elif nodes[n_id]['orientation'] == "+" and nodes[nn_id]['orientation'] == "-": # becomes end-end
					output_svg.add_curve(*end_end)
				else:  # - - so it becomes start-end
					output_svg.add_curve(*start_end)

			else:
				# from end to end
				if nodes[n_id]['orientation'] == "+" and nodes[nn_id]['orientation'] == "+": # stays end-end
					output_svg.add_curve(*end_end)
				elif nodes[n_id]['orientation'] == "-" and nodes[nn_id]['orientation'] == "+": # becomes start-end
					output_svg.add_curve(*start_end)
				elif nodes[n_id]['orientation'] == "+" and nodes[nn_id]['orientation'] == "-": # becomes end-start
					output_svg.add_curve(*end_start)
				else:  # - - so it becomes start_start
					output_svg.add_curve(*start_start)

	output_svg.out_svg(out_file)
	# os.remove("tmp.gfa")
