import os
import sys
import pdb


def process_gaf(gaf):
	info = dict()
	with open(gaf, "r") as infile:
		for l in infile:
			l = l.strip().split()
			seq_name = l[0]
			path = l[5]
			# removing direction and the first comma
			# then splitting for the comma to get a list of nodes
			# where that sequence aligned
			path = path.replace("<", ",").replace(">", ",")[1:].split(",")
			for s in path:
				if s in info:
					info[s].append(seq_name)
				else:
					info[s] = [seq_name]
	return info


def output_coloring(info1, color1, info2, color2, prefix):
	nodes = dict()
	for s in info1.keys():
		if s in nodes:
			nodes[s][0] += len(info1[s])
		else:
			nodes[s] = [len(info1[s]), 0]

	for s in info2.keys():
		if s in nodes:
			nodes[s][1] += len(info2[s])
		else:
			nodes[s] = [0, len(info2[s])]
	
	with open(f"{prefix}_graph_coloring.csv", "w") as outfile:
		outfile.write(f"Name,Colour,sample1,sample2,sample1_sample2\n")
		for s, values in nodes.items():
			if values[0]/sum(values) >= 0.9:
				color = color1
			elif values[1]/sum(values) >= 0.9:
				color = color2
			else:
				color = "gray"

			outfile.write(f"{s},{color},{values[0]},{values[1]},sample1_{values[0]}_samle2_{values[1]}\n")


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("You need to give the <gaf1> <color> <gaf2> <color> <preifx>")
		print("The color needs to be an accepted color by X11")
		sys.exit()

	gaf1 = sys.argv[1]
	color1 = sys.argv[2]
	gaf2 = sys.argv[3]
	color2 = sys.argv[4]
	prefix = sys.argv[5]

	info1 = process_gaf(gaf1)
	info2 = process_gaf(gaf2)
	output_coloring(info1, color1, info2, color2, prefix)
