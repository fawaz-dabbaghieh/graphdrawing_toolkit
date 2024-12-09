import os
import sys
import pysam
from graphdraw.Graph import Graph
from graphdraw.graph_io import *
import pdb


def extract_subgraph(in_gfa, in_bam, chrom, start, end, n_size, output_gfa):
	try:
		bamfile = pysam.AlignmentFile(in_bam, "rb")
	except FileNotFoundError:
		logging.error(f"The file {in_bam} does not exist")
		sys.exit()
	except ValueError:
		logging.error(f"The file {in_bam} is not a bam file")
		sys.exit()

	# reading graph
	graph = Graph(in_gfa)

	bam_iter = bamfile.fetch(chrom, start, end)
	nodes = set()

	outcsv = output_gfa.split(".")[0]
	outcsv += ".csv"
	open_csv = open(outcsv, "w")
	open_csv.write("Name,Colour,Chrom,start_end,align_length\n")
	out_nodes = set()
	for x in bam_iter:
		out_nodes.add(x.qname)
		new_n = x.qname.split("_")[0]
		ref = x.reference_name
		ref_start = x.reference_start
		ref_end = x.reference_end
		open_csv.write(f"{new_n},{'Red'},{ref},{str(ref_start)+'_'+str(ref_end)},{ref_end - ref_start}\n")
		# print(x.qname)
		nodes.add(new_n)
		for n in graph.bfs(new_n, n_size):
			nodes.add(n)
	# pdb.set_trace()
	write_gfa(graph, nodes, output_gfa)


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("you need to give the <gfa> <bam> <chr:start-end> <neighborhood_size> <output.gfa>")
		sys.exit()

	in_gfa = sys.argv[1]
	in_bam = sys.argv[2]
	chromosome, coord = sys.argv[3].split(":")
	start = int(coord.split("-")[0].replace(",", ""))
	end = int(coord.split("-")[1].replace(",", ""))
	n_size = int(sys.argv[4])
	output_gfa = sys.argv[5]

	extract_subgraph(in_gfa, in_bam, chromosome, start, end, n_size, output_gfa)
