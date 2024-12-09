import sys
import pysam
import pdb
import logging


chr_color = {
"chr1":"#000000",
"chr2":"#0000FF",
"chr3":"#A52A2A",
"chr4":"#7FFF00",
"chr5":"#FF7F50",
"chr6":"#B8860B",
"chr7":"#A9A9A9",
"chr8":"#8B008B",
"chr9":"#DC143C",
"chr10":"#8FBC8F",
"chr11":"#FF1493",
"chr12":"#FFD700",
"chr13":"#FF69B4",
"chr14":"#4B0082",
"chr15":"#F0E68C",
"chr16":"#BC8F8F",
"chr17":"#D2B48C",
"chr18":"#FFFF00",
"chr19":"#9ACD32",
"chr20":"#006400",
"chr21":"#DDA0DD",
"chr22":"#98FB98",
"chrX":"#FFEBCD",
"chrY":"#808080",
"chrM":"#20B2AA"}


class Alignment:
	def __init__(self, sam_alignment):
		self.flag = sam_alignment.flag
		self.a_len = sam_alignment.alen
		self.mapq = sam_alignment.mapq
		self.q_name = sam_alignment.qname
		# infer_query_length does not include hard-clipped bases
		self.q_len = sam_alignment.infer_query_length()
		self.q_start = sam_alignment.qstart
		self.q_end = sam_alignment.qend
		self.ref_start = sam_alignment.reference_start
		self.ref_end = sam_alignment.reference_end
		self.ref_name = sam_alignment.reference_name
		if sam_alignment.is_reverse:
			self.orientation = "-"
		else:
			self.orientation = "+"

	# def __eq__(self, other):
	# 	return self.q_name == other.q_name


	def __hash__(self):
		return hash(self.q_name)


	def compare(self, other):
		if not isinstance(other, Alignment):
			print("You need to give the compare function an Alignment object")
			sys.exit()

		# first check if one is not primary alignment
		# if (self.flag in {0, 16}) and (other.flag not in {0, 16}):
		# 	return self
		# elif (self.flag not in {0, 16}) and (other.flag in {0, 16}):
		# 	return other

		if self.a_len > other.a_len:
			return self
		else:
			return other


def sys_exit():
	print("Error Happened! Please check the log file")
	sys.exit()


def filter_and_output(in_bam, out_table, tsv=False):
	if tsv:
		sep = "\t"
	else:
		sep = ","

	try:
		if in_bam.endswith(".bam"):
			alignment_file = pysam.AlignmentFile(in_bam, "rb")
		else:
			alignment_file = pysam.AlignmentFile(in_bam, "r")
	except FileNotFoundError:
		logging.error(f"The file {in_bam} does not exist")
		sys_exit(1)
	except ValueError:
		logging.error(f"The file {in_bam} is not a bam file")
		sys_exit(1)

	"""
	I need to change this to be a bit more smart, also takes into account supplementary alignments
	Or at least have better comparison and not just look at flags 0 and 16
	"""

	"""
	If the output is TSV, I want to output the original size of the node, mapping quality, alignment length
	These information is then used when drawing the SVG plot
	"""
	all_alignments = dict()
	for a in alignment_file:
		if a.flag in {0, 16}:
			if a.qname in all_alignments:
				all_alignments[a.qname] = compare(all_alignments[a.qname], a)
			else:
				all_alignments[a.qname] = Alignment(a)

	with open(out_table, "w") as out_file:
		if tsv:
			out_file.write("node\tchromosome\tstart\tend\torientation\tmapq\n")
		else:
			out_file.write("Name,Colour,Chromosome,al_len,al_start,al_end\n")
		for a in all_alignments.values():
			n_id = a.q_name.split("_")[0]  # node id
			if tsv:
				out_file.write(f"{n_id}{sep}{a.ref_name}{sep}{a.ref_start}{sep}{a.ref_end}{sep}{a.orientation}{sep}{a.mapq}\n")
			else:
				out_file.write(f"{n_id}{sep}{chr_color[a.ref_name]}{sep}{a.ref_name}{sep}{a.a_len}{sep}{a.ref_start}{sep}{a.ref_end}\n")


###################################################################
if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Please give the input sam or bam file and output csv file name")
		sys.exit()

	in_alignment = sys.argv[1]
	out_csv = sys.argv[2]

	filter_and_output(in_alignment, out_csv)
